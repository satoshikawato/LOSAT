//! Compressed amino-acid lookup primitives for the future `blastp-fast` path.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//! Reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h

use super::{ilog2, PV_ARRAY_BTS, PV_ARRAY_MASK};
use crate::algorithm::tblastx::scan::OffsetPair;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::build_matrix_info;
use crate::stats::{compute_std_aa_composition, lookup_protein_params_ungapped};
use crate::utils::matrix::{aa_char_to_ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

// EXPERIMENT (LOSAT_X_BPCLUT; BPCLUT is a strict-set switch): this file's table built flat
// (count, layout, fill), each cell's hits in one slice in NCBI's scan order; a child module, so that
// it uses this file's private builder pieces (x_clut_direct.rs, exactness argument in its header).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
#[path = "x_clut_direct.rs"]
pub mod x_clut_direct;

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:186-201
// ```c
// #define COMPRESSED_HITS_PER_BACKBONE_CELL 4
// #define COMPRESSED_HITS_CELL_MASK 0x03
// #define COMPRESSED_HITS_PER_OVERFLOW_CELL 4
// #define COMPRESSED_OVERFLOW_CELLS_IN_BANK 209710
// #define COMPRESSED_OVERFLOW_MAX_BANKS 1024
// ```
pub const COMPRESSED_HITS_PER_BACKBONE_CELL: usize = 4;
pub const COMPRESSED_HITS_CELL_MASK: usize = 0x03;
pub const COMPRESSED_HITS_PER_OVERFLOW_CELL: usize = 4;
pub const COMPRESSED_OVERFLOW_CELLS_IN_BANK: usize = 209_710;
pub const COMPRESSED_OVERFLOW_MAX_BANKS: usize = 1024;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1303-1306
// ```c
// /* compute a custom score matrix, for use only
//    with lookup table creation. The matrix dimensions
//    are BLASTAA_SIZE x compressed_alphabet_size, and
//    the score entries are scaled up by kMatrixScale */
// ```
pub type CompressedScoreMatrix = [[i32; BLASTAA_SIZE]; BLASTAA_SIZE];

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:203-208
// ```c
// typedef struct CompressedOverflowCell {
//     struct CompressedOverflowCell* next;
//     Int4 query_offsets[COMPRESSED_HITS_PER_OVERFLOW_CELL];
// } CompressedOverflowCell;
// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedOverflowCell {
    pub next: Option<usize>,
    pub query_offsets: [i32; COMPRESSED_HITS_PER_OVERFLOW_CELL],
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:763-775
// ```c
// lookup->overflow_banks[bank_idx] = (CompressedOverflowCell*) malloc(...);
// return lookup->overflow_banks[lookup->curr_overflow_bank] +
//                                lookup->curr_overflow_cell++;
// ```
impl Default for CompressedOverflowCell {
    fn default() -> Self {
        Self {
            next: None,
            query_offsets: [0; COMPRESSED_HITS_PER_OVERFLOW_CELL],
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:210-231
// ```c
// typedef struct CompressedMixedOffsets{
//     Int4 query_offsets[COMPRESSED_HITS_PER_BACKBONE_CELL-2];
//     CompressedOverflowCell* head;
// } CompressedMixedOffsets;
// typedef struct CompressedLookupBackboneCell {
//     Int4 num_used;
//     Int4 query_offset;
//     union { Int4 query_offsets[4]; CompressedMixedOffsets overflow_list; } payload;
// } CompressedLookupBackboneCell;
// ```
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CompressedBackbonePayload {
    Inline([i32; COMPRESSED_HITS_PER_BACKBONE_CELL]),
    Overflow {
        query_offsets: [i32; COMPRESSED_HITS_PER_BACKBONE_CELL - 2],
        head: Option<usize>,
    },
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:787-801
// ```c
// CompressedLookupBackboneCell *backbone_cell = lookup->backbone + index;
// Int4 num_entries = backbone_cell->num_used;
// switch (num_entries) { case 0: ... case 1: ... }
// ```
impl Default for CompressedBackbonePayload {
    fn default() -> Self {
        Self::Inline([0; COMPRESSED_HITS_PER_BACKBONE_CELL])
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:221-231
// ```c
// typedef struct CompressedLookupBackboneCell {
//     Int4 num_used;
//     Int4 query_offset;
//     union { ... } payload;
// } CompressedLookupBackboneCell;
// ```
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompressedLookupBackboneCell {
    pub num_used: i32,
    pub query_offset: i32,
    pub payload: CompressedBackbonePayload,
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:236-267
// ```c
// typedef struct BlastCompressedAaLookupTable {
//     Int4 threshold;
//     Int4 word_length;
//     Int4 alphabet_size;
//     Int4 compressed_alphabet_size;
//     Int4 reciprocal_alphabet_size;
//     ...
//     PV_ARRAY_TYPE *pv;
//     Int4 pv_array_bts;
//     Uint1* compress_table;
//     Int4* scaled_compress_table;
// } BlastCompressedAaLookupTable;
// ```
#[derive(Clone, Debug)]
pub struct BlastCompressedAaLookupTable {
    pub threshold: i32,
    pub word_length: i32,
    pub alphabet_size: i32,
    pub compressed_alphabet_size: i32,
    pub reciprocal_alphabet_size: i32,
    pub longest_chain: i32,
    pub backbone_size: i32,
    pub backbone: Vec<CompressedLookupBackboneCell>,
    pub overflow_cells: Vec<CompressedOverflowCell>,
    pub pv: Vec<u32>,
    pub pv_array_bts: i32,
    pub compress_table: [u8; BLASTAA_SIZE],
    pub scaled_compress_table: [i32; BLASTAA_SIZE],
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:270-278
    // ```c
    // scaled_compress_table = lookup->scaled_compress_table;
    // ...
    // compressed_char = scaled_compress_table[next_char];
    // next_char = s[word_length];
    // ```
    //
    // Rust keeps a byte-wide scan mirror so the hot compressed BLASTP scan can
    // index by `Uint1`-equivalent subject bytes without an Option/bounds branch;
    // entries outside NCBISTDAA remain `-1`, preserving the bad-character flush.
    pub(crate) scaled_compress_table_u8: [i32; 256],
    pub neighbor_matches: i32,
    pub exact_matches: i32,
    /// Whether the words filled NCBI's overflow banks: NCBI then writes past
    /// `overflow_banks` and corrupts its heap, so the table is not used.
    pub overflow_banks_exhausted: bool,
    /// The overflow cells the table can take: NCBI's banks (`NCBI_OVERFLOW_CELL_CAPACITY`),
    /// or no limit for BLASTP's one table of all query batches
    /// (`try_build_blosum62_compressed_lookup`), which tests each batch's own table instead.
    overflow_cell_capacity: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1321-1330
// ```c
//     lookup->overflow_banks = (CompressedOverflowCell **) calloc(
//                                          COMPRESSED_OVERFLOW_MAX_BANKS,
//                                          sizeof(CompressedOverflowCell *));
//     ...
//     lookup->curr_overflow_cell = COMPRESSED_OVERFLOW_CELLS_IN_BANK;
//     lookup->curr_overflow_bank = -1;
// ```
// The cells a table can take before `s_CompressedListGetNewCell` writes past `overflow_banks`.
const NCBI_OVERFLOW_CELL_CAPACITY: usize =
    COMPRESSED_OVERFLOW_CELLS_IN_BANK * COMPRESSED_OVERFLOW_MAX_BANKS;

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:935-952
// ```c
// typedef struct CompressedNeighborInfo {
//     Int4 row_max[BLASTAA_SIZE];
//     Int4 matrixSorted[BLASTAA_SIZE][BLASTAA_SIZE];
//     Uint1 matrixSortedChar[BLASTAA_SIZE][BLASTAA_SIZE];
// } CompressedNeighborInfo;
// ```
#[derive(Clone, Debug)]
pub struct CompressedNeighborInfo {
    row_max: [i32; BLASTAA_SIZE],
    matrix_sorted: CompressedScoreMatrix,
    matrix_sorted_char: [[u8; BLASTAA_SIZE]; BLASTAA_SIZE],
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1289-1320
// ```c
// lookup->word_length = word_size;
// lookup->threshold = (Int4)(kMatrixScale * opt->threshold);
// if (word_size == 6 || word_size == 5) {
//     lookup->compressed_alphabet_size = 15;
//     lookup->reciprocal_alphabet_size = 286331154;
// }
// else {
//     lookup->compressed_alphabet_size = 10;
//     lookup->reciprocal_alphabet_size = 429496730;
// }
// lookup->backbone_size = (Int4)pow(lookup->compressed_alphabet_size,
//                                   word_size) + 1;
// lookup->backbone = calloc(lookup->backbone_size, sizeof(...));
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1292
// ```c
//     lookup->threshold = (Int4)(kMatrixScale * opt->threshold);
// ```
/// The threshold of a compressed table in NCBI's matrix units, for a caller of the builders that
/// take it converted (BLASTX, `blastx/seed.rs`).
pub(crate) fn x_compressed_scaled_threshold(threshold: f64) -> i32 {
    // `(Int4)` as the x86-64 build converts it (`INT_MIN` past `Int4`).
    crate::core::blast_util::ncbi_int4_from_double(100.0 * threshold)
}

impl BlastCompressedAaLookupTable {
    pub fn new(word_size: usize, threshold: f64) -> Option<Self> {
        let (compressed_alphabet_size, reciprocal_alphabet_size) =
            compressed_lookup_parameters(word_size)?;
        let compress_table = build_compressed_translation(compressed_alphabet_size)?;
        let scaled_compress_table =
            build_scaled_compress_table(&compress_table, compressed_alphabet_size, word_size);
        let scaled_compress_table_u8 = build_scaled_compress_table_u8(&scaled_compress_table);
        let backbone_size = compressed_alphabet_size
            .checked_pow(word_size as u32)?
            .checked_add(1)?;

        Some(Self {
            threshold: x_compressed_scaled_threshold(threshold),
            word_length: word_size as i32,
            alphabet_size: BLASTAA_SIZE as i32,
            compressed_alphabet_size: compressed_alphabet_size as i32,
            reciprocal_alphabet_size,
            longest_chain: 0,
            backbone_size: backbone_size as i32,
            backbone: vec![CompressedLookupBackboneCell::default(); backbone_size],
            overflow_cells: Vec::new(),
            pv: Vec::new(),
            pv_array_bts: PV_ARRAY_BTS as i32,
            compress_table,
            scaled_compress_table,
            scaled_compress_table_u8,
            neighbor_matches: 0,
            exact_matches: 0,
            overflow_banks_exhausted: false,
            overflow_cell_capacity: NCBI_OVERFLOW_CELL_CAPACITY,
        })
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:753-775
    // ```c
    // if (lookup->curr_overflow_cell == COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
    //     lookup->overflow_banks[bank_idx] = malloc(...);
    //     lookup->curr_overflow_bank++;
    //     lookup->curr_overflow_cell = 0;
    // }
    // return lookup->overflow_banks[lookup->curr_overflow_bank] +
    //                            lookup->curr_overflow_cell++;
    // ```
    // NCBI reference: c++/src/algo/blast/core/blast_aalookup.c:759-768
    // ```c
    //     if (lookup->curr_overflow_cell ==
    //                         COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
    //         /* need a new bank */
    //         Int4 bank_idx = lookup->curr_overflow_bank + 1;
    //         lookup->overflow_banks[bank_idx] = (CompressedOverflowCell*) malloc(
    //                                            COMPRESSED_OVERFLOW_CELLS_IN_BANK *
    //                                            sizeof(CompressedOverflowCell));
    //         ASSERT(bank_idx < COMPRESSED_OVERFLOW_MAX_BANKS);
    //         ASSERT(lookup->overflow_banks[bank_idx]);
    //         lookup->curr_overflow_bank++;
    // ```
    // The release build has no ASSERT, so after the last bank NCBI writes past
    // `overflow_banks`; the table is marked and the last cell is reused, and the words that
    // are left are not added. The capacity is that of NCBI's banks, except for BLASTP's
    // table of all query batches, which has none (`overflow_cell_capacity`); the cells are
    // taken in the same order either way.
    fn compressed_list_get_new_cell(&mut self) -> usize {
        let cell_index = self.overflow_cells.len();
        let max_cells = self.overflow_cell_capacity;
        if cell_index >= max_cells {
            self.overflow_banks_exhausted = true;
            return max_cells - 1;
        }
        if cell_index == self.overflow_cells.capacity() {
            // No NCBI counterpart: a table with NCBI's capacity grows by one bank, as NCBI
            // allocates its banks; BLASTP's table of all query batches, which has no limit,
            // grows geometrically, so that its cells are not copied once per bank on platforms
            // whose `realloc` copies. Only the capacity differs; the cells are the same.
            let additional = if max_cells > NCBI_OVERFLOW_CELL_CAPACITY {
                cell_index.max(COMPRESSED_OVERFLOW_CELLS_IN_BANK)
            } else {
                COMPRESSED_OVERFLOW_CELLS_IN_BANK.min(max_cells - cell_index)
            };
            self.overflow_cells.reserve_exact(additional);
        }
        self.overflow_cells.push(CompressedOverflowCell::default());
        cell_index
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:782-849
    // ```c
    // switch (num_entries) {
    // case 0:
    //     backbone_cell->query_offset = query_offset;
    //     break;
    // case 1: case 2: case 3: case 4:
    //     backbone_cell->payload.query_offsets[num_entries-1] = query_offset;
    //     break;
    // case 5:
    //     tmp_offset_0 = backbone_cell->payload.query_offsets[0];
    //     tmp_offset_1 = backbone_cell->payload.query_offsets[1];
    //     new_cell = s_CompressedListGetNewCell(lookup);
    //     new_cell->query_offsets[0] = backbone_cell->payload.query_offsets[2];
    //     new_cell->query_offsets[1] = backbone_cell->payload.query_offsets[3];
    //     new_cell->query_offsets[2] = query_offset;
    //     backbone_cell->payload.overflow_list.query_offsets[0] = tmp_offset_0;
    //     backbone_cell->payload.overflow_list.query_offsets[1] = tmp_offset_1;
    //     backbone_cell->payload.overflow_list.head = new_cell;
    //     break;
    // default:
    //     cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK;
    //     if (cell_index == 0) {
    //         new_cell = s_CompressedListGetNewCell(lookup);
    //         new_cell->next = backbone_cell->payload.overflow_list.head;
    //         backbone_cell->payload.overflow_list.head = new_cell;
    //     }
    //     backbone_cell->payload.overflow_list.head->query_offsets[cell_index] =
    //         query_offset;
    // }
    // backbone_cell->num_used++;
    // ```
    pub fn add_word_hit(&mut self, index: usize, query_offset: i32) {
        let num_entries = self.backbone[index].num_used;
        match num_entries {
            0 => {
                self.backbone[index].query_offset = query_offset;
            }
            1..=4 => {
                let CompressedBackbonePayload::Inline(query_offsets) =
                    &mut self.backbone[index].payload
                else {
                    unreachable!("NCBI BLAST overflow list cannot exist before entry 6");
                };
                query_offsets[(num_entries - 1) as usize] = query_offset;
            }
            5 => {
                let query_offsets = {
                    let CompressedBackbonePayload::Inline(query_offsets) =
                        &self.backbone[index].payload
                    else {
                        unreachable!("NCBI BLAST overflow list cannot exist before entry 6");
                    };
                    *query_offsets
                };
                let new_cell = self.compressed_list_get_new_cell();
                self.overflow_cells[new_cell].next = None;
                self.overflow_cells[new_cell].query_offsets[0] = query_offsets[2];
                self.overflow_cells[new_cell].query_offsets[1] = query_offsets[3];
                self.overflow_cells[new_cell].query_offsets[2] = query_offset;
                self.backbone[index].payload = CompressedBackbonePayload::Overflow {
                    query_offsets: [query_offsets[0], query_offsets[1]],
                    head: Some(new_cell),
                };
            }
            _ => {
                let cell_index = usize::try_from(num_entries - 3)
                    .expect("NCBI BLAST compressed hit count must fit in usize")
                    & COMPRESSED_HITS_CELL_MASK;
                let head = match &self.backbone[index].payload {
                    CompressedBackbonePayload::Overflow { head, .. } => *head,
                    CompressedBackbonePayload::Inline(_) => {
                        unreachable!("NCBI BLAST entry 7+ requires overflow storage")
                    }
                };
                let head = if cell_index == 0 {
                    let new_cell = self.compressed_list_get_new_cell();
                    self.overflow_cells[new_cell].next = head;
                    let CompressedBackbonePayload::Overflow { head, .. } =
                        &mut self.backbone[index].payload
                    else {
                        unreachable!("NCBI BLAST entry 7+ requires overflow storage");
                    };
                    *head = Some(new_cell);
                    new_cell
                } else {
                    head.expect("NCBI BLAST compressed overflow head must exist")
                };
                self.overflow_cells[head].query_offsets[cell_index] = query_offset;
            }
        }
        self.backbone[index].num_used += 1;
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:859-909
    // ```c
    // switch (lookup->word_length) {
    // case 5:
    //     index = w[0] + W6p1[w[1]] + W6p2[w[2]] +
    //             W6p3[w[3]] + W6p4[w[4]];
    //     break;
    // case 6:
    //     index = w[0] + W6p1[w[1]] + W6p2[w[2]] +
    //             W6p3[w[3]] + W6p4[w[4]] + W6p5[w[5]];
    //     break;
    // case 7:
    //     index = w[0] + W7p1[w[1]] + ... + W7p6[w[6]];
    //     break;
    // }
    // s_CompressedLookupAddWordHit(lookup, index, query_offset);
    // ```
    pub fn add_encoded(&mut self, word: &[u8], query_offset: i32) {
        let index = compute_encoded_compressed_index(
            self.word_length as usize,
            word,
            self.compressed_alphabet_size as usize,
        );
        self.add_word_hit(index, query_offset);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:919-929
    // ```c
    // index = s_ComputeCompressedIndex(lookup->word_length, w,
    //                                  lookup->compressed_alphabet_size,
    //                                  &skip, lookup);
    // if (skip == 0)
    //     s_CompressedLookupAddWordHit(lookup, index, query_offset);
    // ```
    pub fn add_unencoded(&mut self, word: &[u8], query_offset: i32) -> bool {
        let (index, skip) = compute_compressed_index(
            self.word_length as usize,
            word,
            self.compressed_alphabet_size as usize,
            &self.scaled_compress_table,
        );
        if skip == 0 {
            self.add_word_hit(index as usize, query_offset);
            true
        } else {
            false
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1071-1110
    // ```c
    // lookup->exact_matches++;
    // score = 0;
    // for (i = 0; i < lookup->word_length; i++) {
    //     int c = lookup->compress_table[w[i]];
    //     if (c >= lookup->compressed_alphabet_size)
    //         return;
    //     score += info->matrix[w[i]][c];
    // }
    // if (lookup->threshold == 0 || score < lookup->threshold) {
    //     s_CompressedLookupAddUnencoded(lookup, w, query_offset);
    // }
    // ...
    // score = info->row_max[w[0]];
    // for (i = 1; i < lookup->word_length; i++)
    //     score += info->row_max[w[i]];
    // s_CompressedAddWordHitsCore(info, score, 0);
    // ```
    pub fn add_word_hits(
        &mut self,
        info: &CompressedNeighborInfo,
        matrix: &CompressedScoreMatrix,
        query: &[u8],
        query_offset: i32,
    ) {
        self.exact_matches += 1;
        let word_size = self.word_length as usize;
        let offset = usize::try_from(query_offset)
            .expect("NCBI BLAST compressed query offset must be non-negative");
        let word = &query[offset..offset + word_size];

        let mut self_score = 0i32;
        for &residue in word {
            let compressed = self.compress_table[residue as usize] as usize;
            if compressed >= self.compressed_alphabet_size as usize {
                return;
            }
            self_score += matrix[residue as usize][compressed];
        }

        if self.threshold == 0 || self_score < self.threshold {
            self.add_unencoded(word, query_offset);
        } else {
            self.neighbor_matches -= 1;
        }

        if self.threshold == 0 {
            return;
        }

        let mut score = info.row_max[word[0] as usize];
        for &residue in &word[1..] {
            score += info.row_max[residue as usize];
        }
        let mut subject_word = [0u8; 32];
        self.add_word_hits_core(info, word, &mut subject_word, score, 0, query_offset);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1008-1053
    // ```c
    // score -= info->row_max[currQueryChar];
    // rowSorted = info->matrixSorted[currQueryChar];
    // charSorted = info->matrixSortedChar[currQueryChar];
    // if (current_pos == info->wordsize - 1) {
    //     for (i = 0; i < info->compressed_alphabet_size &&
    //            (score + rowSorted[i] >= info->threshold); i++) {
    //         subject_word[current_pos] = charSorted[i];
    //         s_CompressedLookupAddEncoded(lookup, subject_word, query_offset);
    //     }
    //     return;
    // }
    // for (i = 0; i < info->compressed_alphabet_size &&
    //        (score + rowSorted[i] >= info->threshold); i++) {
    //     subject_word[current_pos] = charSorted[i];
    //     s_CompressedAddWordHitsCore(info, score + rowSorted[i],
    //                                 current_pos + 1);
    // }
    // ```
    fn add_word_hits_core(
        &mut self,
        info: &CompressedNeighborInfo,
        query_word: &[u8],
        subject_word: &mut [u8; 32],
        mut score: i32,
        current_pos: usize,
        query_offset: i32,
    ) {
        if self.overflow_banks_exhausted {
            return;
        }
        let curr_query_char = query_word[current_pos] as usize;
        score -= info.row_max[curr_query_char];
        let row_sorted = &info.matrix_sorted[curr_query_char];
        let char_sorted = &info.matrix_sorted_char[curr_query_char];
        let compressed_alphabet_size = self.compressed_alphabet_size as usize;
        let word_size = self.word_length as usize;

        if current_pos == word_size - 1 {
            for i in 0..compressed_alphabet_size {
                if score + row_sorted[i] < self.threshold {
                    break;
                }
                subject_word[current_pos] = char_sorted[i];
                self.add_encoded(&subject_word[..word_size], query_offset);
                self.neighbor_matches += 1;
            }
            return;
        }

        for i in 0..compressed_alphabet_size {
            if score + row_sorted[i] < self.threshold {
                break;
            }
            subject_word[current_pos] = char_sorted[i];
            self.add_word_hits_core(
                info,
                query_word,
                subject_word,
                score + row_sorted[i],
                current_pos + 1,
                query_offset,
            );
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1137-1176
    // ```c
    // for (i = 0; i < lookup->alphabet_size; i++) {
    //     info.row_max[i] = compressed_matrix[i][0];
    //     for (j = 1; j < lookup->compressed_alphabet_size; j++)
    //         info.row_max[i] = MAX(info.row_max[i], compressed_matrix[i][j]);
    // }
    // s_loadSortedMatrix(&info);
    // for (loc = location; loc; loc = loc->next){
    //     Int4 from = loc->ssr->left;
    //     Int4 to = loc->ssr->right - lookup->word_length + 1;
    //     for (offset = from; offset <= to; offset++){
    //         s_CompressedAddWordHits(&info, query->sequence, offset);
    //     }
    // }
    // ```
    pub fn add_neighboring_words(
        &mut self,
        matrix: &CompressedScoreMatrix,
        query: &[u8],
        locations: &[(i32, i32)],
    ) {
        let info = CompressedNeighborInfo::new(matrix, self.compressed_alphabet_size as usize);
        self.add_location_words(&info, matrix, query, locations, |_| true);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1168-1175
    // ```c
    //     for (loc = location; loc; loc = loc->next){
    //         Int4 from = loc->ssr->left;
    //         Int4 to = loc->ssr->right - lookup->word_length + 1;
    //
    //         for (offset = from; offset <= to; offset++){
    //             s_CompressedAddWordHits(&info, query->sequence, offset);
    //         }
    //     }
    // ```
    // The loop of `add_neighboring_words`. No NCBI counterpart: `after_word` sees the table
    // after the words of each query offset (BLASTP's test of its query batches' own tables)
    // and stops the loop when it returns `false`; it changes no word that is added.
    fn add_location_words(
        &mut self,
        info: &CompressedNeighborInfo,
        matrix: &CompressedScoreMatrix,
        query: &[u8],
        locations: &[(i32, i32)],
        mut after_word: impl FnMut(&Self) -> bool,
    ) -> bool {
        let word_size = self.word_length;
        for &(left, right) in locations {
            let to = right - word_size + 1;
            for offset in left..=to {
                let start = usize::try_from(offset)
                    .expect("NCBI BLAST compressed query offset must be non-negative");
                let end = start + word_size as usize;
                if end <= query.len() {
                    self.add_word_hits(info, matrix, query, offset);
                    if !after_word(self) {
                        return false;
                    }
                }
            }
        }
        true
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:343-384
    // ```c
    // dest[0].qs_offsets.q_off = backbone_cell->query_offset;
    // if (numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL+1) {
    //     for (i = 0; i < numhits-1; i++)
    //         dest[i].qs_offsets.q_off = backbone_cell->payload.query_offsets[i];
    // } else {
    //     first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
    //     dest[1].qs_offsets.q_off =
    //          backbone_cell->payload.overflow_list.query_offsets[0];
    //     dest[2].qs_offsets.q_off =
    //          backbone_cell->payload.overflow_list.query_offsets[1];
    //     for (i = 0; i < first_cell_entries; i++)
    //         dest[i].qs_offsets.q_off = curr_cell->query_offsets[i];
    //     while (curr_cell != NULL) { ... }
    // }
    // ```
    pub fn query_offsets_in_scan_order(&self, index: usize) -> Vec<i32> {
        let cell = &self.backbone[index];
        let numhits = cell.num_used;
        if numhits <= 0 {
            return Vec::new();
        }

        let mut offsets =
            Vec::with_capacity(usize::try_from(numhits).expect("NCBI hit count fits usize"));
        offsets.push(cell.query_offset);
        if numhits <= (COMPRESSED_HITS_PER_BACKBONE_CELL as i32 + 1) {
            let CompressedBackbonePayload::Inline(query_offsets) = &cell.payload else {
                unreachable!("NCBI BLAST small compressed cell uses inline storage");
            };
            for i in 0..(numhits - 1) as usize {
                offsets.push(query_offsets[i]);
            }
            return offsets;
        }

        let CompressedBackbonePayload::Overflow {
            query_offsets,
            head,
        } = &cell.payload
        else {
            unreachable!("NCBI BLAST large compressed cell uses overflow storage");
        };
        offsets.push(query_offsets[0]);
        offsets.push(query_offsets[1]);

        let first_cell_entries = usize::try_from(numhits - 3)
            .expect("NCBI compressed overflow hit count must fit usize")
            & COMPRESSED_HITS_CELL_MASK;
        let mut curr_cell = head.expect("NCBI BLAST compressed overflow head must exist");
        for i in 0..first_cell_entries {
            offsets.push(self.overflow_cells[curr_cell].query_offsets[i]);
        }
        if first_cell_entries != 0 {
            if let Some(next) = self.overflow_cells[curr_cell].next {
                curr_cell = next;
            } else {
                return offsets;
            }
        }
        loop {
            for i in 0..COMPRESSED_HITS_PER_OVERFLOW_CELL {
                offsets.push(self.overflow_cells[curr_cell].query_offsets[i]);
            }
            let Some(next) = self.overflow_cells[curr_cell].next else {
                break;
            };
            curr_cell = next;
        }
        offsets.truncate(numhits as usize);
        offsets
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:343-384
    // ```c
    // dest[0].qs_offsets.q_off = backbone_cell->query_offset;
    // dest[0].qs_offsets.s_off = s_off;
    // if (numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL+1) {
    //     for (i = 0; i < numhits-1; i++) {
    //         dest[i].qs_offsets.q_off = backbone_cell->payload.query_offsets[i];
    //         dest[i].qs_offsets.s_off = s_off;
    //     }
    // } else { ... while (curr_cell != NULL) { ... } }
    // ```
    #[allow(dead_code)]
    fn copy_scan_hits(&self, index: usize, s_off: u32, dest: &mut [OffsetPair]) {
        let cell = &self.backbone[index];
        let numhits =
            usize::try_from(cell.num_used).expect("NCBI BLAST compressed hit count must fit usize");
        if numhits == 0 {
            return;
        }

        dest[0] = OffsetPair {
            q_off: cell.query_offset as u32,
            s_off,
        };
        if numhits <= COMPRESSED_HITS_PER_BACKBONE_CELL + 1 {
            let CompressedBackbonePayload::Inline(query_offsets) = &cell.payload else {
                unreachable!("NCBI BLAST small compressed cell uses inline storage");
            };
            for i in 0..(numhits - 1) {
                dest[i + 1] = OffsetPair {
                    q_off: query_offsets[i] as u32,
                    s_off,
                };
            }
            return;
        }

        let CompressedBackbonePayload::Overflow {
            query_offsets,
            head,
        } = &cell.payload
        else {
            unreachable!("NCBI BLAST large compressed cell uses overflow storage");
        };
        dest[1] = OffsetPair {
            q_off: query_offsets[0] as u32,
            s_off,
        };
        dest[2] = OffsetPair {
            q_off: query_offsets[1] as u32,
            s_off,
        };

        let first_cell_entries = (numhits - 3) & COMPRESSED_HITS_CELL_MASK;
        let mut written = 3usize;
        let mut curr_cell = head.expect("NCBI BLAST compressed overflow head must exist");
        for i in 0..first_cell_entries {
            dest[written] = OffsetPair {
                q_off: self.overflow_cells[curr_cell].query_offsets[i] as u32,
                s_off,
            };
            written += 1;
        }
        if first_cell_entries != 0 {
            let Some(next) = self.overflow_cells[curr_cell].next else {
                return;
            };
            curr_cell = next;
        }
        while written < numhits {
            for i in 0..COMPRESSED_HITS_PER_OVERFLOW_CELL {
                if written >= numhits {
                    break;
                }
                dest[written] = OffsetPair {
                    q_off: self.overflow_cells[curr_cell].query_offsets[i] as u32,
                    s_off,
                };
                written += 1;
            }
            let Some(next) = self.overflow_cells[curr_cell].next else {
                break;
            };
            curr_cell = next;
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:229-384
    // ```c
    // while (s_DetermineScanningOffsets(subject, word_length, word_length, s_range)) {
    //     for(s = s_first; s <= s_last; s++){
    //         index = s_ComputeCompressedIndex(word_length - 1, s, ...);
    //         if(!skip) break;
    //     }
    //     next_char = ((s <= s_last)? s[word_length-1] : 0);
    //     preshift = (Int4)((((Int8)index) * recip) >> 32);
    //     for (; s <= s_last; s++) {
    //         compressed_char = scaled_compress_table[next_char];
    //         next_char = s[word_length];
    //         ...
    //         if (PV_TEST(pv, index, pv_array_bts)) { ... }
    //     }
    // }
    // ```
    #[allow(dead_code)]
    pub(crate) fn scan_subject(
        &self,
        subject: &[u8],
        seq_ranges: &[(i32, i32)],
        offset_pairs: &mut [OffsetPair],
        array_size: i32,
        s_range: &mut [i32; 3],
    ) -> i32 {
        let mut totalhits = 0i32;
        let word_length = usize::try_from(self.word_length)
            .expect("NCBI BLAST compressed word length must fit usize");
        let compressed_alphabet_size = usize::try_from(self.compressed_alphabet_size)
            .expect("NCBI BLAST compressed alphabet size must fit usize");
        let recip = i64::from(self.reciprocal_alphabet_size);
        let pv_array_bts = usize::try_from(self.pv_array_bts)
            .expect("NCBI BLAST compressed PV shift must fit usize");

        while determine_compressed_scanning_offsets(
            seq_ranges,
            self.word_length,
            self.word_length,
            s_range,
        ) {
            let s_first = usize::try_from(s_range[1])
                .expect("NCBI BLAST compressed scan start must be non-negative");
            let s_last = usize::try_from(s_range[2])
                .expect("NCBI BLAST compressed scan end must be non-negative");
            let mut s = s_first;
            let mut index = 0i32;

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:496-500
            // ```c
            //     scan_range[1] = subject->seq_ranges[0].left;
            //     scan_range[2] = subject->seq_ranges[0].right - wordsize;
            //
            //     if (scan_range[2] < scan_range[1])
            //         scan_range[2] = scan_range[1];
            // ```
            // A subject shorter than the word is still primed once, and the priming word
            // reads past the subject into its closing NULLB and beyond. The NULLB has no
            // compressed letter (its scaled_compress_table entry is -1, so skip is set), so
            // the scan finds nothing whatever follows; letters past the subject are read
            // here as NULLB.
            let mut padded_word = [0u8; 8];
            while s <= s_last {
                let end = s + word_length - 1;
                let word = match subject.get(s..end) {
                    Some(word) => word,
                    None => {
                        let available = subject.get(s..).unwrap_or(&[]);
                        let word = &mut padded_word[..word_length - 1];
                        word.fill(0);
                        word[..available.len()].copy_from_slice(available);
                        &padded_word[..word_length - 1]
                    }
                };
                let computed = compute_compressed_index_u8(
                    word_length - 1,
                    word,
                    compressed_alphabet_size,
                    &self.scaled_compress_table_u8,
                );
                index = computed.0;
                let skip = computed.1;
                if skip == 0 {
                    break;
                }
                s += 1;
            }

            let mut next_char = if s <= s_last {
                subject.get(s + word_length - 1).copied().unwrap_or(0)
            } else {
                0
            };
            let mut preshift = compressed_reciprocal_divide(index, recip);
            let scaled_compress_table = &self.scaled_compress_table_u8;

            while s <= s_last {
                let mut compressed_char = scaled_compress_table[usize::from(next_char)];
                next_char = subject.get(s + word_length).copied().unwrap_or(0);

                if compressed_char < 0 {
                    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:278-300
                    // ```c
                    // if(compressed_char < 0){ /* flush (rare) "bad" character(s) */
                    //     preshift = 0;
                    //     s++;
                    //     for(skip = word_length-1; skip && (s <= s_last) ; s++){ ... }
                    //     s--; /*undo the following increment*/
                    //     continue;
                    // }
                    // ```
                    preshift = 0;
                    s += 1;
                    let mut skip = self.word_length - 1;
                    while skip != 0 && s <= s_last {
                        compressed_char = scaled_compress_table[usize::from(next_char)];
                        next_char = subject.get(s + word_length).copied().unwrap_or(0);
                        if compressed_char < 0 {
                            skip = self.word_length - 1;
                            preshift = 0;
                            s += 1;
                            continue;
                        }
                        index = preshift + compressed_char;
                        preshift = compressed_reciprocal_divide(index, recip);
                        skip -= 1;
                        s += 1;
                    }
                    // NCBI's `s--` is immediately balanced by the outer `for`
                    // increment; this `while` loop already holds that resumed
                    // subject offset.
                    continue;
                }

                index = preshift + compressed_char;
                preshift = compressed_reciprocal_divide(index, recip);

                let index_usize = usize::try_from(index)
                    .expect("NCBI BLAST compressed scan index must be non-negative");
                if compressed_pv_test(&self.pv, index_usize, pv_array_bts) {
                    let numhits = self.backbone[index_usize].num_used;
                    if numhits != 0 {
                        if numhits <= array_size - totalhits {
                            let dest_base = usize::try_from(totalhits)
                                .expect("NCBI BLAST compressed total hits fit usize");
                            let dest_end = dest_base
                                + usize::try_from(numhits)
                                    .expect("NCBI BLAST compressed hit count fits usize");
                            self.copy_scan_hits(
                                index_usize,
                                s as u32,
                                &mut offset_pairs[dest_base..dest_end],
                            );
                            totalhits += numhits;
                        } else {
                            // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:373-378
                            // ```c
                            // s_range[1] = (Int4)(s - subject->sequence);
                            // return totalhits;
                            // ```
                            s_range[1] = i32::try_from(s)
                                .expect("NCBI BLAST compressed subject offset fits Int4");
                            return totalhits;
                        }
                    }
                }
                s += 1;
            }

            s_range[1] = i32::try_from(s).expect("NCBI BLAST compressed subject offset fits Int4");
        }

        totalhits
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1182-1242
    // ```c
    // for (i = count = 0; i < lookup->backbone_size; i++) {
    //     if (lookup->backbone[i].num_used)
    //         count++;
    // }
    // pv_array_bts = PV_ARRAY_BTS;
    // if (count <= 0.01 * lookup->backbone_size) {
    //     pv_array_bts += ilog2(lookup->backbone_size / (8 * kTargetPVBytes));
    // }
    // pv = lookup->pv = calloc((lookup->backbone_size >> pv_array_bts) + 1,
    //                          sizeof(PV_ARRAY_TYPE));
    // for (i = 0; i < lookup->backbone_size; i++) {
    //     count = lookup->backbone[i].num_used;
    //     if (count > 0) {
    //         PV_SET(pv, i, pv_array_bts);
    //         longest_chain = MAX(count, longest_chain);
    //     }
    // }
    // lookup->longest_chain = longest_chain;
    // ```
    pub fn finalize(&mut self) {
        let nonempty_count = self
            .backbone
            .iter()
            .filter(|cell| cell.num_used > 0)
            .count();
        let mut pv_array_bts = PV_ARRAY_BTS;
        let backbone_size = self.backbone.len();
        const TARGET_PV_BYTES: usize = 262_144;
        if (nonempty_count as f64) <= 0.01 * (backbone_size as f64) {
            pv_array_bts += ilog2(backbone_size / (8 * TARGET_PV_BYTES));
        }

        self.pv = vec![0u32; (backbone_size >> pv_array_bts) + 1];
        self.pv_array_bts = pv_array_bts as i32;
        let mut longest_chain = 0i32;
        for (index, cell) in self.backbone.iter().enumerate() {
            let count = cell.num_used;
            if count > 0 {
                compressed_pv_set(&mut self.pv, index, pv_array_bts);
                longest_chain = longest_chain.max(count);
            }
        }
        self.longest_chain = longest_chain;
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/masksubj.inl:43-58
// ```c
// while (range[1] > range[2]) {
//     range[0]++;
//     if (range[0] >= subject->num_seq_ranges)
//         return FALSE;
//     range[1] = subject->seq_ranges[range[0]].left +
//                word_length - lut_word_length;
//     range[2] = subject->seq_ranges[range[0]].right - lut_word_length;
// }
// return TRUE;
// ```
#[inline(always)]
#[allow(dead_code)]
fn determine_compressed_scanning_offsets(
    seq_ranges: &[(i32, i32)],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    while range[1] > range[2] {
        range[0] += 1;
        if range[0] >= seq_ranges.len() as i32 {
            return false;
        }
        let (left, right) = seq_ranges[range[0] as usize];
        range[1] = left + word_length - lut_word_length;
        range[2] = right - lut_word_length;
    }
    true
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:272-274
// ```c
// preshift = (Int4)((((Int8)index) * recip) >> 32);
// ```
#[inline(always)]
#[allow(dead_code)]
pub(crate) fn compressed_reciprocal_divide(index: i32, reciprocal_alphabet_size: i64) -> i32 {
    ((i64::from(index) * reciprocal_alphabet_size) >> 32) as i32
}

impl CompressedNeighborInfo {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:969-998
    // ```c
    // for (longChar = 0; longChar < BLASTAA_SIZE; longChar++) {
    //     for (shortChar = 0; shortChar < info->compressed_alphabet_size; shortChar++) {
    //         sortTable[shortChar].diff = info->row_max[longChar] -
    //                                     info->matrix[longChar][shortChar];
    //         sortTable[shortChar].letter = shortChar;
    //     }
    //     qsort(sortTable, info->compressed_alphabet_size,
    //           sizeof(LetterAndScoreDifferencePair), ScoreDifferenceSort);
    //     for (i = 0; i < info->compressed_alphabet_size; i++) {
    //         info->matrixSorted[longChar][i] = info->matrix[longChar][letter];
    //         info->matrixSortedChar[longChar][i] = letter;
    //     }
    // }
    // ```
    pub fn new(matrix: &CompressedScoreMatrix, compressed_alphabet_size: usize) -> Self {
        let mut row_max = [0i32; BLASTAA_SIZE];
        for row in 0..BLASTAA_SIZE {
            row_max[row] = matrix[row][0];
            for col in 1..compressed_alphabet_size {
                row_max[row] = row_max[row].max(matrix[row][col]);
            }
        }

        let mut matrix_sorted = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];
        let mut matrix_sorted_char = [[0u8; BLASTAA_SIZE]; BLASTAA_SIZE];
        for long_char in 0..BLASTAA_SIZE {
            let mut sort_table = [(0i32, 0u8); BLASTAA_SIZE];
            for (short_char, entry) in sort_table
                .iter_mut()
                .enumerate()
                .take(compressed_alphabet_size)
            {
                *entry = (
                    row_max[long_char] - matrix[long_char][short_char],
                    short_char as u8,
                );
            }
            sort_table[..compressed_alphabet_size].sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));
            for (i, &(_, letter)) in sort_table[..compressed_alphabet_size].iter().enumerate() {
                matrix_sorted[long_char][i] = matrix[long_char][letter as usize];
                matrix_sorted_char[long_char][i] = letter;
            }
        }

        Self {
            row_max,
            matrix_sorted,
            matrix_sorted_char,
        }
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1289-1299
// ```c
// if (word_size == 6 || word_size == 5) {
//     lookup->compressed_alphabet_size = 15;
//     lookup->reciprocal_alphabet_size = 286331154;
// }
// else {
//     lookup->compressed_alphabet_size = 10;
//     lookup->reciprocal_alphabet_size = 429496730;
// }
// ```
pub fn compressed_lookup_parameters(word_size: usize) -> Option<(usize, i32)> {
    match word_size {
        5 | 6 => Some((15, 286_331_154)),
        7 => Some((10, 429_496_730)),
        _ => None,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4949-4952
// ```c
// static const char* s_alphabet10 = "IJLMV AST BDENZ KQR G FY P H C W";
// static const char* s_alphabet15 = "ST IJV LM KR EQZ A G BD P N F Y H C W";
// ```
fn compressed_alphabet_string(compressed_alphabet_size: usize) -> Option<&'static [u8]> {
    match compressed_alphabet_size {
        10 => Some(b"IJLMV AST BDENZ KQR G FY P H C W"),
        15 => Some(b"ST IJV LM KR EQZ A G BD P N F Y H C W"),
        _ => None,
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4786-4805
// ```c
// if (isspace(c)) {
//     compressed_letter++;
//     j = 0;
// }
// else if (isalpha(c)) {
//     Int4 aa_letter = AMINOACID_TO_NCBISTDAA[c];
//     table[aa_letter] = compressed_letter;
//     rev_table[compressed_letter][j++] = aa_letter;
//     rev_table[compressed_letter][j] = -1;
// }
// ```
fn build_compressed_reverse_lookup(compressed_alphabet_size: usize) -> Option<Vec<Vec<usize>>> {
    let alphabet = compressed_alphabet_string(compressed_alphabet_size)?;
    let mut reverse = vec![Vec::new(); compressed_alphabet_size];
    let mut compressed_letter = 0usize;

    for &byte in alphabet {
        if byte.is_ascii_whitespace() {
            compressed_letter += 1;
        } else if byte.is_ascii_alphabetic() {
            let aa = aa_char_to_ncbistdaa(byte) as usize;
            if aa < BLASTAA_SIZE && compressed_letter < compressed_alphabet_size {
                reverse[compressed_letter].push(aa);
            }
        }
    }

    Some(reverse)
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4786-4805
// ```c
// for (i = 0; i < BLASTAA_SIZE; i++)
//     table[i] = compressed_alphabet_size;
// ...
// if (isspace(c)) {
//     compressed_letter++;
//     j = 0;
// }
// else if (isalpha(c)) {
//     Int4 aa_letter = AMINOACID_TO_NCBISTDAA[c];
//     table[aa_letter] = compressed_letter;
// }
// ```
pub fn build_compressed_translation(compressed_alphabet_size: usize) -> Option<[u8; BLASTAA_SIZE]> {
    let alphabet = compressed_alphabet_string(compressed_alphabet_size)?;
    let mut table = [compressed_alphabet_size as u8; BLASTAA_SIZE];
    let mut compressed_letter = 0u8;

    for &byte in alphabet {
        if byte.is_ascii_whitespace() {
            compressed_letter = compressed_letter.saturating_add(1);
        } else if byte.is_ascii_alphabetic() {
            let aa = aa_char_to_ncbistdaa(byte) as usize;
            if aa < BLASTAA_SIZE {
                table[aa] = compressed_letter;
            }
        }
    }

    Some(table)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/ncbi_math.c:437-441
// ```c
// long BLAST_Nint(double x)
// {
//    x += (x >= 0. ? 0.5 : -0.5);
//    return (long)x;
// }
// ```
#[inline]
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4818-4858
// ```c
// Blast_ResFreqStdComp(sbp, rfp);
// for (letter = 0; letter < compressed_alphabet_size; letter++) {
//     double prob_sum = 0.;
//     for (i = 0; i < BLASTAA_SIZE; i++) {
//         Int4 aa = rev_table[letter][i];
//         if (aa < 0) break;
//         prob_sum += rfp->prob[aa];
//     }
//     for (i = 0; i < BLASTAA_SIZE; i++) {
//         Int4 aa = rev_table[letter][i];
//         if (aa < 0) break;
//         compressed_prob[aa] = rfp->prob[aa] / prob_sum;
//     }
// }
// ```
fn compressed_probabilities(compressed_alphabet_size: usize) -> Option<[f64; BLASTAA_SIZE]> {
    let reverse = build_compressed_reverse_lookup(compressed_alphabet_size)?;
    let standard_prob = compute_std_aa_composition();
    let mut compressed_prob = [0.0; BLASTAA_SIZE];
    for group in reverse {
        let prob_sum: f64 = group.iter().map(|&aa| standard_prob[aa]).sum();
        if prob_sum == 0.0 {
            continue;
        }
        for aa in group {
            compressed_prob[aa] = standard_prob[aa] / prob_sum;
        }
    }
    Some(compressed_prob)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_stat.c:4886-4937
// ```c
// lambda = RPSfindUngappedLambda(sbp->name);
// matrix_scale_factor /= lambda;
// std_freqs = _PSIMatrixFrequencyRatiosNew(sbp->name);
// ...
// for (q = 0; q < BLASTAA_SIZE; q++) {
//     for (s = 0; s < compressed_alphabet_size; s++) {
//         double val = 0;
//         for (i = 0; i < BLASTAA_SIZE; i++) {
//             Int4 aa = rev_table[s][i];
//             if (aa < 0) break;
//             val += std_freqs->data[q][aa] * compressed_prob[aa];
//         }
//         val = (val < 1e-8) ? min_freq : log(val);
//         scores[q][s] = (Int4)BLAST_Nint(val * matrix_scale_factor);
//     }
// }
// ```
pub fn build_blosum62_compressed_score_matrix(
    compressed_alphabet_size: usize,
) -> Option<CompressedScoreMatrix> {
    let reverse = build_compressed_reverse_lookup(compressed_alphabet_size)?;
    let compressed_prob = compressed_probabilities(compressed_alphabet_size)?;
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, ungapped.lambda).ok()?;
    let matrix_scale_factor = 100.0 / ungapped.lambda;
    let min_freq = BLAST_SCORE_MIN as f64 / matrix_scale_factor;
    let mut scores = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];

    for q in 0..BLASTAA_SIZE {
        for (s, group) in reverse.iter().enumerate() {
            let mut val = 0.0;
            for &aa in group {
                val += matrix_info.start_freq_ratios[q][aa] * compressed_prob[aa];
            }
            let val = if val < 1.0e-8 { min_freq } else { val.ln() };
            scores[q][s] = blast_nint(val * matrix_scale_factor);
        }
    }

    Some(scores)
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1283-1355
// ```c
// new_alphabet = SCompressedAlphabetNew(sbp,
//                             lookup->compressed_alphabet_size,
//                             kMatrixScale);
// ...
// s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                 query, locations);
// s_CompressedLookupFinalize(lookup);
// ```
// `Some(Err(()))` when the words fill NCBI's overflow banks (`s_CompressedListGetNewCell`),
// past which NCBI writes beyond `overflow_banks`.
pub fn build_blosum62_compressed_lookup(
    word_size: usize,
    threshold: f64,
    query: &[u8],
    locations: &[(i32, i32)],
) -> Option<Result<BlastCompressedAaLookupTable, ()>> {
    build_blosum62_compressed_lookup_with_capacity(
        word_size,
        threshold,
        query,
        locations,
        NCBI_OVERFLOW_CELL_CAPACITY,
    )
}

// `build_blosum62_compressed_lookup` with `overflow_cell_capacity` overflow cells (a parameter
// for the tests).
fn build_blosum62_compressed_lookup_with_capacity(
    word_size: usize,
    threshold: f64,
    query: &[u8],
    locations: &[(i32, i32)],
    overflow_cell_capacity: usize,
) -> Option<Result<BlastCompressedAaLookupTable, ()>> {
    let mut lookup = BlastCompressedAaLookupTable::new(word_size, threshold)?;
    lookup.overflow_cell_capacity = overflow_cell_capacity;
    let matrix = build_blosum62_compressed_score_matrix(lookup.compressed_alphabet_size as usize)?;
    lookup.add_neighboring_words(&matrix, query, locations);
    if lookup.overflow_banks_exhausted {
        return Some(Err(()));
    }
    lookup.finalize();
    Some(Ok(lookup))
}

/// The table of `build_blosum62_compressed_lookup` for a threshold already in NCBI's
/// matrix units (`(Int4)(kMatrixScale * opt->threshold)`, converted by the caller as the
/// x86-64 build does), built once for the queries of several of NCBI's tables: `locations`
/// are the lookup segments of all queries in query order, and `tables` the ranges of
/// `locations` (in order) of each table that NCBI builds (BLASTP's query batches). The
/// table has no limit on its overflow cells, and its cells hold their offsets in NCBI's
/// order. `Some(Err(()))` when the words of one of NCBI's tables fill its overflow banks.
pub fn try_build_blosum62_compressed_lookup(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
) -> Option<Result<BlastCompressedAaLookupTable, ()>> {
    try_build_blosum62_compressed_lookup_with_capacity(
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        NCBI_OVERFLOW_CELL_CAPACITY,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-213
// ```c++
//     // 5. Create the lookup table
//     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
//         LookupTableWrap* lut =
//             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
//                                              sbp, lookup_segments_wrap,
//                                              retval->m_InternalData->m_RpsData,
//                                              seqsrc,
//                                              num_threads);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:798-834
// ```c
//     case 5:
//     {
//         ...
//         CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
//         ...
//     default:
//     {
//         ...
//         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
//
//         if (cell_index == 0 ) { /* can't be empty => it's full  */
//
//             /* fetch next vacant cell */
//         	CompressedOverflowCell * new_cell = s_CompressedListGetNewCell(lookup);
// ```
// NCBI builds a table for each query batch; this is one table for all of them, whose
// cells hold the offsets of each batch in the batch's order, one batch after another.
// `overflow_cell_capacity` is that of NCBI's banks (a parameter for the tests).
//
// The test of each of NCBI's tables: a backbone cell with n offsets takes cells(n) overflow
// cells, 0 for n <= 5 and floor(n / 4) from n = 6 on (one cell for the 6th offset, then one
// whenever (num_entries - 3) & 3 == 0). With c offsets of the batches before it, a batch's
// n offsets take cells(c + n) - cells(c) >= cells(n) cells of this table (floor is
// superadditive), so a batch whose words take no more cells of this table than NCBI's
// banks hold cannot fill them in its own table either; this is checked after each query
// offset. A batch past that is tested exactly once: when no words come before it, this
// table so far is the batch's own; otherwise the batch's own table is built with NCBI's
// capacity (`compressed_table_fills_banks`). Either stops as soon as the banks are full.
fn try_build_blosum62_compressed_lookup_with_capacity(
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
    overflow_cell_capacity: usize,
) -> Option<Result<BlastCompressedAaLookupTable, ()>> {
    let mut lookup = BlastCompressedAaLookupTable::new(word_size, 0.0)?;
    lookup.threshold = scaled_threshold;
    lookup.overflow_cell_capacity = usize::MAX;
    let matrix = build_blosum62_compressed_score_matrix(lookup.compressed_alphabet_size as usize)?;
    let info = CompressedNeighborInfo::new(&matrix, lookup.compressed_alphabet_size as usize);
    let mut next = 0usize;
    for table in tables {
        // The words between NCBI's tables (a split query batch, whose chunks have tables
        // of their own) are added without a test.
        lookup.add_location_words(&info, &matrix, query, &locations[next..table.start], |_| {
            true
        });
        let table_locations = &locations[table.clone()];
        let cells_before = lookup.overflow_cells.len();
        let words_before = table.start > 0;
        let mut tested = false;
        let fits = lookup.add_location_words(&info, &matrix, query, table_locations, |lookup| {
            if tested || lookup.overflow_cells.len() - cells_before <= overflow_cell_capacity {
                return true;
            }
            tested = true;
            words_before
                && !compressed_table_fills_banks(
                    word_size,
                    scaled_threshold,
                    &matrix,
                    query,
                    table_locations,
                    overflow_cell_capacity,
                )
        });
        if !fits {
            return Some(Err(()));
        }
        next = table.end;
    }
    lookup.add_location_words(&info, &matrix, query, &locations[next..], |_| true);
    lookup.finalize();
    Some(Ok(lookup))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
// ```c
//     s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
//                                     query, locations);
//     s_CompressedLookupFinalize(lookup);
// ```
// One of NCBI's tables (a query batch's lookup segments), built with NCBI's capacity only to
// see whether its words fill the banks; it stops there. Its offsets are those of the whole
// query, which changes no count of a cell.
fn compressed_table_fills_banks(
    word_size: usize,
    scaled_threshold: i32,
    matrix: &CompressedScoreMatrix,
    query: &[u8],
    locations: &[(i32, i32)],
    overflow_cell_capacity: usize,
) -> bool {
    let mut table = BlastCompressedAaLookupTable::new(word_size, 0.0)
        .expect("the word size of a table already built");
    table.threshold = scaled_threshold;
    table.overflow_cell_capacity = overflow_cell_capacity;
    let info = CompressedNeighborInfo::new(matrix, table.compressed_alphabet_size as usize);
    table.add_location_words(&info, matrix, query, locations, |table| {
        !table.overflow_banks_exhausted
    });
    table.overflow_banks_exhausted
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1071-1110
// ```c
//     lookup->exact_matches++;
//     score = 0;
//     for (i = 0; i < lookup->word_length; i++) {
//         int c = lookup->compress_table[w[i]];
//         if (c >= lookup->compressed_alphabet_size)
//             return;
//         score += info->matrix[w[i]][c];
//     }
//     ...
//     if (lookup->threshold == 0 || score < lookup->threshold) {
//         s_CompressedLookupAddUnencoded(lookup, w, query_offset);
//     }
//     ...
//     if (lookup->threshold == 0)
//         return;
//     ...
//     score = info->row_max[w[0]];
//     for (i = 1; i < lookup->word_length; i++)
//         score += info->row_max[w[i]];
//     s_CompressedAddWordHitsCore(info, score, 0);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053
// ```c
//     score -= info->row_max[currQueryChar];
//     rowSorted = info->matrixSorted[currQueryChar];
//     ...
//     if (current_pos == info->wordsize - 1) {
//         ...
//         for (i = 0; i < info->compressed_alphabet_size &&
//                (score + rowSorted[i] >= info->threshold); i++) {
//             subject_word[current_pos] = charSorted[i];
//             s_CompressedLookupAddEncoded(lookup, subject_word,
//                                          query_offset);
//         ...
//         return;
//     }
//     ...
//     for (i = 0; i < info->compressed_alphabet_size &&
//            (score + rowSorted[i] >= info->threshold); i++) {
//         subject_word[current_pos] = charSorted[i];
//         s_CompressedAddWordHitsCore(info, score + rowSorted[i],
//                                     current_pos + 1);
// ```
// No NCBI counterpart: counts the words of `add_word_hits` (above) without a table, for LOSAT's
// choice of a path (EXPERIMENT LOSAT_X_BPQBLOCK, `blastp/blast_engine.rs`
// `x_block_search_refuses_nothing`); it changes no value NCBI computes.
/// The number of words that `add_word_hits` adds to a table for one query offset: none when a
/// letter of the word has no compressed letter; otherwise the word itself (`add_unencoded`) when
/// there is no threshold or its own score is below it, and with a threshold the neighbors that
/// `add_word_hits_core` enumerates (the same loops, counted). The enumeration bounds the
/// positions after the current one by their row maxima, which a word reaches, so its neighbors
/// are exactly the compressed words whose score reaches the threshold. Neither their number nor
/// the word's own score depends on the order of its letters, so a count is kept per set of
/// letters (`counts`, keyed by the number of each letter) and the letters of a proteome take a
/// few tens of thousands of enumerations.
pub(crate) struct CompressedWordCount {
    threshold: i32,
    word_length: usize,
    compressed_alphabet_size: usize,
    compress_table: [u8; BLASTAA_SIZE],
    matrix: CompressedScoreMatrix,
    info: CompressedNeighborInfo,
    counts: rustc_hash::FxHashMap<u128, u64>,
}

impl CompressedWordCount {
    /// The counts of the table of `build_blosum62_compressed_lookup` with a threshold in NCBI's
    /// matrix units (as `try_build_blosum62_compressed_lookup` takes it).
    pub(crate) fn new(word_size: usize, scaled_threshold: i32) -> Option<Self> {
        let (compressed_alphabet_size, _) = compressed_lookup_parameters(word_size)?;
        let matrix = build_blosum62_compressed_score_matrix(compressed_alphabet_size)?;
        Some(Self {
            threshold: scaled_threshold,
            word_length: word_size,
            compressed_alphabet_size,
            compress_table: build_compressed_translation(compressed_alphabet_size)?,
            info: CompressedNeighborInfo::new(&matrix, compressed_alphabet_size),
            matrix,
            counts: rustc_hash::FxHashMap::default(),
        })
    }

    /// The words that `add_word_hits` adds for `word` (NCBISTDAA letters, `word_length` of them).
    fn words(&mut self, word: &[u8]) -> u64 {
        let mut letters = 0u128;
        for &residue in word {
            let compressed = self.compress_table.get(usize::from(residue));
            if compressed.is_none_or(|&c| usize::from(c) >= self.compressed_alphabet_size) {
                return 0;
            }
            // 3 bits per letter (BLASTAA_SIZE of them) take the up to 7 letters of a word.
            letters += 1u128 << (3 * u32::from(residue));
        }
        if let Some(&words) = self.counts.get(&letters) {
            return words;
        }
        let words = self.count(word);
        self.counts.insert(letters, words);
        words
    }

    /// `words` without the kept counts.
    fn count(&self, word: &[u8]) -> u64 {
        let mut own_score = 0i32;
        for &residue in word {
            let compressed = self.compress_table[usize::from(residue)];
            own_score += self.matrix[usize::from(residue)][usize::from(compressed)];
        }
        let own = u64::from(self.threshold == 0 || own_score < self.threshold);
        if self.threshold == 0 {
            return own;
        }
        let mut score = self.info.row_max[usize::from(word[0])];
        for &residue in &word[1..] {
            score += self.info.row_max[usize::from(residue)];
        }
        own + self.count_neighbors(word, score, 0)
    }

    /// The words that `add_word_hits_core` adds from `current_pos` on.
    fn count_neighbors(&self, word: &[u8], mut score: i32, current_pos: usize) -> u64 {
        let curr_query_char = usize::from(word[current_pos]);
        score -= self.info.row_max[curr_query_char];
        let row_sorted = &self.info.matrix_sorted[curr_query_char];
        let mut neighbors = 0u64;
        for &letter_score in &row_sorted[..self.compressed_alphabet_size] {
            if score + letter_score < self.threshold {
                break;
            }
            neighbors += if current_pos == self.word_length - 1 {
                1
            } else {
                self.count_neighbors(word, score + letter_score, current_pos + 1)
            };
        }
        neighbors
    }

    /// Whether a table whose words are some of the words of `queries`' offsets (each query's
    /// NCBISTDAA letters) takes at most `overflow_cell_capacity` overflow cells: a backbone cell
    /// with n offsets takes cells(n) <= n / 4 of them
    /// (`try_build_blosum62_compressed_lookup_with_capacity`), so such a table takes at most a
    /// quarter of the queries' words. Stops counting once they are more than that.
    fn fits(
        &mut self,
        queries: &mut dyn Iterator<Item = &[u8]>,
        overflow_cell_capacity: usize,
    ) -> bool {
        let limit = (COMPRESSED_HITS_PER_OVERFLOW_CELL as u64)
            .saturating_mul(overflow_cell_capacity as u64);
        let mut words = 0u64;
        for query in queries {
            for word in query.windows(self.word_length) {
                words = words.saturating_add(self.words(word));
                if words > limit {
                    return false;
                }
            }
        }
        true
    }

    /// Whether no table whose words are some of the words of `queries`' offsets can fill NCBI's
    /// overflow banks: a query batch's table, or the table of one of its query chunks (whose
    /// parts are pieces of the batch's queries). `false` when their words are more than four
    /// times the banks' cells, though the table might still fit.
    pub(crate) fn fits_banks(&mut self, queries: &mut dyn Iterator<Item = &[u8]>) -> bool {
        self.fits(queries, NCBI_OVERFLOW_CELL_CAPACITY)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1339-1347
// ```c
// table_scale = iexp(lookup->compressed_alphabet_size, word_size - 1);
// for (i = 0; i < BLASTAA_SIZE; i++) {
//     Uint1 letter = new_alphabet->compress_table[i];
//     if (letter >= lookup->compressed_alphabet_size)
//         lookup->scaled_compress_table[i] = -1;
//     else
//         lookup->scaled_compress_table[i] = table_scale * letter;
// }
// ```
pub fn build_scaled_compress_table(
    compress_table: &[u8; BLASTAA_SIZE],
    compressed_alphabet_size: usize,
    word_size: usize,
) -> [i32; BLASTAA_SIZE] {
    let table_scale = compressed_alphabet_size.pow((word_size - 1) as u32) as i32;
    let mut scaled = [0i32; BLASTAA_SIZE];
    for i in 0..BLASTAA_SIZE {
        let letter = compress_table[i] as usize;
        scaled[i] = if letter >= compressed_alphabet_size {
            -1
        } else {
            table_scale * letter as i32
        };
    }
    scaled
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1339-1347
// ```c
// lookup->scaled_compress_table = (Int4 *)malloc(
//                                    BLASTAA_SIZE * sizeof(Int4));
// ...
// lookup->scaled_compress_table[i] = table_scale * letter;
// ```
//
// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:270-278
// ```c
// scaled_compress_table = lookup->scaled_compress_table;
// ...
// compressed_char = scaled_compress_table[next_char];
// ```
fn build_scaled_compress_table_u8(scaled_compress_table: &[i32; BLASTAA_SIZE]) -> [i32; 256] {
    let mut table = [-1i32; 256];
    table[..BLASTAA_SIZE].copy_from_slice(scaled_compress_table);
    table
}

// NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:859-909
// ```c
// index = w[0] + W6p1[w[1]] + W6p2[w[2]] + W6p3[w[3]] + W6p4[w[4]];
// ...
// index =  w[0] + W7p1[w[1]] + W7p2[w[2]] + W7p3[w[3]] +
//          W7p4[w[4]] + W7p5[w[5]] + W7p6[w[6]];
// ```
pub fn compute_encoded_compressed_index(
    word_size: usize,
    word: &[u8],
    compressed_alphabet_size: usize,
) -> usize {
    debug_assert!(word.len() >= word_size);
    let a = compressed_alphabet_size;
    let a2 = a * a;
    let a3 = a2 * a;
    let a4 = a3 * a;
    match word_size {
        5 => {
            usize::from(word[0])
                + usize::from(word[1]) * a
                + usize::from(word[2]) * a2
                + usize::from(word[3]) * a3
                + usize::from(word[4]) * a4
        }
        6 => {
            usize::from(word[0])
                + usize::from(word[1]) * a
                + usize::from(word[2]) * a2
                + usize::from(word[3]) * a3
                + usize::from(word[4]) * a4
                + usize::from(word[5]) * a4 * a
        }
        7 => {
            usize::from(word[0])
                + usize::from(word[1]) * a
                + usize::from(word[2]) * a2
                + usize::from(word[3]) * a3
                + usize::from(word[4]) * a4
                + usize::from(word[5]) * a4 * a
                + usize::from(word[6]) * a4 * a * a
        }
        _ => {
            let mut index = 0usize;
            let mut scale = 1usize;
            for &letter in word.iter().take(word_size) {
                index += letter as usize * scale;
                scale *= compressed_alphabet_size;
            }
            index
        }
    }
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:304-322
// ```c
// *skip = 0;
// for(i = 0; i < wordsize; i++) {
//     Int4 ch = scaled_compress_table[word[i]];
//     if (ch < 0){
//         *skip = i + 2;
//         ch = 0;
//     }
//     index = index / compressed_alphabet_size +  ch;
// }
// return index;
// ```
pub fn compute_compressed_index(
    word_size: usize,
    word: &[u8],
    compressed_alphabet_size: usize,
    scaled_compress_table: &[i32; BLASTAA_SIZE],
) -> (i32, i32) {
    let mut skip = 0i32;
    let mut index = 0i32;
    for (i, &residue) in word.iter().take(word_size).enumerate() {
        let mut ch = scaled_compress_table
            .get(residue as usize)
            .copied()
            .unwrap_or(-1);
        if ch < 0 {
            skip = i32::try_from(i + 2).expect("NCBI BLAST compressed skip must fit in Int4");
            ch = 0;
        }
        index = index / compressed_alphabet_size as i32 + ch;
    }
    (index, skip)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h:304-322
// ```c
// *skip = 0;
// for(i = 0; i < wordsize; i++) {
//     Int4 ch = scaled_compress_table[word[i]];
//     if (ch < 0){
//         *skip = i + 2;
//         ch = 0;
//     }
//     index = index / compressed_alphabet_size +  ch;
// }
// return index;
// ```
#[inline(always)]
pub(crate) fn compute_compressed_index_u8(
    word_size: usize,
    word: &[u8],
    compressed_alphabet_size: usize,
    scaled_compress_table: &[i32; 256],
) -> (i32, i32) {
    let mut skip = 0i32;
    let mut index = 0i32;
    for (i, &residue) in word.iter().take(word_size).enumerate() {
        let mut ch = scaled_compress_table[usize::from(residue)];
        if ch < 0 {
            skip = i32::try_from(i + 2).expect("NCBI BLAST compressed skip must fit in Int4");
            ch = 0;
        }
        index = index / compressed_alphabet_size as i32 + ch;
    }
    (index, skip)
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:49-57
// ```c
// #define PV_SET(lookup, index, pv_array_bts) \
//     lookup[(index) >> (pv_array_bts)] |= (1U << ((index) & PV_ARRAY_MASK))
// #define PV_TEST(lookup, index, pv_array_bts) \
//     (lookup[(index) >> (pv_array_bts)] & (1U << ((index) & PV_ARRAY_MASK)))
// ```
#[inline(always)]
pub fn compressed_pv_test(pv: &[u32], index: usize, pv_array_bts: usize) -> bool {
    (pv[index >> pv_array_bts] & (1u32 << (index & PV_ARRAY_MASK))) != 0
}

// NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_lookup.h:49-57
// ```c
// #define PV_SET(lookup, index, pv_array_bts) \
//     lookup[(index) >> (pv_array_bts)] |= (1U << ((index) & PV_ARRAY_MASK))
// ```
#[inline(always)]
pub fn compressed_pv_set(pv: &mut [u32], index: usize, pv_array_bts: usize) {
    pv[index >> pv_array_bts] |= 1u32 << (index & PV_ARRAY_MASK);
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::matrix::ncbistdaa;

    // A small deterministic generator (xorshift64) for the tests of the query batches' tables.
    struct XorShift(u64);

    impl XorShift {
        fn below(&mut self, n: usize) -> usize {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            (x % n as u64) as usize
        }
    }

    // Queries of `lengths` letters drawn from `letters`, concatenated as BLASTP's lookup query
    // (each query followed by a sentinel), with one lookup segment per query.
    fn concatenated_queries(
        rng: &mut XorShift,
        lengths: &[usize],
        letters: &[u8],
    ) -> (Vec<u8>, Vec<(i32, i32)>) {
        let mut query = Vec::new();
        let mut locations = Vec::new();
        for &length in lengths {
            let left = query.len() as i32;
            for _ in 0..length {
                query.push(letters[rng.below(letters.len())]);
            }
            locations.push((left, left + length as i32 - 1));
            query.push(0);
        }
        (query, locations)
    }

    // A table of `locations` with no limit on its overflow cells (not finalized).
    fn unlimited_table(
        threshold: i32,
        query: &[u8],
        locations: &[(i32, i32)],
    ) -> BlastCompressedAaLookupTable {
        let mut table = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
        table.threshold = threshold;
        table.overflow_cell_capacity = usize::MAX;
        let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
        table.add_neighboring_words(&matrix, query, locations);
        table
    }

    fn assert_same_table(lhs: &BlastCompressedAaLookupTable, rhs: &BlastCompressedAaLookupTable) {
        assert!(lhs.backbone == rhs.backbone, "backbone cells differ");
        assert!(
            lhs.overflow_cells == rhs.overflow_cells,
            "overflow cells differ"
        );
        assert_eq!(lhs.pv, rhs.pv);
        assert_eq!(lhs.pv_array_bts, rhs.pv_array_bts);
        assert_eq!(lhs.longest_chain, rhs.longest_chain);
        assert_eq!(lhs.neighbor_matches, rhs.neighbor_matches);
        assert_eq!(lhs.exact_matches, rhs.exact_matches);
        for index in 0..lhs.backbone.len() {
            if lhs.backbone[index].num_used > 5 {
                assert_eq!(
                    lhs.query_offsets_in_scan_order(index),
                    rhs.query_offsets_in_scan_order(index)
                );
            }
        }
    }

    // The queries of the tests below: few letters and a low threshold (11, in NCBI's matrix
    // units), so that words repeat and cells overflow.
    const FEW_LETTERS: [u8; 4] = [ncbistdaa::A, ncbistdaa::L, ncbistdaa::S, ncbistdaa::K];
    const LOW_THRESHOLD: i32 = 1100;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-213
    // ```c++
    //     if ( !retval->m_QuerySplitter->IsQuerySplit() ) {
    //         LookupTableWrap* lut =
    //             CSetupFactory::CreateLookupTable(query_data, opts_memento.get(),
    // ```
    // (a) The table of all batches takes more cells than the capacity, no batch's own table
    // does: the table is built, and it is the table without a limit.
    #[test]
    fn compressed_batches_within_capacity_give_the_unlimited_table() {
        let mut rng = XorShift(0x9e37_79b9_7f4a_7c15);
        let (query, locations) =
            concatenated_queries(&mut rng, &[300, 250, 400, 150, 350, 300], &FEW_LETTERS);
        let tables = [0..2, 2..3, 3..6];
        let batch_cells: Vec<usize> = tables
            .iter()
            .map(|table| {
                unlimited_table(LOW_THRESHOLD, &query, &locations[table.clone()])
                    .overflow_cells
                    .len()
            })
            .collect();
        let capacity = *batch_cells.iter().max().unwrap();
        let mut unlimited = unlimited_table(LOW_THRESHOLD, &query, &locations);
        assert!(capacity > 0);
        assert!(unlimited.overflow_cells.len() > capacity);
        unlimited.finalize();

        let table = try_build_blosum62_compressed_lookup_with_capacity(
            5,
            LOW_THRESHOLD,
            &query,
            &locations,
            &tables,
            capacity,
        )
        .expect("word size 5")
        .expect("no batch fills the capacity");
        assert_same_table(&table, &unlimited);

        // Locations outside NCBI's tables (a split batch) are in the table, untested.
        let table = try_build_blosum62_compressed_lookup_with_capacity(
            5,
            LOW_THRESHOLD,
            &query,
            &locations,
            &[0..2, 3..6],
            capacity,
        )
        .expect("word size 5")
        .expect("no batch fills the capacity");
        assert_same_table(&table, &unlimited);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
    // ```c
    //     if (lookup->curr_overflow_cell ==
    //                         COMPRESSED_OVERFLOW_CELLS_IN_BANK) {
    //         /* need a new bank */
    //         Int4 bank_idx = lookup->curr_overflow_bank + 1;
    //         lookup->overflow_banks[bank_idx] = (CompressedOverflowCell*) malloc(
    // ```
    // (b) One batch's own table needs more cells than the capacity: refused, whether the
    // batch is the first (its table so far is the whole table) or a later one.
    #[test]
    fn compressed_batch_past_capacity_is_refused() {
        let mut rng = XorShift(0x2545_f491_4f6c_dd1d);
        let (query, locations) = concatenated_queries(&mut rng, &[100, 400, 60, 80], &FEW_LETTERS);
        let small = [0..1, 2..3, 3..4]
            .iter()
            .map(|table| {
                unlimited_table(LOW_THRESHOLD, &query, &locations[table.clone()])
                    .overflow_cells
                    .len()
            })
            .max()
            .unwrap();
        let large = unlimited_table(LOW_THRESHOLD, &query, &locations[1..2])
            .overflow_cells
            .len();
        assert!(small > 0 && large > small + 1);
        let build = |tables: &[std::ops::Range<usize>], capacity: usize| {
            try_build_blosum62_compressed_lookup_with_capacity(
                5,
                LOW_THRESHOLD,
                &query,
                &locations,
                tables,
                capacity,
            )
            .expect("word size 5")
        };
        // A later batch: its own table is built for the test.
        assert!(build(&[0..1, 1..2, 2..3, 3..4], large - 1).is_err());
        assert!(build(&[0..1, 1..2, 2..3, 3..4], large).is_ok());
        assert!(build(&[1..2, 2..4], large - 1).is_err());
        // The first batch: the table so far is its own.
        let first = unlimited_table(LOW_THRESHOLD, &query, &locations[0..2])
            .overflow_cells
            .len();
        let rest = unlimited_table(LOW_THRESHOLD, &query, &locations[2..4])
            .overflow_cells
            .len();
        assert!(first >= large && rest < first - 1);
        assert!(build(&[0..2, 2..4], first - 1).is_err());
        assert!(build(&[0..2, 2..4], first).is_ok());
        // The large query in a split batch (no table of its own) is not tested.
        assert!(build(&[0..1, 2..3, 3..4], small).is_ok());
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:825-841
    // ```c
    //     default:
    //     {
    //     	/* continue with existing overflow list */
    //         ...
    //         Int4 cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
    //
    //         if (cell_index == 0 ) { /* can't be empty => it's full  */
    // ```
    // (c) On random queries: a batch's own table takes cells(n) = 0 (n <= 5) or floor(n / 4)
    // cells per backbone cell, the batch takes at least as many cells of the table of all
    // batches (the shortcut is sound), and the build refuses exactly when some batch's own
    // table needs more cells than the capacity.
    #[test]
    fn compressed_batch_shortcut_agrees_with_exact_test() {
        let mut rng = XorShift(0x1234_5678_9abc_def1);
        let mut shortcut_undecided = 0usize;
        for round in 0..12 {
            let letters: &[u8] = if round % 3 == 2 {
                &[
                    ncbistdaa::A,
                    ncbistdaa::L,
                    ncbistdaa::S,
                    ncbistdaa::K,
                    ncbistdaa::E,
                    ncbistdaa::G,
                    ncbistdaa::X,
                ]
            } else {
                &FEW_LETTERS
            };
            let threshold = [1100, 1600, 2000, 2400][rng.below(4)];
            let lengths: Vec<usize> = (0..2 + rng.below(6)).map(|_| 5 + rng.below(400)).collect();
            let (query, locations) = concatenated_queries(&mut rng, &lengths, letters);
            // Batches of consecutive queries; some are left out (split batches).
            let mut tables = Vec::new();
            let mut start = 0;
            while start < locations.len() {
                let end = (start + 1 + rng.below(3)).min(locations.len());
                if rng.below(4) != 0 {
                    tables.push(start..end);
                }
                start = end;
            }
            let mut whole = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
            whole.threshold = threshold;
            whole.overflow_cell_capacity = usize::MAX;
            let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
            let mut next = 0;
            let mut exact = Vec::new();
            for table in &tables {
                whole.add_neighboring_words(&matrix, &query, &locations[next..table.start]);
                let before = whole.overflow_cells.len();
                whole.add_neighboring_words(&matrix, &query, &locations[table.clone()]);
                let shortcut = whole.overflow_cells.len() - before;
                let own = unlimited_table(threshold, &query, &locations[table.clone()]);
                let cells: usize = own
                    .backbone
                    .iter()
                    .map(|cell| match cell.num_used {
                        n if n <= 5 => 0,
                        n => n as usize / 4,
                    })
                    .sum();
                assert_eq!(own.overflow_cells.len(), cells);
                assert!(shortcut >= cells, "round {round}: {shortcut} < {cells}");
                exact.push((cells, shortcut));
                next = table.end;
            }
            whole.add_neighboring_words(&matrix, &query, &locations[next..]);
            whole.finalize();

            let mut capacities = vec![1usize, whole.overflow_cells.len() + 1];
            for &(cells, shortcut) in &exact {
                capacities.extend([cells, cells + 1, shortcut, shortcut + 1]);
                if cells > 1 {
                    capacities.push(cells - 1);
                }
            }
            capacities.retain(|&capacity| capacity > 0);
            for capacity in capacities {
                let refused = exact.iter().any(|&(cells, _)| cells > capacity);
                shortcut_undecided += exact
                    .iter()
                    .filter(|&&(cells, shortcut)| cells <= capacity && shortcut > capacity)
                    .count();
                let built = try_build_blosum62_compressed_lookup_with_capacity(
                    5, threshold, &query, &locations, &tables, capacity,
                )
                .expect("word size 5");
                match built {
                    Err(()) => assert!(refused, "round {round}, capacity {capacity}"),
                    Ok(table) => {
                        assert!(!refused, "round {round}, capacity {capacity}");
                        assert_same_table(&table, &whole);
                    }
                }
            }
        }
        // Batches that the shortcut leaves to the exact test, and that pass it, were seen.
        assert!(shortcut_undecided > 0);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1071-1110
    // ```c
    //     lookup->exact_matches++;
    //     score = 0;
    //     for (i = 0; i < lookup->word_length; i++) {
    //         int c = lookup->compress_table[w[i]];
    //         if (c >= lookup->compressed_alphabet_size)
    //             return;
    // ```
    // (d) `CompressedWordCount` gives the number of words that `add_word_hits` adds for each
    // offset, on random words of all NCBISTDAA letters (some without a compressed letter), with
    // and without a threshold, and the count of a word does not depend on the order of its
    // letters (the counts are kept per set of letters).
    #[test]
    fn compressed_word_count_is_the_words_of_each_offset() {
        let mut rng = XorShift(0x0bad_cafe_1234_5678);
        let standard = (*b"ACDEFGHIKLMNPQRSTVWY").map(aa_char_to_ncbistdaa);
        let mut letters: Vec<u8> = (1..BLASTAA_SIZE as u8).collect();
        for _ in 0..3 {
            letters.extend_from_slice(&standard);
        }
        let matrix = build_blosum62_compressed_score_matrix(15).expect("compressed matrix");
        let info = CompressedNeighborInfo::new(&matrix, 15);
        for threshold in [0, 1100, 1600, 1930, 2000, 2400] {
            let mut counter = CompressedWordCount::new(5, threshold).expect("word size 5");
            let mut table = BlastCompressedAaLookupTable::new(5, 0.0).expect("compressed lookup");
            table.threshold = threshold;
            table.overflow_cell_capacity = usize::MAX;
            let query: Vec<u8> = (0..40).map(|_| letters[rng.below(letters.len())]).collect();
            let mut words_before = 0u64;
            for offset in 0..=query.len() - 5 {
                table.add_word_hits(&info, &matrix, &query, offset as i32);
                let words_after: u64 = table.backbone.iter().map(|cell| cell.num_used as u64).sum();
                let added = words_after - words_before;
                words_before = words_after;
                let word = &query[offset..offset + 5];
                assert_eq!(
                    counter.words(word),
                    added,
                    "threshold {threshold}, offset {offset}"
                );
                if word
                    .iter()
                    .all(|&residue| usize::from(table.compress_table[usize::from(residue)]) < 15)
                {
                    let mut reversed = word.to_vec();
                    reversed.reverse();
                    assert_eq!(counter.count(&reversed), added, "reversed, offset {offset}");
                }
            }
        }
    }

    // (e) The count clears a table exactly when a quarter of its words fits the capacity, and a
    // table it clears is not refused (`try_build_blosum62_compressed_lookup_with_capacity`).
    #[test]
    fn compressed_word_count_clears_only_tables_within_the_capacity() {
        let mut rng = XorShift(0x5851_f42d_4c95_7f2d);
        let (query, locations) = concatenated_queries(&mut rng, &[120, 300, 80, 200], &FEW_LETTERS);
        let queries: Vec<&[u8]> = locations
            .iter()
            .map(|&(left, right)| &query[left as usize..=right as usize])
            .collect();
        let mut counter = CompressedWordCount::new(5, LOW_THRESHOLD).expect("word size 5");
        let words: u64 = queries
            .iter()
            .flat_map(|query| query.windows(5))
            .map(|word| counter.words(word))
            .sum();
        let table = unlimited_table(LOW_THRESHOLD, &query, &locations);
        let table_words: u64 = table.backbone.iter().map(|cell| cell.num_used as u64).sum();
        assert_eq!(words, table_words);
        let cells = table.overflow_cells.len();
        let cleared = words.div_ceil(4) as usize;
        assert!(cells > 0 && cells <= cleared);
        assert!(counter.fits(&mut queries.iter().copied(), cleared));
        assert!(!counter.fits(&mut queries.iter().copied(), cleared - 1));
        let table = try_build_blosum62_compressed_lookup_with_capacity(
            5,
            LOW_THRESHOLD,
            &query,
            &locations,
            &[0..locations.len()],
            cleared,
        )
        .expect("word size 5");
        assert!(table.is_ok());
        // The public test uses NCBI's capacity.
        assert!(counter.fits_banks(&mut queries.iter().copied()));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:400-405
    // ```c
    //           else
    //               /* not enough space in the destination array */
    //           {
    //               s_range[1] = (Int4)(s - subject->sequence);
    //               return totalhits;
    //           }
    // ```
    // (f) BLASTP's block search (EXPERIMENT LOSAT_X_BPQBLOCK, `blastp/x_query_blocks.rs`): the
    // table of some queries (a block) gives, at each subject offset, the hits of those queries
    // that the table of all queries gives, not always in the same order (cells of more than 5
    // offsets); and a scan gives the same hits in the same order with any offset-pair array
    // that takes its longest cell, wherever it stops and resumes.
    #[test]
    fn compressed_block_table_scans_the_whole_tables_hits_per_subject_offset() {
        let mut rng = XorShift(0x6a09_e667_f3bc_c908);
        let (query, locations) =
            concatenated_queries(&mut rng, &[150, 90, 200, 60, 120], &FEW_LETTERS);
        let mut whole = unlimited_table(LOW_THRESHOLD, &query, &locations);
        whole.finalize();
        let block = 1..3;
        let mut part = unlimited_table(LOW_THRESHOLD, &query, &locations[block.clone()]);
        part.finalize();
        let block_offsets = locations[block.start].0..=locations[block.end - 1].1;
        let mut subject: Vec<u8> = (0..400)
            .map(|_| FEW_LETTERS[rng.below(FEW_LETTERS.len())])
            .collect();
        subject.push(0);
        let scan = |table: &BlastCompressedAaLookupTable, array_size: i32| {
            let seq_ranges = [(0, subject.len() as i32 - 1)];
            let mut s_range = [0, 0, seq_ranges[0].1 - 5];
            let mut offset_pairs = vec![OffsetPair::default(); array_size as usize];
            let mut hits = Vec::new();
            let mut stops = 0;
            while s_range[1] <= s_range[2] {
                let count = table.scan_subject(
                    &subject,
                    &seq_ranges,
                    &mut offset_pairs,
                    array_size,
                    &mut s_range,
                );
                hits.extend(
                    offset_pairs[..count as usize]
                        .iter()
                        .map(|pair| (pair.s_off, pair.q_off as i32)),
                );
                stops += 1;
            }
            (hits, stops)
        };
        let (whole_hits, _) = scan(&whole, 4096 + whole.longest_chain);
        let (part_hits, _) = scan(&part, 4096 + part.longest_chain);
        let whole_in_block: Vec<(u32, i32)> = whole_hits
            .iter()
            .copied()
            .filter(|(_, q_off)| block_offsets.contains(q_off))
            .collect();
        assert_ne!(
            part_hits, whole_in_block,
            "some subject offset's order differs"
        );
        let mut sorted_part = part_hits.clone();
        sorted_part.sort_unstable();
        let mut sorted_whole = whole_in_block.clone();
        sorted_whole.sort_unstable();
        assert_eq!(sorted_part, sorted_whole);
        assert!(part_hits.windows(2).all(|pair| pair[0].0 <= pair[1].0));
        for array_size in [
            part.longest_chain,
            part.longest_chain + 7,
            3 * part.longest_chain,
        ] {
            let (hits, stops) = scan(&part, array_size);
            assert!(stops > 1, "array {array_size}: the scan stops and resumes");
            assert_eq!(hits, part_hits, "array {array_size}");
        }
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:782-849
    // ```c
    // case 5:
    //     new_cell->query_offsets[0] = backbone_cell->payload.query_offsets[2];
    //     new_cell->query_offsets[1] = backbone_cell->payload.query_offsets[3];
    //     new_cell->query_offsets[2] = query_offset;
    // ...
    // default:
    //     cell_index = (num_entries -3) & COMPRESSED_HITS_CELL_MASK ;
    //     if (cell_index == 0 ) {
    //         new_cell->next = backbone_cell->payload.overflow_list.head;
    //         backbone_cell->payload.overflow_list.head = new_cell;
    //     }
    // ```
    // `(Int4)(kMatrixScale * opt->threshold)` past `Int4` is `INT_MIN` (every word is a
    // neighbor); words past the overflow banks make the table `Err` (NCBI writes past its banks).
    #[test]
    fn a_threshold_past_int4_takes_every_word_up_to_the_banks() {
        assert_eq!(
            BlastCompressedAaLookupTable::new(5, 19.3)
                .unwrap()
                .threshold,
            1930
        );
        for threshold in [21_474_836.48, 1e10, f64::INFINITY] {
            let lookup = BlastCompressedAaLookupTable::new(5, threshold).unwrap();
            assert_eq!(lookup.threshold, i32::MIN, "{threshold}");
        }
        let query = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0];
        let words = 15usize.pow(5);
        let lookup = build_blosum62_compressed_lookup(5, 1e10, &query, &[(1, 6)])
            .unwrap()
            .expect("two offsets take no overflow cell");
        let offsets: i32 = lookup.backbone.iter().map(|cell| cell.num_used).sum();
        assert_eq!(offsets as usize, 2 * words);
        assert!(lookup.overflow_cells.is_empty());
        // Seven offsets take one overflow cell in each backbone cell.
        assert!(
            build_blosum62_compressed_lookup_with_capacity(5, 1e10, &query, &[(1, 11)], 1000)
                .unwrap()
                .is_err()
        );
    }

    #[test]
    fn compressed_word_hit_overflow_matches_ncbi_scan_order() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        for offset in 10..18 {
            lookup.add_word_hit(42, offset);
        }

        assert_eq!(lookup.backbone[42].num_used, 8);
        assert_eq!(lookup.overflow_cells.len(), 2);
        assert_eq!(
            lookup.query_offsets_in_scan_order(42),
            vec![10, 11, 12, 17, 13, 14, 15, 16]
        );
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1182-1242
    // ```c
    // PV_SET(pv, i, pv_array_bts);
    // longest_chain = MAX(count, longest_chain);
    // ```
    #[test]
    fn compressed_finalize_sets_pv_and_longest_chain() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        lookup.add_word_hit(7, 100);
        lookup.add_word_hit(33, 200);
        lookup.add_word_hit(33, 201);
        lookup.finalize();

        assert_eq!(lookup.longest_chain, 2);
        assert!(compressed_pv_test(
            &lookup.pv,
            7,
            lookup.pv_array_bts as usize
        ));
        assert!(compressed_pv_test(
            &lookup.pv,
            33,
            lookup.pv_array_bts as usize
        ));
        assert!(!compressed_pv_test(
            &lookup.pv,
            34,
            lookup.pv_array_bts as usize
        ));
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:859-929
    // ```c
    // s_CompressedLookupAddEncoded(lookup, subject_word, query_offset);
    // ...
    // if (skip == 0)
    //     s_CompressedLookupAddWordHit(lookup, index, query_offset);
    // ```
    #[test]
    fn compressed_encoded_and_unencoded_indices_match() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        let word = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
        ];
        let encoded_word: Vec<u8> = word
            .iter()
            .map(|&residue| lookup.compress_table[residue as usize])
            .collect();
        let encoded_index = compute_encoded_compressed_index(5, &encoded_word, 15);
        let (unencoded_index, skip) =
            compute_compressed_index(5, &word, 15, &lookup.scaled_compress_table);

        assert_eq!(skip, 0);
        assert_eq!(encoded_index as i32, unencoded_index);
        lookup.add_unencoded(&word, 321);
        assert_eq!(lookup.query_offsets_in_scan_order(encoded_index), vec![321]);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:322-384
    // ```c
    // if (PV_TEST(pv, index, pv_array_bts)) {
    //     Int4 s_off = (Int4)(s - subject->sequence);
    //     ...
    //     dest[0].qs_offsets.q_off = backbone_cell->query_offset;
    //     dest[0].qs_offsets.s_off = s_off;
    //     ...
    //     totalhits += numhits;
    // }
    // ```
    #[test]
    fn compressed_scan_subject_emits_ncbi_subject_major_order() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        let word = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
        ];
        lookup.add_unencoded(&word, 10);
        lookup.add_unencoded(&word, 11);
        lookup.finalize();

        let subject = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            0,
        ];
        let seq_ranges = [(0, 10)];
        let mut s_range = [0, 0, 5];
        let mut offset_pairs = vec![OffsetPair::default(); 8];
        let hits = lookup.scan_subject(&subject, &seq_ranges, &mut offset_pairs, 8, &mut s_range);

        assert_eq!(hits, 4);
        assert_eq!(
            offset_pairs[..hits as usize]
                .iter()
                .map(|pair| (pair.q_off, pair.s_off))
                .collect::<Vec<_>>(),
            vec![(10, 0), (11, 0), (10, 5), (11, 5)]
        );
        assert_eq!(s_range[1], 6);
    }

    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aascan.c:373-378
    // ```c
    // else {
    //     s_range[1] = (Int4)(s - subject->sequence);
    //     return totalhits;
    // }
    // ```
    #[test]
    fn compressed_scan_subject_small_buffer_resumes_without_splitting_subject_offset() {
        let mut lookup = BlastCompressedAaLookupTable::new(5, 19.3).expect("compressed lookup");
        let word = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
        ];
        for query_offset in 10..16 {
            lookup.add_unencoded(&word, query_offset);
        }
        lookup.finalize();

        let subject = [
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            ncbistdaa::X,
            ncbistdaa::S,
            ncbistdaa::T,
            ncbistdaa::A,
            ncbistdaa::G,
            ncbistdaa::W,
            0,
        ];
        let seq_ranges = [(0, 11)];
        let mut s_range = [0, 0, 6];
        let mut offset_pairs = vec![OffsetPair::default(); 6];
        let mut observed = Vec::new();

        while s_range[1] <= s_range[2] {
            let hits =
                lookup.scan_subject(&subject, &seq_ranges, &mut offset_pairs, 6, &mut s_range);
            assert_eq!(hits, 6);
            observed.extend(
                offset_pairs[..hits as usize]
                    .iter()
                    .map(|pair| (pair.q_off, pair.s_off)),
            );
        }

        assert_eq!(
            observed,
            vec![
                (10, 0),
                (11, 0),
                (12, 0),
                (13, 0),
                (14, 0),
                (15, 0),
                (10, 6),
                (11, 6),
                (12, 6),
                (13, 6),
                (14, 6),
                (15, 6),
            ]
        );
    }
}
