//! Compressed amino-acid lookup primitives for the future `blastp-fast` path.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c
//! Reference: ncbi-blast/c++/include/algo/blast/core/blast_aalookup.h

use super::{ilog2, PV_ARRAY_BTS, PV_ARRAY_MASK};
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::build_matrix_info;
use crate::stats::{compute_std_aa_composition, lookup_protein_params_ungapped};
use crate::utils::matrix::{aa_char_to_ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

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
    pub neighbor_matches: i32,
    pub exact_matches: i32,
}

#[derive(Clone, Debug)]
pub struct CompressedNeighborInfo {
    row_max: [i32; BLASTAA_SIZE],
    matrix_sorted: Vec<Vec<i32>>,
    matrix_sorted_char: Vec<Vec<u8>>,
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
impl BlastCompressedAaLookupTable {
    pub fn new(word_size: usize, threshold: f64) -> Option<Self> {
        let (compressed_alphabet_size, reciprocal_alphabet_size) =
            compressed_lookup_parameters(word_size)?;
        let compress_table = build_compressed_translation(compressed_alphabet_size)?;
        let scaled_compress_table =
            build_scaled_compress_table(&compress_table, compressed_alphabet_size, word_size);
        let backbone_size = compressed_alphabet_size
            .checked_pow(word_size as u32)?
            .checked_add(1)?;

        Some(Self {
            threshold: (100.0 * threshold) as i32,
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
            neighbor_matches: 0,
            exact_matches: 0,
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
    fn compressed_list_get_new_cell(&mut self) -> usize {
        let cell_index = self.overflow_cells.len();
        assert!(
            cell_index < COMPRESSED_OVERFLOW_CELLS_IN_BANK * COMPRESSED_OVERFLOW_MAX_BANKS,
            "NCBI BLAST compressed overflow bank capacity exceeded"
        );
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
                let CompressedBackbonePayload::Inline(query_offsets) =
                    self.backbone[index].payload.clone()
                else {
                    unreachable!("NCBI BLAST overflow list cannot exist before entry 6");
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
        matrix: &[Vec<i32>],
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
        matrix: &[Vec<i32>],
        query: &[u8],
        locations: &[(i32, i32)],
    ) {
        let info = CompressedNeighborInfo::new(matrix, self.compressed_alphabet_size as usize);
        let word_size = self.word_length;
        for &(left, right) in locations {
            let to = right - word_size + 1;
            for offset in left..=to {
                let start = usize::try_from(offset)
                    .expect("NCBI BLAST compressed query offset must be non-negative");
                let end = start + word_size as usize;
                if end <= query.len() {
                    self.add_word_hits(&info, matrix, query, offset);
                }
            }
        }
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
    pub fn new(matrix: &[Vec<i32>], compressed_alphabet_size: usize) -> Self {
        let mut row_max = [0i32; BLASTAA_SIZE];
        for row in 0..BLASTAA_SIZE {
            row_max[row] = matrix[row][0];
            for col in 1..compressed_alphabet_size {
                row_max[row] = row_max[row].max(matrix[row][col]);
            }
        }

        let mut matrix_sorted = vec![vec![0i32; compressed_alphabet_size]; BLASTAA_SIZE];
        let mut matrix_sorted_char = vec![vec![0u8; compressed_alphabet_size]; BLASTAA_SIZE];
        for long_char in 0..BLASTAA_SIZE {
            let mut sort_table: Vec<(i32, u8)> = (0..compressed_alphabet_size)
                .map(|short_char| {
                    (
                        row_max[long_char] - matrix[long_char][short_char],
                        short_char as u8,
                    )
                })
                .collect();
            sort_table.sort_by(|lhs, rhs| lhs.0.cmp(&rhs.0));
            for (i, &(_, letter)) in sort_table.iter().enumerate() {
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
) -> Option<Vec<Vec<i32>>> {
    let reverse = build_compressed_reverse_lookup(compressed_alphabet_size)?;
    let compressed_prob = compressed_probabilities(compressed_alphabet_size)?;
    let ungapped = lookup_protein_params_ungapped(ScoringMatrix::Blosum62);
    let matrix_info = build_matrix_info(ScoringMatrix::Blosum62, ungapped.lambda).ok()?;
    let matrix_scale_factor = 100.0 / ungapped.lambda;
    let min_freq = BLAST_SCORE_MIN as f64 / matrix_scale_factor;
    let mut scores = vec![vec![0i32; compressed_alphabet_size]; BLASTAA_SIZE];

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
pub fn build_blosum62_compressed_lookup(
    word_size: usize,
    threshold: f64,
    query: &[u8],
    locations: &[(i32, i32)],
) -> Option<BlastCompressedAaLookupTable> {
    let mut lookup = BlastCompressedAaLookupTable::new(word_size, threshold)?;
    let matrix = build_blosum62_compressed_score_matrix(lookup.compressed_alphabet_size as usize)?;
    lookup.add_neighboring_words(&matrix, query, locations);
    lookup.finalize();
    Some(lookup)
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
    let mut index = 0usize;
    let mut scale = 1usize;
    for &letter in word.iter().take(word_size) {
        index += letter as usize * scale;
        scale *= compressed_alphabet_size;
    }
    index
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
}
