//! Compressed amino-acid lookup primitives for the future `blastp-fast` path.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c
//! NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h

use super::{ilog2, PV_ARRAY_BTS, PV_ARRAY_MASK};
use crate::algorithm::tblastx::scan::OffsetPair;
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::build_matrix_info;
use crate::stats::{compute_std_aa_composition, lookup_protein_params_ungapped};
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::{aa_char_to_ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

// EXPERIMENT (LOSAT_X_BPCLUT; BPCLUT is a strict-set switch): this file's table built flat
// (count, layout, fill), each cell's hits in one slice in NCBI's scan order; a child module, so that
// it uses this file's private builder pieces (x_clut_direct.rs, exactness argument in its header).
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1352-1354
#[path = "x_clut_direct.rs"]
pub mod x_clut_direct;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:186-201
pub const COMPRESSED_HITS_PER_BACKBONE_CELL: usize = 4;
pub const COMPRESSED_HITS_CELL_MASK: usize = 0x03;
pub const COMPRESSED_HITS_PER_OVERFLOW_CELL: usize = 4;
pub const COMPRESSED_OVERFLOW_CELLS_IN_BANK: usize = 209_710;
pub const COMPRESSED_OVERFLOW_MAX_BANKS: usize = 1024;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1306
pub type CompressedScoreMatrix = [[i32; BLASTAA_SIZE]; BLASTAA_SIZE];

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:203-208
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CompressedOverflowCell {
    pub next: Option<usize>,
    pub query_offsets: [i32; COMPRESSED_HITS_PER_OVERFLOW_CELL],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:763-775
impl Default for CompressedOverflowCell {
    fn default() -> Self {
        Self {
            next: None,
            query_offsets: [0; COMPRESSED_HITS_PER_OVERFLOW_CELL],
        }
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:210-231
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CompressedBackbonePayload {
    Inline([i32; COMPRESSED_HITS_PER_BACKBONE_CELL]),
    Overflow {
        query_offsets: [i32; COMPRESSED_HITS_PER_BACKBONE_CELL - 2],
        head: Option<usize>,
    },
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:787-801
impl Default for CompressedBackbonePayload {
    fn default() -> Self {
        Self::Inline([0; COMPRESSED_HITS_PER_BACKBONE_CELL])
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:221-231
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct CompressedLookupBackboneCell {
    pub num_used: i32,
    pub query_offset: i32,
    pub payload: CompressedBackbonePayload,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:236-267
// `alphabet_size` and `backbone_size` are not stored: nothing reads them.
#[derive(Clone, Debug)]
pub struct BlastCompressedAaLookupTable {
    pub threshold: i32,
    pub word_length: i32,
    pub compressed_alphabet_size: i32,
    pub reciprocal_alphabet_size: i32,
    pub longest_chain: i32,
    pub backbone: Vec<CompressedLookupBackboneCell>,
    pub overflow_cells: Vec<CompressedOverflowCell>,
    pub pv: Vec<u32>,
    pub pv_array_bts: i32,
    pub compress_table: [u8; BLASTAA_SIZE],
    pub scaled_compress_table: [i32; BLASTAA_SIZE],
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:269-291
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
// The cells a table can take before `s_CompressedListGetNewCell` writes past `overflow_banks`.
const NCBI_OVERFLOW_CELL_CAPACITY: usize =
    COMPRESSED_OVERFLOW_CELLS_IN_BANK * COMPRESSED_OVERFLOW_MAX_BANKS;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:935-952
#[derive(Clone, Debug)]
pub struct CompressedNeighborInfo {
    row_max: [i32; BLASTAA_SIZE],
    matrix_sorted: CompressedScoreMatrix,
    matrix_sorted_char: [[u8; BLASTAA_SIZE]; BLASTAA_SIZE],
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1289-1320
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1292
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
            compressed_alphabet_size: compressed_alphabet_size as i32,
            reciprocal_alphabet_size,
            longest_chain: 0,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:753-775
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:759-768
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:782-849
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:859-909
    pub fn add_encoded(&mut self, word: &[u8], query_offset: i32) {
        let index = compute_encoded_compressed_index(
            self.word_length as usize,
            word,
            self.compressed_alphabet_size as usize,
        );
        self.add_word_hit(index, query_offset);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:919-929
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1071-1121
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1008-1053
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1137-1176
    #[cfg(test)]
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:343-384
    #[cfg(test)]
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:343-384
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:229-384
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
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:278-311
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
                    // NCBI's decrement of `s` is immediately balanced by the outer `for`
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
                            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:373-378
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1182-1246
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:284
#[inline(always)]
#[allow(dead_code)]
pub(crate) fn compressed_reciprocal_divide(index: i32, reciprocal_alphabet_size: i64) -> i32 {
    ((i64::from(index) * reciprocal_alphabet_size) >> 32) as i32
}

impl CompressedNeighborInfo {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:969-998
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1289-1300
pub fn compressed_lookup_parameters(word_size: usize) -> Option<(usize, i32)> {
    match word_size {
        5 | 6 => Some((15, 286_331_154)),
        7 => Some((10, 429_496_730)),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4949-4952
fn compressed_alphabet_string(compressed_alphabet_size: usize) -> Option<&'static [u8]> {
    match compressed_alphabet_size {
        10 => Some(b"IJLMV AST BDENZ KQR G FY P H C W"),
        15 => Some(b"ST IJV LM KR EQZ A G BD P N F Y H C W"),
        _ => None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4786-4805
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4786-4805
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
#[inline]
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4818-4858
// Always inlined: BLOSUM62's builder keeps the code it had as the only caller.
#[inline(always)]
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4886-4937
pub fn build_blosum62_compressed_score_matrix(
    compressed_alphabet_size: usize,
) -> Option<CompressedScoreMatrix> {
    compressed_score_matrix_of(ScoringMatrix::Blosum62, compressed_alphabet_size)
}

/// `build_blosum62_compressed_score_matrix` for the search matrix `matrix`: its ungapped lambda
/// (row 0 of NCBI's table, `RPSfindUngappedLambda`) and frequency ratios
/// (`_PSIMatrixFrequencyRatiosNew`). `None` for IDENTITY, which has no frequency ratios.
pub fn build_compressed_score_matrix(
    matrix: ScoringMatrix,
    compressed_alphabet_size: usize,
) -> Option<CompressedScoreMatrix> {
    compressed_score_matrix_of(matrix, compressed_alphabet_size)
}

// The body of both builders above, so that BLOSUM62's keeps its constant matrix.
#[inline(always)]
fn compressed_score_matrix_of(
    matrix: ScoringMatrix,
    compressed_alphabet_size: usize,
) -> Option<CompressedScoreMatrix> {
    let reverse = build_compressed_reverse_lookup(compressed_alphabet_size)?;
    let compressed_prob = compressed_probabilities(compressed_alphabet_size)?;
    let ungapped = lookup_protein_params_ungapped(matrix);
    let matrix_info = build_matrix_info(matrix, ungapped.lambda).ok()?;
    let matrix_scale_factor = 100.0 / ungapped.lambda;
    let min_freq = BLAST_SCORE_MIN as f64 / matrix_scale_factor;
    let mut scores = [[0i32; BLASTAA_SIZE]; BLASTAA_SIZE];

    for q in 0..BLASTAA_SIZE {
        for (s, group) in reverse.iter().enumerate() {
            let mut val = 0.0;
            for &aa in group {
                val += matrix_info.start_freq_ratios[q][aa] * compressed_prob[aa];
            }
            let val = if val < 1.0e-8 {
                min_freq
            } else {
                val.glibc_ln()
            };
            scores[q][s] = blast_nint(val * matrix_scale_factor);
        }
    }

    Some(scores)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1283-1355
// `Some(Err(()))` when the words fill NCBI's overflow banks (`s_CompressedListGetNewCell`),
// past which NCBI writes beyond `overflow_banks`.
#[cfg(test)]
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
#[cfg(test)]
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
/// matrix units (`opt->threshold` times `kMatrixScale`, cast to `Int4`, converted by the caller as the
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
    try_build_compressed_lookup_with_capacity(
        Blosum62Compressed,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        NCBI_OVERFLOW_CELL_CAPACITY,
    )
}

/// `try_build_blosum62_compressed_lookup` with the words of the search matrix `matrix`
/// (`build_compressed_score_matrix`). `None` also for IDENTITY, whose table NCBI cannot build.
pub fn try_build_compressed_lookup(
    matrix: ScoringMatrix,
    word_size: usize,
    scaled_threshold: i32,
    query: &[u8],
    locations: &[(i32, i32)],
    tables: &[std::ops::Range<usize>],
) -> Option<Result<BlastCompressedAaLookupTable, ()>> {
    try_build_compressed_lookup_with_capacity(
        matrix,
        word_size,
        scaled_threshold,
        query,
        locations,
        tables,
        NCBI_OVERFLOW_CELL_CAPACITY,
    )
}

/// The search matrix whose compressed score matrix scores a table's words: BLOSUM62
/// (`Blosum62Compressed`, `build_blosum62_compressed_score_matrix`) or another `ScoringMatrix`
/// (`build_compressed_score_matrix`).
pub trait CompressedScores: Copy {
    fn score_matrix(self, compressed_alphabet_size: usize) -> Option<CompressedScoreMatrix>;
}

/// BLOSUM62, the default matrix, as a type of its own: the builders' BLOSUM62 instances.
#[derive(Clone, Copy)]
pub struct Blosum62Compressed;

impl CompressedScores for Blosum62Compressed {
    #[inline(always)]
    fn score_matrix(self, compressed_alphabet_size: usize) -> Option<CompressedScoreMatrix> {
        build_blosum62_compressed_score_matrix(compressed_alphabet_size)
    }
}

impl CompressedScores for ScoringMatrix {
    fn score_matrix(self, compressed_alphabet_size: usize) -> Option<CompressedScoreMatrix> {
        build_compressed_score_matrix(self, compressed_alphabet_size)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-213
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:798-834
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
fn try_build_compressed_lookup_with_capacity<S: CompressedScores>(
    scores: S,
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
    let matrix = scores.score_matrix(lookup.compressed_alphabet_size as usize)?;
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
// One of NCBI's tables (a query batch's lookup segments), built with NCBI's capacity only to
// see whether its words fill the banks; it stops there. Its offsets are those of the whole
// query, which changes no count of a cell.
// Always inlined: BLOSUM62's builder keeps the code it had as the only caller.
#[inline(always)]
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1071-1121
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1019-1053
// No NCBI counterpart: counts the words of `add_word_hits` (above) without a table, for LOSAT's
// choice of a path and of the blocks' ends (EXPERIMENT LOSAT_X_BPQBLOCK, `blastp/blast_engine.rs`
// `x_block_search_refuses_nothing` and `x_word_blocks`); it changes no value NCBI computes.
/// The number of words that `add_word_hits` adds to a table for one query offset: none when a
/// letter of the word has no compressed letter; otherwise the word itself (`add_unencoded`) when
/// there is no threshold or its own score is below it, and with a threshold the neighbors that
/// `add_word_hits_core` enumerates (the same loops, counted). The enumeration bounds the
/// positions after the current one by their row maxima, which a word reaches, so its neighbors
/// are exactly the compressed words whose score reaches the threshold. Neither their number nor
/// the word's own score depends on the order of its letters, so a count is kept per set of
/// letters (`counts`, at the set's rank among the sets of `word_length` letters that have a
/// compressed letter) and the letters of a proteome take a few tens of thousands of enumerations.
/// The kept counts are shared: threads that count words with one counter keep them for each other
/// (two threads may count a set both; they keep the same count).
pub(crate) struct CompressedWordCount {
    threshold: i32,
    word_length: usize,
    compressed_alphabet_size: usize,
    compress_table: [u8; BLASTAA_SIZE],
    matrix: CompressedScoreMatrix,
    info: CompressedNeighborInfo,
    /// Each residue's index among the residues that have a compressed letter; `u8::MAX` for the
    /// others.
    letter_index: [u8; BLASTAA_SIZE],
    /// `binomial[n][k]` is C(n, k), for n below the residues with a compressed letter plus
    /// `word_length` and k up to `word_length`.
    binomial: Vec<[u32; X_MOST_WORD_LETTERS + 1]>,
    /// Per rank of a set of letters, its count plus 1; 0 while it is not counted.
    counts: Box<[std::sync::atomic::AtomicU32]>,
}

// No NCBI counterpart: the longest word of `compressed_lookup_parameters`.
const X_MOST_WORD_LETTERS: usize = 7;

impl CompressedWordCount {
    /// The counts of the table of `build_blosum62_compressed_lookup` with a threshold in NCBI's
    /// matrix units (as `try_build_blosum62_compressed_lookup` takes it).
    pub(crate) fn new(word_size: usize, scaled_threshold: i32) -> Option<Self> {
        let (compressed_alphabet_size, _) = compressed_lookup_parameters(word_size)?;
        if word_size > X_MOST_WORD_LETTERS {
            return None;
        }
        let matrix = build_blosum62_compressed_score_matrix(compressed_alphabet_size)?;
        let compress_table = build_compressed_translation(compressed_alphabet_size)?;
        let mut letter_index = [u8::MAX; BLASTAA_SIZE];
        let mut letters = 0usize;
        for (index, &compressed) in letter_index.iter_mut().zip(&compress_table) {
            if usize::from(compressed) < compressed_alphabet_size {
                *index = letters as u8;
                letters += 1;
            }
        }
        // A set of `word_size` of the `letters` letters, sorted, a_0 <= ... <= a_{w-1}, has the
        // rank sum C(a_i + i, i + 1) (the combinatorial number system of the increasing a_i + i,
        // below letters + w - 1), one of C(letters + w - 1, w).
        let mut binomial = vec![[0u32; X_MOST_WORD_LETTERS + 1]; letters + word_size];
        for n in 0..binomial.len() {
            binomial[n][0] = 1;
            for k in 1..=word_size.min(n) {
                binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
            }
        }
        let sets = binomial[letters + word_size - 1][word_size] as usize;
        Some(Self {
            threshold: scaled_threshold,
            word_length: word_size,
            compressed_alphabet_size,
            compress_table,
            info: CompressedNeighborInfo::new(&matrix, compressed_alphabet_size),
            matrix,
            letter_index,
            binomial,
            counts: (0..sets)
                .map(|_| std::sync::atomic::AtomicU32::new(0))
                .collect(),
        })
    }

    /// [`Self::new`] for the table of the search matrix `matrix` (`try_build_compressed_lookup`):
    /// the compressed scores of another matrix than BLOSUM62 (`build_compressed_score_matrix`) in
    /// place of BLOSUM62's. BLASTP's block search counts the words of its matrix (item SC2b).
    pub(crate) fn new_of(
        matrix: ScoringMatrix,
        word_size: usize,
        scaled_threshold: i32,
    ) -> Option<Self> {
        let mut count = Self::new(word_size, scaled_threshold)?;
        if matrix != ScoringMatrix::Blosum62 {
            count.matrix = build_compressed_score_matrix(matrix, count.compressed_alphabet_size)?;
            count.info = CompressedNeighborInfo::new(&count.matrix, count.compressed_alphabet_size);
        }
        Some(count)
    }

    /// The rank of the set of `word`'s letters; `None` when a letter has no compressed letter.
    fn rank(&self, word: &[u8]) -> Option<usize> {
        let mut sorted = [0u8; X_MOST_WORD_LETTERS];
        for (index, &residue) in sorted.iter_mut().zip(word) {
            *index = *self.letter_index.get(usize::from(residue))?;
            if *index == u8::MAX {
                return None;
            }
        }
        let sorted = &mut sorted[..self.word_length];
        sorted.sort_unstable();
        Some(
            sorted
                .iter()
                .enumerate()
                .map(|(i, &index)| self.binomial[usize::from(index) + i][i + 1] as usize)
                .sum(),
        )
    }

    /// The words that `add_word_hits` adds for `word` (NCBISTDAA letters, `word_length` of them).
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
        // At most 15^6 + 1 words (15 letters at word size 6; 15^5 + 1 at 5, 10^7 + 1 at 7).
        if let Some(words_and_one) = u32::try_from(words).ok().and_then(|w| w.checked_add(1)) {
            kept.store(words_and_one, Ordering::Relaxed);
        }
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

    /// The words of `queries`' offsets (each query's NCBISTDAA letters), the most that a table
    /// whose words are some of them holds; `None` once they are more than
    /// `overflow_cell_capacity` overflow cells can take. A backbone cell with n offsets takes
    /// cells(n) <= n / 4 of them (`try_build_compressed_lookup_with_capacity`), so such a table
    /// takes at most a quarter of the queries' words. Stops counting once they are more than four
    /// times the capacity.
    fn words_within(
        &self,
        queries: &mut dyn Iterator<Item = &[u8]>,
        overflow_cell_capacity: usize,
    ) -> Option<u64> {
        let limit = (COMPRESSED_HITS_PER_OVERFLOW_CELL as u64)
            .saturating_mul(overflow_cell_capacity as u64);
        let mut words = 0u64;
        for query in queries {
            for word in query.windows(self.word_length) {
                words = words.saturating_add(self.words(word));
                if words > limit {
                    return None;
                }
            }
        }
        Some(words)
    }

    /// The words of `queries`' offsets when no table whose words are some of them can fill NCBI's
    /// overflow banks: a query batch's table, or the table of one of its query chunks (whose
    /// parts are pieces of the batch's queries). `None` when their words are more than four
    /// times the banks' cells, though the table might still fit.
    pub(crate) fn batch_words(&self, queries: &mut dyn Iterator<Item = &[u8]>) -> Option<u64> {
        self.words_within(queries, NCBI_OVERFLOW_CELL_CAPACITY)
    }

    /// The cells of NCBI's backbone of the table that a word's index can reach: the compressed
    /// alphabet's size to the power of the word length. NCBI's `backbone_size` (blast_aalookup.c
    /// above) is one more, a last cell that no word's index reaches.
    pub(crate) fn backbone_cells(&self) -> u64 {
        (self.compressed_alphabet_size as u64).pow(self.word_length as u32)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1339-1347
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

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1337-1347
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:269-290
fn build_scaled_compress_table_u8(scaled_compress_table: &[i32; BLASTAA_SIZE]) -> [i32; 256] {
    let mut table = [-1i32; 256];
    table[..BLASTAA_SIZE].copy_from_slice(scaled_compress_table);
    table
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:859-909
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:304-325
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:304-325
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

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:49-57
#[inline(always)]
pub fn compressed_pv_test(pv: &[u32], index: usize, pv_array_bts: usize) -> bool {
    (pv[index >> pv_array_bts] & (1u32 << (index & PV_ARRAY_MASK))) != 0
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_lookup.h:49-57
#[inline(always)]
pub fn compressed_pv_set(pv: &mut [u32], index: usize, pv_array_bts: usize) {
    pv[index >> pv_array_bts] |= 1u32 << (index & PV_ARRAY_MASK);
}

#[cfg(test)]
#[path = "compressed_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "compressed_lookup_tests.rs"]
mod compressed_lookup_tests;
