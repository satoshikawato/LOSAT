//! The stages of PSSM creation: the port of NCBI's `algo/blast/core/blast_psi_priv.c`
//! (598d8ae6), in the order `PSICreatePssmWithDiagnostics` runs them (`mod.rs`). Integer fields
//! keep NCBI's types (`Uint4` sizes wrap as NCBI's do); the floating-point operations are NCBI's,
//! in NCBI's order, with the C library's `log`, `exp` and `pow`.

use super::{
    KarlinBlk, PsiDiagnostics, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiMsaCell, PsiScoreBlk,
};
use crate::stats::karlin_calc::{compute_karlin_params_ungapped, ScoreFreqProfile};
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:55-61
const K_PSI_NEAR_IDENTICAL: f64 = 0.94;
const K_PSI_IDENTICAL: f64 = 1.0;
const K_QUERY_INDEX: usize = 0;
const K_EPSILON: f64 = 0.0001;
const K_PSI_SCALE_FACTOR: i32 = 200;
const K_POSIT_SCALING_PERCENT: f64 = 0.05;
const K_POSIT_SCALING_NUM_ITERATIONS: u32 = 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:550
const EFFECTIVE_ALPHABET: usize = 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2065-2066
const MAX_IND_OBSERVATIONS: usize = 400;
const PSEUDO_MAX: f64 = 1000000.0;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/ncbi_math.h:161
const NCBIMATH_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:122
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

const GAP: u8 = ncbistdaa::GAP;
const X: u8 = ncbistdaa::X;
const STAR: u8 = ncbistdaa::STOP;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
fn blast_nint(x: f64) -> i64 {
    let x = x + if x >= 0.0 { 0.5 } else { -0.5 };
    x as i64
}

/// `_PSIPackedMsa`: the input rows and whether each takes part in the PSSM.
pub(super) struct PackedMsa {
    query_length: usize,
    /// `dimensions->num_seqs`: the aligned sequences, without the query.
    num_seqs: usize,
    data: Vec<Vec<PsiMsaCell>>,
    use_sequence: Vec<bool>,
}

impl PackedMsa {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:159-178
    pub(super) fn new(msa: &PsiMsa) -> Self {
        Self {
            query_length: msa.query_length,
            num_seqs: msa.rows.len() - 1,
            data: msa.rows.clone(),
            use_sequence: vec![true; msa.rows.len()],
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:209-224
    // The count includes the query, so `_PSIMsaNew`'s `num_seqs` is one more than the
    // sequences it keeps: its last row stays empty (`Msa::new`).
    fn number_of_aligned_seqs(&self) -> usize {
        self.use_sequence.iter().filter(|&&used| used).count()
    }
}

/// `_PSIMsaCell`: a cell with the extents of its aligned region.
#[derive(Clone, Copy, Default)]
struct MsaCell {
    letter: u8,
    is_aligned: bool,
    left: i32,
    right: i32,
}

/// `_PSIMsa`.
pub(super) struct Msa {
    pub(super) query_length: usize,
    /// `dimensions->num_seqs` (`_PSIPackedMsaGetNumberOfAlignedSeqs`); `cell` has one more row.
    pub(super) num_seqs: usize,
    cell: Vec<Vec<MsaCell>>,
    pub(super) query: Vec<u8>,
    residue_counts: Vec<[u32; BLASTAA_SIZE]>,
    num_matching_seqs: Vec<u32>,
}

impl Msa {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:325-385
    // `_PSIAllocateMatrix` allocates with `calloc`: the rows not copied hold zeros.
    pub(super) fn new(packed: &PackedMsa) -> Self {
        let query_length = packed.query_length;
        let num_seqs = packed.number_of_aligned_seqs();
        let mut cell = vec![vec![MsaCell::default(); query_length]; num_seqs + 1];
        let mut ss = 0;
        for (s, row) in packed.data.iter().enumerate() {
            if !packed.use_sequence[s] {
                continue;
            }
            for (p, packed_cell) in row.iter().enumerate() {
                cell[ss][p] = MsaCell {
                    letter: packed_cell.letter,
                    is_aligned: packed_cell.is_aligned,
                    left: -1,
                    right: query_length as i32,
                };
            }
            ss += 1;
        }
        let query = packed.data[K_QUERY_INDEX]
            .iter()
            .map(|cell| cell.letter)
            .collect();
        let mut msa = Self {
            query_length,
            num_seqs,
            cell,
            query,
            residue_counts: vec![[0; BLASTAA_SIZE]; query_length],
            num_matching_seqs: vec![0; query_length],
        };
        msa.update_position_counts();
        msa
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1011-1024
    fn update_position_counts(&mut self) {
        // Reset the data in case this function is called again after the initial counts (the
        // structure group's compatibility mode).
        self.num_matching_seqs.fill(0);
        for counts in &mut self.residue_counts {
            *counts = [0; BLASTAA_SIZE];
        }
        for row in &self.cell {
            for (p, pos) in row.iter().enumerate() {
                if pos.is_aligned {
                    let res = pos.letter as usize;
                    if res >= BLASTAA_SIZE {
                        continue;
                    }
                    self.residue_counts[p][res] += 1;
                    self.num_matching_seqs[p] += 1;
                }
            }
        }
    }

    /// `_PSIStructureGroupCustomization` of the `msa`: the query row is made of gaps that are not
    /// aligned (`ignore_msa_master`), and the counts are made again.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:799-808
    pub(super) fn structure_group_customization(&mut self) {
        for cell in &mut self.cell[K_QUERY_INDEX] {
            cell.letter = 0;
            cell.is_aligned = false;
        }
        self.update_position_counts();
    }
}

/// `_PSIAlignedBlock`.
pub(super) struct AlignedBlock {
    /// `pos_extnt[i].left`, `pos_extnt[i].right` (`SSeqRange`, `Int4`).
    pos_extnt: Vec<(i32, i32)>,
    /// `size` (`Uint4`).
    size: Vec<u32>,
}

impl AlignedBlock {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:511-528
    pub(super) fn new(query_length: usize) -> Self {
        Self {
            pos_extnt: vec![(-1, query_length as i32); query_length],
            size: vec![0; query_length],
        }
    }
}

/// `_PSISequenceWeights`, without the field that only the diagnostics psiblast does not ask for
/// keep (`independent_observations`).
pub(super) struct SequenceWeights {
    row_sigma: Vec<f64>,
    norm_seq_weights: Vec<f64>,
    sigma: Vec<f64>,
    match_weights: Vec<[f64; BLASTAA_SIZE]>,
    pub(super) std_prob: [f64; BLASTAA_SIZE],
    gapless_column_weights: Vec<f64>,
    pos_distinct_distrib: Vec<[i32; EFFECTIVE_ALPHABET + 1]>,
    pos_num_participating: Vec<i32>,
}

impl SequenceWeights {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:564-616
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1341-1346
    pub(super) fn new(num_seqs: usize, query_length: usize) -> Self {
        Self {
            row_sigma: vec![0.0; num_seqs + 1],
            norm_seq_weights: vec![0.0; num_seqs + 1],
            sigma: vec![0.0; query_length],
            match_weights: vec![[0.0; BLASTAA_SIZE]; query_length],
            std_prob: crate::stats::karlin_calc::compute_std_aa_composition(),
            gapless_column_weights: vec![0.0; query_length],
            pos_distinct_distrib: vec![[0; EFFECTIVE_ALPHABET + 1]; query_length + 1],
            pos_num_participating: vec![0; query_length + 1],
        }
    }
}

/// `_PSIInternalPssmData`.
pub(super) struct InternalPssm {
    pssm: Vec<[i32; BLASTAA_SIZE]>,
    scaled_pssm: Vec<[i32; BLASTAA_SIZE]>,
    freq_ratios: Vec<[f64; BLASTAA_SIZE]>,
    pseudocounts: Vec<f64>,
}

impl InternalPssm {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:434-460
    pub(super) fn new(query_length: usize) -> Self {
        Self {
            pssm: vec![[0; BLASTAA_SIZE]; query_length],
            scaled_pssm: vec![[0; BLASTAA_SIZE]; query_length],
            freq_ratios: vec![[0.0; BLASTAA_SIZE]; query_length],
            pseudocounts: vec![0.0; query_length],
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Purge of biased segments
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:947-988
pub(super) fn purge_biased_segments(msa: &mut PackedMsa) {
    for s in K_QUERY_INDEX + 1..msa.num_seqs + 1 {
        purge_similar_alignments(msa, K_QUERY_INDEX, s, K_PSI_IDENTICAL);
    }
    for i in 1..msa.num_seqs + 1 {
        let mut j = 1;
        while i + j < msa.num_seqs + 1 {
            purge_similar_alignments(msa, j, i + j, K_PSI_NEAR_IDENTICAL);
            j += 1;
        }
    }
}

/// `_EPSIPurgeFsmState`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum PurgeState {
    Counting,
    Resting,
}

/// `_PSIAlignmentTraits`.
#[derive(Default)]
struct AlignmentTraits {
    start: u32,
    effective_length: u32,
    n_x_residues: u32,
    n_identical: u32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1066-1072
fn reset_alignment_traits(traits: &mut AlignmentTraits, position: u32) {
    *traits = AlignmentTraits {
        start: position,
        ..AlignmentTraits::default()
    };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1076-1112
fn handle_neither_aligned(
    traits: &AlignmentTraits,
    state: &mut PurgeState,
    msa: &mut PackedMsa,
    seq_index: usize,
    max_percent_identity: f64,
) {
    if *state == PurgeState::Counting {
        if traits.effective_length > 0 {
            let percent_identity = traits.n_identical as f64 / traits.effective_length as f64;
            if percent_identity >= max_percent_identity {
                let align_stop = traits.start + traits.effective_length + traits.n_x_residues;
                purge_aligned_region(msa, seq_index, traits.start as usize, align_stop as usize);
            }
        }
        *state = PurgeState::Resting;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1185-1252
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1118-1183
// The cells are read as they are when each position is reached: a purge rewrites the cells of
// `seq_index2` that the loop has passed.
fn purge_similar_alignments(
    msa: &mut PackedMsa,
    seq_index1: usize,
    seq_index2: usize,
    max_percent_identity: f64,
) {
    if seq_index1 == seq_index2 || !msa.use_sequence[seq_index1] || !msa.use_sequence[seq_index2] {
        return;
    }
    let mut state = PurgeState::Counting;
    let mut traits = AlignmentTraits::default();
    reset_alignment_traits(&mut traits, 0);
    for p in 0..msa.query_length {
        let seq1 = msa.data[seq_index1][p];
        let seq2 = msa.data[seq_index2][p];
        let pos1_aligned = if seq_index1 == K_QUERY_INDEX {
            false
        } else {
            seq1.is_aligned
        };
        let pos2_aligned = seq2.is_aligned;
        if !pos1_aligned && !pos2_aligned {
            handle_neither_aligned(&traits, &mut state, msa, seq_index2, max_percent_identity);
        } else {
            let neither_is_x = seq1.letter != X && seq2.letter != X;
            if neither_is_x {
                match state {
                    PurgeState::Counting => traits.effective_length += 1,
                    PurgeState::Resting => {
                        reset_alignment_traits(&mut traits, p as u32);
                        traits.effective_length = 1;
                        state = PurgeState::Counting;
                    }
                }
            } else if state == PurgeState::Counting {
                traits.n_x_residues += 1;
            }
            if neither_is_x
                && (pos2_aligned && seq1.is_aligned)
                && seq1.letter == seq2.letter
                && state == PurgeState::Counting
            {
                traits.n_identical += 1;
            }
        }
    }
    handle_neither_aligned(&traits, &mut state, msa, seq_index2, max_percent_identity);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2762-2806
fn purge_aligned_region(msa: &mut PackedMsa, seq_index: usize, start: usize, stop: usize) {
    if seq_index == 0 || seq_index > msa.num_seqs + 1 || stop > msa.query_length {
        return;
    }
    let row = &mut msa.data[seq_index];
    for cell in row.iter_mut().take(stop).skip(start) {
        *cell = PsiMsaCell {
            letter: 0,
            is_aligned: false,
        };
    }
    if !row.iter().any(|cell| cell.is_aligned) {
        msa.use_sequence[seq_index] = false;
    }
}

// ---------------------------------------------------------------------------------------------
// Validation
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:827-859
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:703-797
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:683-696
///
/// `ignore_unaligned_positions` (`PSIBlastOptions`, set by `CPsiBlastInputClustalW`) skips the
/// check of the aligned columns:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:836-848
pub(super) fn validate_msa(msa: &Msa, ignore_unaligned_positions: bool) -> Result<(), PsiError> {
    for row in &msa.cell {
        if let Some(cell) = row.iter().find(|cell| cell.is_aligned) {
            if cell.letter == GAP {
                return Err(PsiError::StartingGap);
            }
        }
    }
    for row in &msa.cell {
        if let Some(cell) = row.iter().rev().find(|cell| cell.is_aligned) {
            if cell.letter == GAP {
                return Err(PsiError::EndingGap);
            }
        }
    }
    for p in 0..(if ignore_unaligned_positions {
        0
    } else {
        msa.query_length
    }) {
        let mut found_aligned_sequence = false;
        let mut found_non_gap_residue = false;
        for row in &msa.cell {
            if row[p].is_aligned {
                found_aligned_sequence = true;
                if row[p].letter != GAP {
                    found_non_gap_residue = true;
                    break;
                }
            }
        }
        if !found_aligned_sequence {
            return Err(PsiError::UnalignedColumn);
        }
        if !found_non_gap_residue {
            return Err(PsiError::ColumnOfGaps);
        }
    }
    for p in 0..msa.query_length {
        if msa.cell[K_QUERY_INDEX][p].letter == GAP || msa.query[p] == GAP {
            return Err(PsiError::GapInQuery);
        }
    }
    if msa.num_seqs == 0 {
        return Err(PsiError::NoAlignedSeqs);
    }
    Ok(())
}

/// `_PSIValidateMSA_StructureGroup` of the `msa`: only the check that sequences take part
/// (`s_PSIValidateParticipatingSequences`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:810-819
pub(super) fn validate_msa_structure_group(msa: &Msa) -> Result<(), PsiError> {
    if msa.num_seqs == 0 {
        return Err(PsiError::NoAlignedSeqs);
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// Alignment blocks
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1306-1313
pub(super) fn compute_alignment_blocks(msa: &mut Msa, aligned_blocks: &mut AlignedBlock) {
    for s in K_QUERY_INDEX + 1..msa.num_seqs + 1 {
        get_left_extents(&mut msa.cell[s]);
        get_right_extents(&mut msa.cell[s]);
        compute_position_extents(&msa.cell[s], aligned_blocks);
    }
    compute_aligned_region_lengths(msa, aligned_blocks);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1330-1349
fn get_left_extents(row: &mut [MsaCell]) {
    if row.is_empty() {
        return;
    }
    if row[0].is_aligned && row[0].letter != GAP {
        row[0].left = 0;
    }
    for curr in 1..row.len() {
        if !row[curr].is_aligned {
            continue;
        }
        row[curr].left = if row[curr - 1].is_aligned {
            row[curr - 1].left
        } else {
            curr as i32
        };
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1364-1383
fn get_right_extents(row: &mut [MsaCell]) {
    let Some(last) = row.len().checked_sub(1) else {
        return;
    };
    if row[last].is_aligned && row[last].letter != GAP {
        row[last].right = last as i32;
    }
    for curr in (0..last).rev() {
        if !row[curr].is_aligned {
            continue;
        }
        row[curr].right = if row[curr + 1].is_aligned {
            row[curr + 1].right
        } else {
            curr as i32
        };
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1401-1417
// `PSI_IGNORE_GAPS_IN_COLUMNS` is not defined (blast_psi_priv.c:50).
fn compute_position_extents(row: &[MsaCell], aligned_blocks: &mut AlignedBlock) {
    for (i, cell) in row.iter().enumerate() {
        if cell.is_aligned {
            let extent = &mut aligned_blocks.pos_extnt[i];
            extent.0 = extent.0.max(cell.left);
            extent.1 = extent.1.min(cell.right);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1431-1459
fn compute_aligned_region_lengths(msa: &Msa, aligned_blocks: &mut AlignedBlock) {
    let query_length = msa.query_length;
    for i in 0..query_length {
        let (left, right) = aligned_blocks.pos_extnt[i];
        aligned_blocks.size[i] = right.wrapping_sub(left).wrapping_add(1) as u32;
    }
    for i in 0..query_length {
        if msa.query[i] != X {
            continue;
        }
        for idx in 0..i {
            if aligned_blocks.pos_extnt[idx].1 as u32 >= i as u32 && msa.query[idx] != X {
                aligned_blocks.size[idx] = aligned_blocks.size[idx].wrapping_sub(1);
            }
        }
        for idx in (i + 1..query_length).rev() {
            if aligned_blocks.pos_extnt[idx].0 as u32 <= i as u32 && msa.query[idx] != X {
                aligned_blocks.size[idx] = aligned_blocks.size[idx].wrapping_sub(1);
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Sequence weights
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1580-1639
// `last_calc_pos` is set to `pos` just before the test, so it always differs from one less than
// `pos`: the test always holds and the weights are computed at every position: the `else` branch is
// never taken. Without the structure group compatibility mode, `kExpectedNumMatchingSeqs` is 1.
pub(super) fn compute_sequence_weights(
    msa: &Msa,
    aligned_blocks: &AlignedBlock,
    nsg_compatibility_mode: bool,
    seq_weights: &mut SequenceWeights,
) -> Result<(), PsiError> {
    let expected_num_matching_seqs = expected_matching_seqs(nsg_compatibility_mode);
    let mut aligned_seqs: Vec<usize> = Vec::with_capacity(msa.num_seqs + 1);
    for pos in 0..msa.query_length {
        if aligned_blocks.size[pos] == 0 || msa.num_matching_seqs[pos] <= expected_num_matching_seqs
        {
            continue;
        }
        get_aligned_sequences_for_position(msa, pos, &mut aligned_seqs);
        if aligned_seqs.len() as u32 <= expected_num_matching_seqs {
            continue;
        }
        seq_weights.norm_seq_weights.fill(0.0);
        seq_weights.row_sigma.fill(0.0);
        calculate_normalized_sequence_weights(msa, aligned_blocks, pos, &aligned_seqs, seq_weights);
        seq_weights.pos_num_participating[pos] = aligned_seqs.len() as i32;
        calculate_match_weights(msa, pos, &aligned_seqs, seq_weights);
    }
    check_sequence_weights(msa, seq_weights, nsg_compatibility_mode)?;
    spread_gap_weights(msa, seq_weights, nsg_compatibility_mode);
    check_sequence_weights(msa, seq_weights, nsg_compatibility_mode)
}

/// `kExpectedNumMatchingSeqs` of `_PSIComputeSequenceWeights`, `_PSISpreadGapWeights` and
/// `_PSICheckSequenceWeights`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1566
fn expected_matching_seqs(nsg_compatibility_mode: bool) -> u32 {
    if nsg_compatibility_mode {
        0
    } else {
        1
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1915-1928
fn get_aligned_sequences_for_position(msa: &Msa, position: usize, aligned: &mut Vec<usize>) {
    aligned.clear();
    for (i, row) in msa.cell.iter().enumerate() {
        if row[position].is_aligned {
            aligned.push(i);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1785-1871
fn calculate_normalized_sequence_weights(
    msa: &Msa,
    aligned_blocks: &AlignedBlock,
    position: usize,
    aligned_seqs: &[usize],
    seq_weights: &mut SequenceWeights,
) {
    let (left, right) = aligned_blocks.pos_extnt[position];
    let mut distinct_residues_found = false;
    let mut sigma: u32 = 0;
    let mut i = left as u32;
    while i <= right as u32 {
        let column = i as usize;
        let mut residue_counts_for_column = [0u32; BLASTAA_SIZE];
        let mut num_distinct_residues_for_column: u32 = 0;
        let mut num_local_std_letters: usize = 0;
        for &seq_idx in aligned_seqs {
            let residue = msa.cell[seq_idx][column].letter;
            if residue_counts_for_column[residue as usize] == 0 {
                num_distinct_residues_for_column += 1;
                if residue != GAP && residue != X {
                    num_local_std_letters += 1;
                }
            }
            residue_counts_for_column[residue as usize] += 1;
        }
        sigma = sigma.wrapping_add(num_distinct_residues_for_column);
        let num_local_std_letters = num_local_std_letters.min(EFFECTIVE_ALPHABET);
        seq_weights.pos_distinct_distrib[position][num_local_std_letters] += 1;
        if num_distinct_residues_for_column > 1 {
            distinct_residues_found = true;
        }
        for &seq_idx in aligned_seqs {
            let residue = msa.cell[seq_idx][column].letter;
            seq_weights.row_sigma[seq_idx] += 1.0
                / (residue_counts_for_column[residue as usize] * num_distinct_residues_for_column)
                    as f64;
        }
        i += 1;
    }
    seq_weights.sigma[position] = sigma as f64;
    if distinct_residues_found {
        let mut weight_sum = 0.0;
        let extent = (right - left + 1) as f64;
        for &seq_idx in aligned_seqs {
            seq_weights.norm_seq_weights[seq_idx] = seq_weights.row_sigma[seq_idx] / extent;
            weight_sum += seq_weights.norm_seq_weights[seq_idx];
        }
        for &seq_idx in aligned_seqs {
            seq_weights.norm_seq_weights[seq_idx] /= weight_sum;
        }
    } else {
        for &seq_idx in aligned_seqs {
            seq_weights.norm_seq_weights[seq_idx] = 1.0 / aligned_seqs.len() as f64;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1888-1900
fn calculate_match_weights(
    msa: &Msa,
    position: usize,
    aligned_seqs: &[usize],
    seq_weights: &mut SequenceWeights,
) {
    for &seq_idx in aligned_seqs {
        let residue = msa.cell[seq_idx][position].letter;
        seq_weights.match_weights[position][residue as usize] +=
            seq_weights.norm_seq_weights[seq_idx];
        if residue != GAP {
            seq_weights.gapless_column_weights[position] += seq_weights.norm_seq_weights[seq_idx];
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1945-1963
fn spread_gap_weights(msa: &Msa, seq_weights: &mut SequenceWeights, nsg_compatibility_mode: bool) {
    let expected_num_matching_seqs = expected_matching_seqs(nsg_compatibility_mode);
    for pos in 0..msa.query_length {
        if msa.num_matching_seqs[pos] <= expected_num_matching_seqs
            || msa.cell[K_QUERY_INDEX][pos].letter == X
        {
            continue;
        }
        for res in 0..BLASTAA_SIZE {
            if seq_weights.std_prob[res] > K_EPSILON {
                seq_weights.match_weights[pos][res] +=
                    seq_weights.match_weights[pos][GAP as usize] * seq_weights.std_prob[res];
            }
        }
        seq_weights.match_weights[pos][GAP as usize] = 0.0;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1990-2016
// `!(0.99..=1.01).contains(..)` is NCBI's two comparisons for every number but NaN, which
// NCBI's lets through and this rejects. The weights are finite (sums of finite products,
// normalised by a positive total in `_PSICalculateNormalizedSequenceWeights`), so no NaN
// reaches this check.
fn check_sequence_weights(
    msa: &Msa,
    seq_weights: &SequenceWeights,
    nsg_compatibility_mode: bool,
) -> Result<(), PsiError> {
    let expected_num_matching_seqs = expected_matching_seqs(nsg_compatibility_mode);
    for pos in 0..msa.query_length {
        if msa.num_matching_seqs[pos] <= expected_num_matching_seqs
            || msa.cell[K_QUERY_INDEX][pos].letter == X
        {
            continue;
        }
        let mut running_total = 0.0;
        for weight in &seq_weights.match_weights[pos] {
            running_total += *weight;
        }
        if !(0.99..=1.01).contains(&running_total) {
            return Err(PsiError::BadSeqWeights);
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------------------------
// Frequency ratios
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2093-2180
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2084-2085
pub(super) fn compute_freq_ratios(
    msa: &Msa,
    seq_weights: &mut SequenceWeights,
    sbp: &PsiScoreBlk,
    aligned_blocks: &AlignedBlock,
    pseudo_count: i32,
    nsg_compatibility_mode: bool,
    internal: &mut InternalPssm,
) -> Result<(), PsiError> {
    const K_ZERO_OBS_PSEUDO: f64 = 30.0;
    let expno = initialize_exp_num_observations(sbp.background);
    for p in 0..msa.query_length {
        let query_is_x = msa.cell[K_QUERY_INDEX][p].letter == X;
        let mut column_counts = 0.0;
        let mut observations = 0.0;
        if !query_is_x {
            observations = effective_observations(
                aligned_blocks,
                seq_weights,
                p,
                msa.query_length as i32,
                &expno,
            );
            column_counts = if pseudo_count == 0 {
                column_specific_pseudocounts(seq_weights, p, sbp.background, observations)
            } else {
                pseudo_count as f64
            };
        }
        let pseudo_weight = if column_counts >= PSEUDO_MAX {
            observations = 0.0;
            K_ZERO_OBS_PSEUDO
        } else {
            column_counts
        };
        for r in 0..BLASTAA_SIZE {
            if query_is_x || seq_weights.std_prob[r] <= K_EPSILON {
                internal.freq_ratios[p][r] = 0.0;
                continue;
            }
            let beta = pseudo_weight;
            internal.pseudocounts[p] = beta;
            let mut pseudo = 0.0;
            for i in 0..BLASTAA_SIZE {
                if sbp.matrix[r][i] != BLAST_SCORE_MIN {
                    pseudo += seq_weights.match_weights[p][i] * sbp.freq_ratios[r][i];
                }
            }
            pseudo *= beta;
            let numerator =
                observations * seq_weights.match_weights[p][r] / seq_weights.std_prob[r] + pseudo;
            let denominator = observations + beta;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2162-2166
            if nsg_compatibility_mode && denominator == 0.0 {
                return Err(PsiError::Unknown);
            }
            let q_over_p_estimate = numerator / denominator;
            internal.freq_ratios[p][r] = q_over_p_estimate * seq_weights.std_prob[r];
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3068-3084
// `expno[MAX_IND_OBSERVATIONS]` is never written nor read (`s_effectiveObservations` stops at
// `MAX_IND_OBSERVATIONS - 1`).
fn initialize_exp_num_observations(background: &[f64; 20]) -> [f64; MAX_IND_OBSERVATIONS] {
    let mut expno = [0.0; MAX_IND_OBSERVATIONS];
    for (j, value) in expno.iter_mut().enumerate().skip(1) {
        let mut weighted_sum = 0.0;
        for &probability in background.iter().take(EFFECTIVE_ALPHABET) {
            weighted_sum += (j as f64 * (1.0 - probability).glibc_ln()).glibc_exp();
        }
        *value = EFFECTIVE_ALPHABET as f64 - weighted_sum;
    }
    expno
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3135-3175
// The distribution of a column whose extent was set holds one count per column of the extent
// (`_PSICalculateNormalizedSequenceWeights`), at least `halfNumColumns`, so `k` stays at 0 or
// above. The loop's `k >= 0` stands for NCBI's `ASSERT` of the same test, which only debug builds
// check: a release build of NCBI would read before the row there; LOSAT stops instead (not
// reached).
fn effective_observations(
    align_blk: &AlignedBlock,
    seq_weights: &SequenceWeights,
    column: usize,
    query_length: i32,
    expno: &[f64; MAX_IND_OBSERVATIONS],
) -> f64 {
    let (left, right) = align_blk.pos_extnt[column];
    if left < 0 {
        return 0.0;
    }
    if right >= query_length {
        return 0.0;
    }
    let half_num_columns = 1.max((right - left + 2) / 2);
    let mut k = EFFECTIVE_ALPHABET as i32;
    let mut columns_accounted_for: i32 = 0;
    let mut total_distinct_counts: i32 = 0;
    while columns_accounted_for < half_num_columns && k >= 0 {
        let count = seq_weights.pos_distinct_distrib[column][k as usize];
        total_distinct_counts += count * k;
        columns_accounted_for += count;
        if columns_accounted_for > half_num_columns {
            total_distinct_counts -= (columns_accounted_for - half_num_columns) * k;
            columns_accounted_for = half_num_columns;
        }
        k -= 1;
    }
    let ave_distinct_aa = total_distinct_counts as f64 / columns_accounted_for as f64;
    let mut i = 1;
    while i < MAX_IND_OBSERVATIONS && expno[i] <= ave_distinct_aa {
        i += 1;
    }
    let mut indep = if i == MAX_IND_OBSERVATIONS {
        i as f64
    } else {
        i as f64 - (expno[i] - ave_distinct_aa) / (expno[i] - expno[i - 1])
    };
    let participating = seq_weights.pos_num_participating[column] as f64;
    // MIN(a, b) and MAX(a, b) as NCBI's macros compare.
    indep = if indep < participating {
        indep
    } else {
        participating
    };
    if 0.0 > indep - 1.0 {
        0.0
    } else {
        indep - 1.0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2972-3004
const CHAR_ORDER: [usize; EFFECTIVE_ALPHABET] = [
    1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3086-3118
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3014-3060
fn column_specific_pseudocounts(
    seq_weights: &SequenceWeights,
    column: usize,
    background: &[f64; 20],
    observations: f64,
) -> f64 {
    const K_PSEUDO_MULT: f64 = 500.0;
    const K_PSEUDO_NUMERATOR: f64 = 0.0457;
    const K_PSEUDO_EXPONENT: f64 = 0.8;
    const K_PSEUDO_SMALL_INITIAL: f64 = 5.5;
    const K_POS_EPSILON: f64 = 0.0001;
    let mut intermediate_sums = [0.0; EFFECTIVE_ALPHABET];
    let mut overall_sum = 0.0;
    for c in 0..EFFECTIVE_ALPHABET {
        let initial = seq_weights.match_weights[column][CHAR_ORDER[c]];
        intermediate_sums[c] = initial * observations + background[c] * K_PSEUDO_SMALL_INITIAL;
        overall_sum += intermediate_sums[c];
    }
    let mut relative_entropy = 0.0;
    for c in 0..EFFECTIVE_ALPHABET {
        let adjusted = intermediate_sums[c] / overall_sum;
        if adjusted > K_POS_EPSILON {
            relative_entropy += adjusted * (adjusted / background[c]).glibc_ln();
        }
    }
    if relative_entropy < K_POS_EPSILON {
        relative_entropy = K_POS_EPSILON;
    }
    let pseudo_denominator = relative_entropy.glibc_pow(K_PSEUDO_EXPONENT);
    let alpha = K_PSEUDO_NUMERATOR / pseudo_denominator;
    if alpha < 1.0 - K_POS_EPSILON {
        K_PSEUDO_MULT * alpha / (1.0 - alpha)
    } else {
        PSEUDO_MAX
    }
}

// ---------------------------------------------------------------------------------------------
// Scores
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2409-2469
pub(super) fn convert_freq_ratios_to_pssm(
    internal: &mut InternalPssm,
    query: &[u8],
    sbp: &PsiScoreBlk,
    std_probs: &[f64; BLASTAA_SIZE],
) {
    let ideal_lambda = sbp.ideal_lambda;
    for (i, &residue) in query.iter().enumerate() {
        let k_residue = residue as usize;
        let mut is_unaligned_column = true;
        for j in 0..BLASTAA_SIZE {
            let mut q_over_p_estimate = 0.0;
            if std_probs[j] > K_EPSILON {
                q_over_p_estimate = internal.freq_ratios[i][j] / std_probs[j];
            }
            if is_unaligned_column && q_over_p_estimate != 0.0 {
                is_unaligned_column = false;
            }
            internal.scaled_pssm[i][j] = if q_over_p_estimate == 0.0 || std_probs[j] < K_EPSILON {
                BLAST_SCORE_MIN
            } else {
                let tmp = q_over_p_estimate.glibc_ln() / ideal_lambda;
                blast_nint(K_PSI_SCALE_FACTOR as f64 * tmp) as i32
            };
            if (j == X as usize || j == STAR as usize)
                && sbp.matrix[k_residue][X as usize] != BLAST_SCORE_MIN
            {
                internal.scaled_pssm[i][j] = sbp.matrix[k_residue][j] * K_PSI_SCALE_FACTOR;
            }
        }
        if is_unaligned_column {
            for j in 0..BLASTAA_SIZE {
                internal.pssm[i][j] = sbp.matrix[k_residue][j];
                let ratio = sbp.freq_ratios[k_residue][j];
                internal.scaled_pssm[i][j] = if ratio != 0.0 {
                    let tmp = (K_PSI_SCALE_FACTOR * sbp.bit_scale_factor) as f64 * ratio.glibc_ln()
                        / NCBIMATH_LN2;
                    blast_nint(tmp) as i32
                } else {
                    BLAST_SCORE_MIN
                };
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2509-2593
/// Returns `kbp_psi[0]` after the scaling: the block of its last `_PSIUpdateLambdaK`.
pub(super) fn scale_matrix(
    query: &[u8],
    std_probs: &[f64; BLASTAA_SIZE],
    internal: &mut InternalPssm,
    sbp: &PsiScoreBlk,
) -> Result<KarlinBlk, PsiError> {
    let ideal_lambda = sbp.ideal_lambda;
    let mut first_time = true;
    let mut factor = 1.0;
    let mut factor_low = 1.0;
    let mut factor_high = 1.0;
    let mut too_high = true;
    let mut kbp_psi;
    loop {
        rescale(internal, factor);
        kbp_psi = update_lambda_k(&internal.pssm, query, std_probs);
        let new_lambda = kbp_psi.lambda;
        if new_lambda > ideal_lambda {
            if first_time {
                factor_high = 1.0 + K_POSIT_SCALING_PERCENT;
                factor = factor_high;
                factor_low = 1.0;
                too_high = true;
                first_time = false;
            } else {
                if !too_high {
                    break;
                }
                factor_high += factor_high - 1.0;
                factor = factor_high;
            }
        } else if new_lambda > 0.0 {
            if first_time {
                factor_high = 1.0;
                factor_low = 1.0 - K_POSIT_SCALING_PERCENT;
                factor = factor_low;
                too_high = false;
                first_time = false;
            } else {
                if too_high {
                    break;
                }
                factor_low += factor_low - 1.0;
                factor = factor_low;
            }
        } else {
            return Err(PsiError::PositiveAvgScore);
        }
    }
    for _ in 0..K_POSIT_SCALING_NUM_ITERATIONS {
        factor = (factor_high + factor_low) / 2.0;
        rescale(internal, factor);
        kbp_psi = update_lambda_k(&internal.pssm, query, std_probs);
        let new_lambda = kbp_psi.lambda;
        if new_lambda > ideal_lambda {
            factor_low = factor;
        } else {
            factor_high = factor;
        }
    }
    Ok(kbp_psi)
}

/// The rescaling of `_PSIScaleMatrix`'s loops (cited above `scale_matrix`).
fn rescale(internal: &mut InternalPssm, factor: f64) {
    for (pssm, scaled) in internal.pssm.iter_mut().zip(&internal.scaled_pssm) {
        for (score, &scaled) in pssm.iter_mut().zip(scaled) {
            *score = if scaled != BLAST_SCORE_MIN {
                blast_nint(factor * scaled as f64 / K_PSI_SCALE_FACTOR as f64) as i32
            } else {
                BLAST_SCORE_MIN
            };
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2732-2753
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2699-2734
/// `kbp_psi[0]` after `_PSIUpdateLambdaK`: Lambda, K and H, all -1 when any of them cannot be
/// computed. The gapped PSI K that the call also sets depends on this block only
/// (`gapped_psi_k`), so it is computed once, from the last block.
fn update_lambda_k(
    pssm: &[[i32; BLASTAA_SIZE]],
    query: &[u8],
    std_probs: &[f64; BLASTAA_SIZE],
) -> KarlinBlk {
    let sfp = compute_score_probabilities(pssm, query, std_probs);
    match compute_karlin_params_ungapped(&sfp) {
        Ok(params) => KarlinBlk {
            lambda: params.lambda,
            k: params.k,
            h: params.h,
        },
        Err(_) => KarlinBlk {
            lambda: -1.0,
            k: -1.0,
            h: -1.0,
        },
    }
}

/// `sbp->kbp_gap_psi[0]->K` after `_PSIUpdateLambdaK` (cited above `update_lambda_k`).
pub(super) fn gapped_psi_k(kbp_psi: KarlinBlk, sbp: &PsiScoreBlk) -> f64 {
    kbp_psi.k * sbp.gap_std_k / sbp.ideal_k
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2100-2119
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2699-2704
/// Whether `Blast_ScoreFreqNew` makes the score frequencies of `pssm` in
/// `_PSIComputeScoreProbabilities` (the range of the 20 standard letters' scores over the
/// positions that are not X, without `BLAST_SCORE_MIN` and `BLAST_SCORE_MAX`, passes
/// `BlastScoreChk`: a negative smallest score and a positive largest one). Where it does not,
/// `Blast_KarlinBlkUngappedCalc` returns 1 and leaves `kbp_psi[0]` as it was.
fn factor_one_score_range_is_valid(pssm: &[[i32; BLASTAA_SIZE]], query: &[u8]) -> bool {
    let (mut min_score, mut max_score) = (BLAST_SCORE_MAX, BLAST_SCORE_MIN);
    for (_, column) in query.iter().zip(pssm).filter(|(&residue, _)| residue != X) {
        for &r in &STD_ALPHABET {
            let score = column[r as usize];
            if score > BLAST_SCORE_MIN && score < BLAST_SCORE_MAX {
                min_score = min_score.min(score);
                max_score = max_score.max(score);
            }
        }
    }
    // Inside (BLAST_SCORE_MIN, BLAST_SCORE_MAX) the range is below BLAST_SCORE_RANGE_MAX.
    min_score < 0 && max_score > 0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:355-378
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi.c:399-414
/// `PSICreatePssmFromFrequencyRatios` without IMPALA's scaling (`CPsiBlastInputFreqRatios`'s
/// default): the scores of the PSSM of `freq_ratios` (28 per position) and `kbp_psi[0]` of its
/// scaling.
pub(super) fn pssm_from_freq_ratios(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
) -> Result<(Vec<[i32; BLASTAA_SIZE]>, KarlinBlk), PsiError> {
    let std_probs = crate::stats::karlin_calc::compute_std_aa_composition();
    let mut internal = InternalPssm::new(query.len());
    internal.freq_ratios.copy_from_slice(freq_ratios);
    convert_freq_ratios_to_pssm(&mut internal, query, sbp, &std_probs);
    let ungapped = scale_matrix(query, &std_probs, &mut internal, sbp)?;
    Ok((internal.into_scores(), ungapped))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:715-743
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_posit.c:357-386
/// `s_ScalePosMatrix` of the position-based composition-based statistics: the scaled PSSM
/// (`posPrivateMatrix`, the redo's `startMatrix`) and `kbp_psi[0]` of the factor-1 PSSM
/// (`posMatrix`), or `None` where `_PSIUpdateLambdaK` leaves `kbp_psi[0]` as it was
/// (`factor_one_score_range_is_valid`). `posMatrix` and `posPrivateMatrix` start zeroed
/// (`_PSIAllocateMatrix`, `Nlm_Int4MatrixNew`); `impalaScaleMatrix` runs without its binary
/// search, so `factor` stays 1.0, and leaves the cells of `BLAST_SCORE_MIN` as the conversion left
/// them.
pub(super) fn scale_pos_matrix(
    query: &[u8],
    freq_ratios: &[[f64; BLASTAA_SIZE]],
    sbp: &PsiScoreBlk,
    scaling_factor: f64,
) -> (Vec<[i32; BLASTAA_SIZE]>, Option<KarlinBlk>) {
    let std_probs = crate::stats::karlin_calc::compute_std_aa_composition();
    let mut internal = InternalPssm::new(query.len());
    internal.freq_ratios.copy_from_slice(freq_ratios);
    convert_freq_ratios_to_pssm(&mut internal, query, sbp, &std_probs);
    let factor = 1.0;
    for (pssm, scaled) in internal.pssm.iter_mut().zip(&internal.scaled_pssm) {
        for (score, &scaled) in pssm.iter_mut().zip(scaled) {
            if scaled != BLAST_SCORE_MIN {
                *score = blast_nint(scaled as f64 * factor / K_PSI_SCALE_FACTOR as f64) as i32;
            }
        }
    }
    let kbp_psi = factor_one_score_range_is_valid(&internal.pssm, query)
        .then(|| update_lambda_k(&internal.pssm, query, &std_probs));
    let scalefactor = scaling_factor / K_PSI_SCALE_FACTOR as f64;
    for row in internal.scaled_pssm.iter_mut() {
        for score in row.iter_mut() {
            if *score != BLAST_SCORE_MIN {
                *score = blast_nint(*score as f64 * factor * scalefactor) as i32;
            }
        }
    }
    (internal.scaled_pssm, kbp_psi)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1795-1815
// `Blast_GetStdAlphabet` (blast_stat.c:1863-1883) gives these letters in NCBIstdaa, in this order.
const STD_ALPHABET: [u8; 20] = [
    ncbistdaa::A,
    ncbistdaa::C,
    ncbistdaa::D,
    ncbistdaa::E,
    ncbistdaa::F,
    ncbistdaa::G,
    ncbistdaa::H,
    ncbistdaa::I,
    ncbistdaa::K,
    ncbistdaa::L,
    ncbistdaa::M,
    ncbistdaa::N,
    ncbistdaa::P,
    ncbistdaa::Q,
    ncbistdaa::R,
    ncbistdaa::S,
    ncbistdaa::T,
    ncbistdaa::V,
    ncbistdaa::W,
    ncbistdaa::Y,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2675-2727
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2102-2104, 2118-2119 and
// 2703-2704
// NCBI's `Blast_ScoreFreqNew` returns NULL when the PSSM's smallest score is 0 or more or its
// largest 0 or less; `Blast_KarlinBlkUngappedCalc` then leaves `kbp_psi[0]` as it was (the
// previous Lambda). The scaled scores of a BLOSUM62 round (the only matrix LOSAT's psiblast
// accepts) keep negative and positive values at every factor the scaling tries, so a round does
// not reach that branch; a PSSM file's own ratios can (every ratio at most 1e-4 and not 0): the
// position-based redo models it (`scale_pos_matrix`, `factor_one_score_range_is_valid`, item
// PS4b), and for a PSSM of ratios only NCBI's bracket search of `_PSIScaleMatrix` then never ends
// (LOSAT refuses that PSSM, `pssm_from_freq_ratios`'s error).
fn compute_score_probabilities(
    pssm: &[[i32; BLASTAA_SIZE]],
    query: &[u8],
    std_probs: &[f64; BLASTAA_SIZE],
) -> ScoreFreqProfile {
    let effective_length = query.iter().filter(|&&residue| residue != X).count() as u32;
    let mut min_score = BLAST_SCORE_MAX;
    let mut max_score = BLAST_SCORE_MIN;
    let scores = || {
        query
            .iter()
            .zip(pssm)
            .filter(|(&residue, _)| residue != X)
            .flat_map(|(_, column)| STD_ALPHABET.iter().map(move |&r| (r, column[r as usize])))
            .filter(|&(_, score)| score > BLAST_SCORE_MIN && score < BLAST_SCORE_MAX)
    };
    for (_, score) in scores() {
        max_score = max_score.max(score);
        min_score = min_score.min(score);
    }
    let mut sprob = vec![0.0; (max_score - min_score + 1).max(0) as usize];
    for (residue, score) in scores() {
        sprob[(score - min_score) as usize] +=
            std_probs[residue as usize] / effective_length as f64;
    }
    let mut score_avg = 0.0;
    for (offset, probability) in sprob.iter().enumerate() {
        score_avg += (min_score + offset as i32) as f64 * probability;
    }
    ScoreFreqProfile::from_ncbi_parts(min_score, sprob, min_score, max_score, score_avg)
}

// ---------------------------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2826-2898
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2360-2378
// Each diagnostic is computed only when `request` asks for it (`PSIDiagnosticsResponseNew`
// allocates only those, blast_psi.c:618-714).
pub(super) fn save_diagnostics(
    msa: &Msa,
    aligned_block: &AlignedBlock,
    seq_weights: &SequenceWeights,
    internal: &InternalPssm,
    request: PsiDiagnosticsRequest,
) -> PsiDiagnostics {
    let information_content = request.information_content.then(|| {
        internal
            .freq_ratios
            .iter()
            .map(|ratios| {
                let mut info_sum = 0.0;
                for r in 0..BLASTAA_SIZE {
                    if seq_weights.std_prob[r] > K_EPSILON {
                        let q_over_p_estimate = ratios[r] / seq_weights.std_prob[r];
                        if q_over_p_estimate > K_EPSILON {
                            info_sum += ratios[r] * q_over_p_estimate.glibc_ln() / NCBIMATH_LN2;
                        }
                    }
                }
                info_sum
            })
            .collect()
    });
    let gapless_column_weights = request.gapless_column_weights.then(|| {
        (0..msa.query_length)
            .map(|p| {
                if msa.num_matching_seqs[p] > 1 && msa.cell[0][p].letter != X {
                    let weight = seq_weights.gapless_column_weights[p] / internal.pseudocounts[p];
                    weight * (seq_weights.sigma[p] / aligned_block.size[p] as f64 - 1.0)
                } else {
                    0.0
                }
            })
            .collect()
    });
    PsiDiagnostics {
        information_content,
        weighted_residue_freqs: request
            .weighted_residue_frequencies
            .then(|| seq_weights.match_weights.clone()),
        frequency_ratios: request
            .frequency_ratios
            .then(|| internal.freq_ratios.clone()),
        gapless_column_weights,
        sigma: request.sigma.then(|| seq_weights.sigma.clone()),
        interval_sizes: request.interval_sizes.then(|| aligned_block.size.clone()),
        num_matching_seqs: request
            .num_matching_seqs
            .then(|| msa.num_matching_seqs.clone()),
    }
}

impl InternalPssm {
    /// `internal_pssm->pssm`, which `s_PSISavePssm` copies (`_PSICopyMatrix_int`,
    /// blast_psi.c:451-452).
    pub(super) fn into_scores(self) -> Vec<[i32; BLASTAA_SIZE]> {
        self.pssm
    }
}

#[cfg(test)]
#[path = "stages_tests.rs"]
mod tests;
