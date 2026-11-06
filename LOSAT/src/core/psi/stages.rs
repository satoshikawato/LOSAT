//! The stages of PSSM creation: the port of NCBI's `algo/blast/core/blast_psi_priv.c`
//! (598d8ae6), in the order `PSICreatePssmWithDiagnostics` runs them (`mod.rs`). Integer fields
//! keep NCBI's types (`Uint4` sizes wrap as NCBI's do); the floating-point operations are NCBI's,
//! in NCBI's order, with the C library's `log`, `exp` and `pow`.

use super::{
    KarlinBlk, PsiDiagnostics, PsiDiagnosticsRequest, PsiError, PsiMsa, PsiMsaCell, PsiScoreBlk,
};
use crate::stats::karlin_calc::{compute_karlin_params_ungapped, ScoreFreqProfile};
use crate::utils::matrix::{ncbistdaa, BLASTAA_SIZE, BLAST_SCORE_MIN};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:55-61
// ```c
// const double kPSINearIdentical = 0.94;
// const double kPSIIdentical = 1.0;
// const unsigned int kQueryIndex = 0;
// const double kEpsilon = 0.0001;
// const int kPSIScaleFactor = 200;
// const double kPositScalingPercent = 0.05;
// const Uint4 kPositScalingNumIterations = 10;
// ```
const K_PSI_NEAR_IDENTICAL: f64 = 0.94;
const K_PSI_IDENTICAL: f64 = 1.0;
const K_QUERY_INDEX: usize = 0;
const K_EPSILON: f64 = 0.0001;
const K_PSI_SCALE_FACTOR: i32 = 200;
const K_POSIT_SCALING_PERCENT: f64 = 0.05;
const K_POSIT_SCALING_NUM_ITERATIONS: u32 = 10;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:550
// ```c
// #define EFFECTIVE_ALPHABET 20 /**< size of alphabet used for pseudocounts calculations */
// ```
const EFFECTIVE_ALPHABET: usize = 20;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2065-2066
// ```c
// #define MAX_IND_OBSERVATIONS  400 /**< max number of independent observation for pseudocount calculation */
// #define PSEUDO_MAX 1000000 /**< effective infinity */
// ```
const MAX_IND_OBSERVATIONS: usize = 400;
const PSEUDO_MAX: f64 = 1000000.0;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/ncbi_math.h:161
// ```c
// #define NCBIMATH_LN2	0.69314718055994530941723212145818
// ```
const NCBIMATH_LN2: f64 = 0.693_147_180_559_945_309_417_232_121_458_18;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_stat.h:122
// ```c
// #define BLAST_SCORE_MAX INT2_MAX   /**< maximum allowed score (for one letter comparison). */
// ```
const BLAST_SCORE_MAX: i32 = i16::MAX as i32;

const GAP: u8 = ncbistdaa::GAP;
const X: u8 = ncbistdaa::X;
const STAR: u8 = ncbistdaa::STOP;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
// ```c
// long BLAST_Nint(double x)
// {
//    x += (x >= 0. ? 0.5 : -0.5);
//    return (long)x;
// }
// ```
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
    // ```c
    //     /* Copy the multiple sequence alignment data */
    //     for (s = 0; s < msa->dimensions->num_seqs + 1; s++) {
    //         for (p = 0; p < msa->dimensions->query_length; p++) {
    //             ASSERT(msa->data[s][p].letter <= BLASTAA_SIZE);
    //             retval->data[s][p].letter = msa->data[s][p].letter;
    //             retval->data[s][p].is_aligned = msa->data[s][p].is_aligned;
    //         }
    //     }
    //     ...
    //     for (s = 0; s < msa->dimensions->num_seqs + 1; s++) {
    //         retval->use_sequence[s] = TRUE;
    //     }
    // ```
    pub(super) fn new(msa: &PsiMsa) -> Self {
        Self {
            query_length: msa.query_length,
            num_seqs: msa.rows.len() - 1,
            data: msa.rows.clone(),
            use_sequence: vec![true; msa.rows.len()],
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:209-224
    // ```c
    // _PSIPackedMsaGetNumberOfAlignedSeqs(const _PSIPackedMsa* msa)
    // {
    //     ...
    //     for (i = 0; i < msa->dimensions->num_seqs + 1; i++) {
    //         if (msa->use_sequence[i]) {
    //             retval++;
    //         }
    //     }
    // ```
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
    // ```c
    //     retval->dimensions->num_seqs = _PSIPackedMsaGetNumberOfAlignedSeqs(msa);
    //
    //     retval->cell = (_PSIMsaCell**)
    //         _PSIAllocateMatrix(retval->dimensions->num_seqs + 1,
    //                            retval->dimensions->query_length,
    //                            sizeof(_PSIMsaCell));
    //     ...
    //         for (s = 0; s < msa->dimensions->num_seqs + 1; s++) {
    //
    //             if ( !msa->use_sequence[s] ) {
    //                 continue;
    //             }
    //
    //             for (p = 0; p < retval->dimensions->query_length; p++) {
    //                 retval->cell[ss][p].letter = msa->data[s][p].letter;
    //                 retval->cell[ss][p].is_aligned = msa->data[s][p].is_aligned;
    //                 retval->cell[ss][p].extents.left = -1;
    //                 retval->cell[ss][p].extents.right =
    //                     msa->dimensions->query_length;
    //             }
    //             ss++;
    //         }
    //     ...
    //     for (p = 0; p < retval->dimensions->query_length; p++) {
    //         ASSERT(IS_residue(msa->data[kQueryIndex][p].letter));
    //         retval->query[p] = msa->data[kQueryIndex][p].letter;
    //     }
    //     ...
    //     _PSIUpdatePositionCounts(retval);
    // ```
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
    // ```c
    //     for (s = 0; s < kNumberOfSeqs; s++) {
    //         _PSIMsaCell* pos = msa->cell[s];   /* entry in MSA matrix */
    //
    //         for (p = 0; p < kQueryLength; p++, pos++) {
    //             if (pos->is_aligned) {
    //                 const Uint1 res = pos->letter;
    //                 if (res >= msa->alphabet_size) {
    //                     continue;
    //                 }
    //                 msa->residue_counts[p][res]++;
    //                 msa->num_matching_seqs[p]++;
    //             }
    //         }
    //     }
    // ```
    fn update_position_counts(&mut self) {
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
    // ```c
    //     retval->size = (Uint4*) calloc(query_length, sizeof(Uint4));
    //     ...
    //     for (i = 0; i < query_length; i++) {
    //         retval->pos_extnt[i].left = -1;
    //         retval->pos_extnt[i].right = query_length;
    //     }
    // ```
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
    // ```c
    //     retval->row_sigma = (double*) calloc(dimensions->num_seqs + 1,
    //                                          sizeof(double));
    //     ...
    //     retval->norm_seq_weights = (double*) calloc(dimensions->num_seqs + 1,
    //                                                 sizeof(double));
    //     ...
    //     retval->sigma = (double*) calloc(dimensions->query_length, sizeof(double));
    //     ...
    //     retval->match_weights = (double**)
    //         _PSIAllocateMatrix(dimensions->query_length,
    //                            sbp->alphabet_size,
    //                            sizeof(double));
    //     ...
    //     retval->std_prob = BLAST_GetStandardAaProbabilities();
    //     ...
    //     retval->gapless_column_weights = (double*)
    //         calloc(dimensions->query_length, sizeof(double));
    //     ...
    //     retval->posDistinctDistrib = (int**)
    //         _PSIAllocateMatrix(1+dimensions->query_length,
    //                            1+EFFECTIVE_ALPHABET,
    //                            sizeof(int));
    //     ...
    //     retval->posNumParticipating = (int*) calloc(1+dimensions->query_length, sizeof(int));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1341-1346
    // ```c
    //     standard_probabilities = Blast_ResFreqNew(&sbp);
    //     Blast_ResFreqStdComp(&sbp, standard_probabilities);
    //
    //     for (i = 0; i < (Uint4) sbp.alphabet_size; i++) {
    //         retval[i] = standard_probabilities->prob[i];
    //     }
    // ```
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:435-460
    // ```c
    //     retval->ncols = query_length;
    //     retval->nrows = alphabet_size;
    //
    //     retval->pssm = (int**) _PSIAllocateMatrix(retval->ncols,
    //                                               retval->nrows,
    //                                               sizeof(int));
    //     ...
    //     retval->pseudocounts = (double*) calloc(query_length, sizeof(double));
    // ```
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
// ```c
// _PSIPurgeBiasedSegments(_PSIPackedMsa* msa)
// {
//     ...
//     s_PSIPurgeSelfHits(msa);
//     s_PSIPurgeNearIdenticalAlignments(msa);
//     ...
// s_PSIPurgeSelfHits(_PSIPackedMsa* msa)
// {
//     ...
//     for (s = kQueryIndex + 1; s < msa->dimensions->num_seqs + 1; s++) {
//         s_PSIPurgeSimilarAlignments(msa, kQueryIndex, s, kPSIIdentical);
//     }
// }
// ...
//     for (i = 1; i < msa->dimensions->num_seqs + 1; i++) {
//         for (j = 1; (i + j) < msa->dimensions->num_seqs + 1; j++) {
//             /* N.B.: The order of comparison of sequence pairs is deliberate,
//              * tests on real data indicated that this approach allowed more
//              * sequences to be purged */
//             s_PSIPurgeSimilarAlignments(msa, j, (i + j), kPSINearIdentical);
//         }
//     }
// ```
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
// ```c
// _PSIResetAlignmentTraits(_PSIAlignmentTraits* traits, Uint4 position)
// {
//     ASSERT(traits);
//     memset((void*) traits, 0, sizeof(_PSIAlignmentTraits));
//     traits->start = position;
// }
// ```
fn reset_alignment_traits(traits: &mut AlignmentTraits, position: u32) {
    *traits = AlignmentTraits {
        start: position,
        ..AlignmentTraits::default()
    };
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1076-1112
// ```c
//     switch (*state) {
//     case eCounting:
//         /* Purge aligned region if max_percent_identity is exceeded */
//         {
//             if (traits->effective_length > 0) {
//                 const double percent_identity =
//                     ((double)traits->n_identical) / traits->effective_length;
//                 if (percent_identity >= max_percent_identity) {
//                     const unsigned int align_stop =
//                         traits->start + traits->effective_length +
//                         traits->n_x_residues;
//                     int rv = _PSIPurgeAlignedRegion(msa, seq_index,
//                                                     traits->start, align_stop);
//                     ...
//         *state = eResting;
//         break;
//
//     case eResting:
//         /* No-op */
//         break;
// ```
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
// ```c
//     if ( seq_index1 == seq_index2 ||
//          !msa->use_sequence[seq_index1] ||
//          !msa->use_sequence[seq_index2] ) {
//         return;
//     }
//
//     _PSIResetAlignmentTraits(&traits, p);
//     ...
//     for (p = 0; p < kQueryLength; p++, seq1++, seq2++) {
//         Boolean kPos1Aligned = (seq_index1 == kQueryIndex ?
//                                       FALSE : seq1->is_aligned);
//         Boolean kPos2Aligned = seq2->is_aligned;
//
//         /* if neither position is aligned */
//         if (!kPos1Aligned && !kPos2Aligned) {
//             _handleNeitherAligned(&traits, &state, msa, seq_index2,
//                                   max_percent_identity);
//         } else {
//             Boolean neither_is_X =
//                seq1->letter != kXResidue && seq2->letter != kXResidue;
//
//             /* Either one of the is aligned*/
//             if (neither_is_X) {
//                 _handleEitherAlignedNeitherX(&traits, &state, p);
//             } else {    /* either is X */
//                 _handleEitherAlignedEitherX(&traits, &state);
//             }
//
//             if (neither_is_X && (kPos2Aligned && seq1->is_aligned) &&
//                 (seq1->letter == seq2->letter)) {
//                 _handleBothAlignedSameResidueNoX(&traits, &state);
//             }
//         }
//     }
//     _handleNeitherAligned(&traits, &state, msa, seq_index2,
//                           max_percent_identity);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1118-1183
// ```c
// _handleBothAlignedSameResidueNoX(...)
//     case eCounting:
//         traits->n_identical++;
// ...
// _handleEitherAlignedEitherX(...)
//     case eCounting:
//         traits->n_x_residues++;
// ...
// _handleEitherAlignedNeitherX(...)
//     case eCounting:
//         traits->effective_length++;
//         break;
//
//     case eResting:
//         _PSIResetAlignmentTraits(traits, position);
//         traits->effective_length = 1;   /* count this residue */
//         *state = eCounting;
//         break;
// ```
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
// ```c
// s_PSIDiscardIfUnused(_PSIPackedMsa* msa, unsigned int seq_index)
// {
//     ...
//     for (i = 0; i < msa->dimensions->query_length; i++) {
//         if (msa->data[seq_index][i].is_aligned) {
//             contains_aligned_regions = TRUE;
//             break;
//         }
//     }
//
//     if ( !contains_aligned_regions ) {
//         msa->use_sequence[seq_index] = FALSE;
//     }
// }
//
// int
// _PSIPurgeAlignedRegion(_PSIPackedMsa* msa,
//                        unsigned int seq_index,
//                        unsigned int start,
//                        unsigned int stop)
// {
//     ...
//     if ( !msa ||
//         seq_index == 0 ||
//         seq_index > msa->dimensions->num_seqs + 1 ||
//         stop > msa->dimensions->query_length) {
//         return PSIERR_BADPARAM;
//     }
//
//     sequence_position = msa->data[seq_index];
//     for (i = start; i < stop; i++) {
//         sequence_position[i].letter = 0;
//         sequence_position[i].is_aligned = FALSE;
//     }
//
//     s_PSIDiscardIfUnused(msa, seq_index);
// ```
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
// ```c
//     retval = s_PSIValidateNoFlankingGaps(msa);
//     ...
//     if ( !ignore_unaligned_positions ) {
//         retval = s_PSIValidateAlignedColumns(msa);
//     ...
//     retval = s_PSIValidateNoGapsInQuery(msa);
//     ...
//     retval = s_PSIValidateParticipatingSequences(msa);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:703-797
// ```c
//     /* Look for starting gaps in alignments */
//     for (s = 0; s < msa->dimensions->num_seqs + 1; s++) {
//         /* find the first aligned residue */
//         for (p = 0; p < (Int4) msa->dimensions->query_length; p++) {
//             if (msa->cell[s][p].is_aligned) {
//                 if (msa->cell[s][p].letter == kGapResidue) {
//                     return PSIERR_STARTINGGAP;
//                 } else {
//                     break;
//     ...
//     /* Look for ending gaps in alignments */
//     ...
//         for (p = msa->dimensions->query_length - 1; p >= 0; p--) {
//             if (msa->cell[s][p].is_aligned) {
//                 if (msa->cell[s][p].letter == kGapResidue) {
//                     return PSIERR_ENDINGGAP;
//     ...
//     for (p = 0; p < msa->dimensions->query_length; p++) {
//         Boolean found_aligned_sequence = FALSE;
//         Boolean found_non_gap_residue = FALSE;
//         for (s = kQueryIndex; s < msa->dimensions->num_seqs + 1; s++) {
//             if (msa->cell[s][p].is_aligned) {
//                 found_aligned_sequence = TRUE;
//                 if (msa->cell[s][p].letter != kGapResidue) {
//                     found_non_gap_residue = TRUE;
//                     break;
//         ...
//         if ( !found_aligned_sequence ) {
//             return PSIERR_UNALIGNEDCOLUMN;
//         }
//         if ( !found_non_gap_residue ) {
//             return PSIERR_COLUMNOFGAPS;
//         }
//     ...
//     return msa->dimensions->num_seqs ? PSI_SUCCESS : PSIERR_NOALIGNEDSEQS;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:683-696
// ```c
//     for (p = 0; p < msa->dimensions->query_length; p++) {
//         if (msa->cell[kQueryIndex][p].letter == kGapResidue ||
//             msa->query[p] == kGapResidue) {
//             return PSIERR_GAPINQUERY;
//         }
//     }
// ```
pub(super) fn validate_msa(msa: &Msa) -> Result<(), PsiError> {
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
    for p in 0..msa.query_length {
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

// ---------------------------------------------------------------------------------------------
// Alignment blocks
// ---------------------------------------------------------------------------------------------

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1306-1313
// ```c
//     /* no need to compute extents for query sequence */
//     for (s = kQueryIndex + 1; s < msa->dimensions->num_seqs + 1; s++) {
//         _PSIGetLeftExtents(msa, s);
//         _PSIGetRightExtents(msa, s);
//         _PSIComputePositionExtents(msa, s, aligned_blocks);
//     }
//
//     _PSIComputeAlignedRegionLengths(msa, aligned_blocks);
// ```
pub(super) fn compute_alignment_blocks(msa: &mut Msa, aligned_blocks: &mut AlignedBlock) {
    for s in K_QUERY_INDEX + 1..msa.num_seqs + 1 {
        get_left_extents(&mut msa.cell[s]);
        get_right_extents(&mut msa.cell[s]);
        compute_position_extents(&msa.cell[s], aligned_blocks);
    }
    compute_aligned_region_lengths(msa, aligned_blocks);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1330-1349
// ```c
//     if (sequence_position[prev].is_aligned &&
//         sequence_position[prev].letter != kGapResidue) {
//         sequence_position[prev].extents.left = prev;
//     }
//
//     for (curr = prev + 1; curr < msa->dimensions->query_length;
//          curr++, prev++) {
//
//         if ( !sequence_position[curr].is_aligned ) {
//             continue;
//         }
//
//         if (sequence_position[prev].is_aligned) {
//             sequence_position[curr].extents.left =
//                 sequence_position[prev].extents.left;
//         } else {
//             sequence_position[curr].extents.left = curr;
//         }
//     }
// ```
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
// ```c
//     last = msa->dimensions->query_length - 1;
//
//     if (sequence_position[last].is_aligned &&
//         sequence_position[last].letter != kGapResidue) {
//         sequence_position[last].extents.right = last;
//     }
//
//     for (curr = last - 1; curr >= 0; curr--, last--) {
//
//         if ( !sequence_position[curr].is_aligned ) {
//             continue;
//         }
//
//         if (sequence_position[last].is_aligned) {
//             sequence_position[curr].extents.right =
//                 sequence_position[last].extents.right;
//         } else {
//             sequence_position[curr].extents.right = curr;
//         }
//     }
// ```
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
// ```c
//     for (i = 0; i < msa->dimensions->query_length; i++) {
// #ifdef PSI_IGNORE_GAPS_IN_COLUMNS
//         ...
// #else
//         if (sequence_position[i].is_aligned) {
// #endif
//             aligned_blocks->pos_extnt[i].left =
//                 MAX(aligned_blocks->pos_extnt[i].left,
//                     sequence_position[i].extents.left);
//             aligned_blocks->pos_extnt[i].right =
//                 MIN(aligned_blocks->pos_extnt[i].right,
//                     sequence_position[i].extents.right);
//         }
//     }
// ```
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
// ```c
//     for (i = 0; i < kQueryLength; i++) {
//         aligned_blocks->size[i] = aligned_blocks->pos_extnt[i].right -
//                                    aligned_blocks->pos_extnt[i].left + 1;
//         ...
//     }
//
//     /* Do not include X's in aligned region lengths */
//     for (i = 0; i < kQueryLength; i++) {
//
//         if (msa->query[i] == kXResidue) {
//
//             Uint4 idx = 0;
//             for (idx = 0; idx < i; idx++) {
//                 if ((Uint4)aligned_blocks->pos_extnt[idx].right >= i &&
//                     msa->query[idx] != kXResidue) {
//                     aligned_blocks->size[idx]--;
//                 }
//             }
//             for (idx = msa->dimensions->query_length - 1; idx > i; idx--) {
//                 if ((Uint4)aligned_blocks->pos_extnt[idx].left <= i &&
//                     msa->query[idx] != kXResidue) {
//                     aligned_blocks->size[idx]--;
//                 }
//             }
//         }
//
//     }
// ```
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
// ```c
//     for (pos = 0; pos < kQueryLength; pos++) {
//
//         /* ignore positions of no interest */
//         if (aligned_blocks->size[pos] == 0 ||
//             msa->num_matching_seqs[pos] <= kExpectedNumMatchingSeqs) {
//             continue;
//         }
//
//         DynamicUint4Array_Copy(prev_pos_aligned_seqs, aligned_seqs);
//         _PSIGetAlignedSequencesForPosition(msa, pos, aligned_seqs);
//         ASSERT(msa->num_matching_seqs[pos] == aligned_seqs->num_used);
//         if (aligned_seqs->num_used <= kExpectedNumMatchingSeqs) {
//             continue;
//         }
//         last_calc_pos = pos;
//
//         if (last_calc_pos != pos - 1 ||
//             !DynamicUint4Array_AreEqual(aligned_seqs, prev_pos_aligned_seqs)) {
//
//             memset((void*)seq_weights->norm_seq_weights, 0,
//                    sizeof(double)*(msa->dimensions->num_seqs+1));
//             memset((void*)seq_weights->row_sigma, 0,
//                    sizeof(double)*(msa->dimensions->num_seqs+1));
//
//             _PSICalculateNormalizedSequenceWeights(msa, aligned_blocks, pos,
//                                                    aligned_seqs, seq_weights);
//         } else {
//             ...
//         }
//         seq_weights->posNumParticipating[pos] = aligned_seqs->num_used;
//
//         /* Uses seq_weights->norm_seq_weights to populate match_weights */
//         _PSICalculateMatchWeights(msa, pos, aligned_seqs, seq_weights);
//     }
//     ...
//     retval = _PSICheckSequenceWeights(msa, seq_weights,
//                                       nsg_compatibility_mode);
//     ...
// #ifndef PSI_IGNORE_GAPS_IN_COLUMNS
//     _PSISpreadGapWeights(msa, seq_weights, nsg_compatibility_mode);
//     retval = _PSICheckSequenceWeights(msa, seq_weights,
//                                       nsg_compatibility_mode);
// #endif
// ```
// `last_calc_pos` is set to `pos` just before the test, so `last_calc_pos != pos - 1` always
// holds and the weights are computed at every position: the `else` branch is never taken.
// Without the structure group compatibility mode, `kExpectedNumMatchingSeqs` is 1.
pub(super) fn compute_sequence_weights(
    msa: &Msa,
    aligned_blocks: &AlignedBlock,
    seq_weights: &mut SequenceWeights,
) -> Result<(), PsiError> {
    const K_EXPECTED_NUM_MATCHING_SEQS: u32 = 1;
    let mut aligned_seqs: Vec<usize> = Vec::with_capacity(msa.num_seqs + 1);
    for pos in 0..msa.query_length {
        if aligned_blocks.size[pos] == 0
            || msa.num_matching_seqs[pos] <= K_EXPECTED_NUM_MATCHING_SEQS
        {
            continue;
        }
        get_aligned_sequences_for_position(msa, pos, &mut aligned_seqs);
        if aligned_seqs.len() as u32 <= K_EXPECTED_NUM_MATCHING_SEQS {
            continue;
        }
        seq_weights.norm_seq_weights.fill(0.0);
        seq_weights.row_sigma.fill(0.0);
        calculate_normalized_sequence_weights(msa, aligned_blocks, pos, &aligned_seqs, seq_weights);
        seq_weights.pos_num_participating[pos] = aligned_seqs.len() as i32;
        calculate_match_weights(msa, pos, &aligned_seqs, seq_weights);
    }
    check_sequence_weights(msa, seq_weights)?;
    spread_gap_weights(msa, seq_weights);
    check_sequence_weights(msa, seq_weights)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1915-1928
// ```c
//     aligned_sequences->num_used = 0;    /* reset the array */
//
//     for (i = 0; i < msa->dimensions->num_seqs + 1; i++) {
// #ifdef PSI_IGNORE_GAPS_IN_COLUMNS
//         ...
// #else
//         if (msa->cell[i][position].is_aligned) {
// #endif
//             DynamicUint4Array_Append(aligned_sequences, i);
//         }
//     }
// ```
fn get_aligned_sequences_for_position(msa: &Msa, position: usize, aligned: &mut Vec<usize>) {
    aligned.clear();
    for (i, row) in msa.cell.iter().enumerate() {
        if row[position].is_aligned {
            aligned.push(i);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:1786-1871
// ```c
//     for (i = (Uint4)aligned_blocks->pos_extnt[position].left;
//          i <= (Uint4)aligned_blocks->pos_extnt[position].right; i++) {
//
//         Uint4 residue_counts_for_column[BLASTAA_SIZE] = { 0 };
//         Uint4 num_distinct_residues_for_column = 0;
//         Uint4 num_local_std_letters = 0;
//         ...
//         for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//             const Uint4 kSeqIdx = aligned_seqs->data[asi];
//             const Uint1 kResidue = msa->cell[kSeqIdx][i].letter;
//
//             if (residue_counts_for_column[kResidue] == 0) {
//                 num_distinct_residues_for_column++;
//                 if (kResidue != kGapResidue && kResidue != kXResidue)
//                     num_local_std_letters++;
//             }
//             residue_counts_for_column[kResidue]++;
//         }
//
//         sigma += num_distinct_residues_for_column;
//         num_local_std_letters = MIN(num_local_std_letters,EFFECTIVE_ALPHABET);
//         seq_weights->posDistinctDistrib[position][num_local_std_letters]++;
//         if (num_distinct_residues_for_column > 1) {
//             distinct_residues_found = TRUE;
//         }
//
//         for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//             const Uint4 seq_idx = aligned_seqs->data[asi];
//             const Uint1 residue = msa->cell[seq_idx][i].letter;
//             seq_weights->row_sigma[seq_idx] +=
//                 (1.0 / (double)
//                  (residue_counts_for_column[residue] *
//                   num_distinct_residues_for_column) );
//         }
//     }
//
//     /* Save sigma for this position */
//     seq_weights->sigma[position] = sigma;
//
//     if (distinct_residues_found) {
//         double weight_sum = 0.0;
//
//         for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//             const Uint4 seq_idx = aligned_seqs->data[asi];
//             seq_weights->norm_seq_weights[seq_idx] =
//                 seq_weights->row_sigma[seq_idx] /
//                 (aligned_blocks->pos_extnt[position].right -
//                  aligned_blocks->pos_extnt[position].left + 1);
//             weight_sum += seq_weights->norm_seq_weights[seq_idx];
//         }
//
//         /* Normalize */
//         for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//             const Uint4 seq_idx = aligned_seqs->data[asi];
//             seq_weights->norm_seq_weights[seq_idx] /= weight_sum;
//         }
//
//     } else {
//         for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//             const Uint4 seq_idx = aligned_seqs->data[asi];
//             seq_weights->norm_seq_weights[seq_idx] =
//                 (1.0/(double) aligned_seqs->num_used);
//         }
//     }
// ```
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
// ```c
//     for (asi = 0; asi < aligned_seqs->num_used; asi++) {
//         const Uint4 seq_idx = aligned_seqs->data[asi];
//         const Uint1 residue = msa->cell[seq_idx][position].letter;
//
//         seq_weights->match_weights[position][residue] +=
//             seq_weights->norm_seq_weights[seq_idx];
//
//         /* Collected for diagnostics information, not used elsewhere */
//         if (residue != kGapResidue) {
//             seq_weights->gapless_column_weights[position] +=
//              seq_weights->norm_seq_weights[seq_idx];
//         }
//     }
// ```
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
// ```c
//     for (pos = 0; pos < msa->dimensions->query_length; pos++) {
//
//         if (msa->num_matching_seqs[pos] <= kExpectedNumMatchingSeqs ||
//             msa->cell[kQueryIndex][pos].letter == kXResidue) {
//             continue;
//         }
//
//         /* Disperse method of spreading the gap weights */
//         for (res = 0; res < msa->alphabet_size; res++) {
//             if (seq_weights->std_prob[res] > kEpsilon) {
//                 seq_weights->match_weights[pos][res] +=
//                     (seq_weights->match_weights[pos][kGapResidue] *
//                      seq_weights->std_prob[res]);
//             }
//         }
//         seq_weights->match_weights[pos][kGapResidue] = 0.0;
//
//     }
// ```
fn spread_gap_weights(msa: &Msa, seq_weights: &mut SequenceWeights) {
    for pos in 0..msa.query_length {
        if msa.num_matching_seqs[pos] <= 1 || msa.cell[K_QUERY_INDEX][pos].letter == X {
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
// ```c
//     for (pos = 0; pos < msa->dimensions->query_length; pos++) {
//
//         double running_total = 0.0;
//         Uint4 residue = 0;
//
//         if (msa->num_matching_seqs[pos] <= kExpectedNumMatchingSeqs ||
//             msa->cell[kQueryIndex][pos].letter == kXResidue) {
//             ...
//             continue;
//         }
//
//         for (residue = 0; residue < msa->alphabet_size; residue++) {
//             running_total += seq_weights->match_weights[pos][residue];
//         }
//
//         if (running_total < 0.99 || running_total > 1.01) {
//             return PSIERR_BADSEQWEIGHTS;
//         }
// ```
// `!(0.99..=1.01).contains(..)` is NCBI's two comparisons for every number but NaN, which
// NCBI's lets through and this rejects. The weights are finite (sums of finite products,
// normalised by a positive total in `_PSICalculateNormalizedSequenceWeights`), so no NaN
// reaches this check.
fn check_sequence_weights(msa: &Msa, seq_weights: &SequenceWeights) -> Result<(), PsiError> {
    for pos in 0..msa.query_length {
        if msa.num_matching_seqs[pos] <= 1 || msa.cell[K_QUERY_INDEX][pos].letter == X {
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
// ```c
//     s_initializeExpNumObservations(&(expno[0]),  backgroundProbabilities);
//
//     for (p = 0; p < msa->dimensions->query_length; p++) {
//         double columnCounts = 0.0; /*column-specific pseudocounts*/
//         double observations = 0.0;
//         double pseudoWeight; /*multiplier for pseudocounts term*/
//         if (msa->cell[kQueryIndex][p].letter != kXResidue)
//         {
//            observations = s_effectiveObservations(aligned_blocks, seq_weights, p, msa->dimensions->query_length, expno);
//            ...
//            if (0 == pseudo_count)
//               columnCounts = s_columnSpecificPseudocounts(seq_weights, p, backgroundProbabilities, observations);
//            else
//               columnCounts = pseudo_count;
//         }
//         if (columnCounts >= PSEUDO_MAX) {
//            pseudoWeight = kZeroObsPseudo;
//            observations = 0;
//         }
//         else {
//            pseudoWeight = columnCounts;
//         }
//
//         for (r = 0; r < msa->alphabet_size; r++) {
//
//             if (msa->cell[kQueryIndex][p].letter == kXResidue ||
//                 seq_weights->std_prob[r] <= kEpsilon) {
//                 internal_pssm->freq_ratios[p][r] = 0.0;
//             } else {
//                 ...
//                 const double kBeta = pseudoWeight;
//                 ...
//                 internal_pssm->pseudocounts[p] = kBeta;
//
//                 for (i = 0; i < msa->alphabet_size; i++) {
//                     if (sbp->matrix->data[r][i] != BLAST_SCORE_MIN) {
//                         pseudo += (seq_weights->match_weights[p][i] *
//                                    freq_ratios->data[r][i]);
//                     }
//                 }
//                 pseudo *= kBeta;
//
//                 numerator =
//                     (observations * seq_weights->match_weights[p][r] /
//                      seq_weights->std_prob[r])
//                     + pseudo;
//
//                 denominator = observations + kBeta;
//                 ...
//                 qOverPEstimate = numerator/denominator;
//
//                 internal_pssm->freq_ratios[p][r] = qOverPEstimate *
//                     seq_weights->std_prob[r];
//
//             }
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2084-2085
// ```c
//     const double kZeroObsPseudo = 30.0; /*arbitrary constant to use for columns with
//                              zero observations in actual data (ZERO_OBS_PSEUDO in posit.c) */
// ```
pub(super) fn compute_freq_ratios(
    msa: &Msa,
    seq_weights: &mut SequenceWeights,
    sbp: &PsiScoreBlk,
    aligned_blocks: &AlignedBlock,
    pseudo_count: i32,
    internal: &mut InternalPssm,
) {
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
            let q_over_p_estimate = numerator / denominator;
            internal.freq_ratios[p][r] = q_over_p_estimate * seq_weights.std_prob[r];
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3068-3084
// ```c
//    expno[0] = 0;
//    for (j=1;j<MAX_IND_OBSERVATIONS;++j) {
//      weighted_sum = 0;
//      for (k=0;k<EFFECTIVE_ALPHABET;++k)
//        weighted_sum += exp(j*log(1.0-backgroundProbabilities[k]));
//      expno[j] = EFFECTIVE_ALPHABET-weighted_sum;
//    }
// ```
// `expno[MAX_IND_OBSERVATIONS]` is never written nor read (`s_effectiveObservations` stops at
// `MAX_IND_OBSERVATIONS - 1`).
fn initialize_exp_num_observations(background: &[f64; 20]) -> [f64; MAX_IND_OBSERVATIONS] {
    let mut expno = [0.0; MAX_IND_OBSERVATIONS];
    for (j, value) in expno.iter_mut().enumerate().skip(1) {
        let mut weighted_sum = 0.0;
        for &probability in background.iter().take(EFFECTIVE_ALPHABET) {
            weighted_sum += (j as f64 * (1.0 - probability).ln()).exp();
        }
        *value = EFFECTIVE_ALPHABET as f64 - weighted_sum;
    }
    expno
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3135-3175
// ```c
//  if (align_blk->pos_extnt[columnNumber].left < 0)
//    return(0);
//  if (align_blk->pos_extnt[columnNumber].right >= queryLength)
//    return(0);
//
//  halfNumColumns = MAX(1,(align_blk->pos_extnt[columnNumber].right -
//                          align_blk->pos_extnt[columnNumber].left+2)/2);
//  k = EFFECTIVE_ALPHABET;
//  columnsAccountedFor = 0;
//  totalDistinctCounts = 0;
//  while (columnsAccountedFor < halfNumColumns) {
//    ASSERT(k >= 0);
//    totalDistinctCounts += (seq_weights->posDistinctDistrib[columnNumber][k] *k);
//    columnsAccountedFor += seq_weights->posDistinctDistrib[columnNumber][k];
//    if (columnsAccountedFor > halfNumColumns) {
//      totalDistinctCounts -=
//        ((columnsAccountedFor - halfNumColumns) * k);
//      columnsAccountedFor = halfNumColumns;
//    }
//    k--;
//  }
//  aveDistinctAA = ((double) totalDistinctCounts)/
//    ((double) columnsAccountedFor);
//
//  for (i=1;i<MAX_IND_OBSERVATIONS && expno[i]<=aveDistinctAA;++i);
//  indep = (i==MAX_IND_OBSERVATIONS) ? i :
//    i-(expno[i]-aveDistinctAA)/(expno[i]-expno[i-1]);
//  indep = MIN(indep, seq_weights->posNumParticipating[columnNumber]);
//  indep = MAX(0,indep - 1);
//  return(indep);
// ```
// The distribution of a column whose extent was set holds one count per column of the extent
// (`_PSICalculateNormalizedSequenceWeights`), at least `halfNumColumns`, so `k` stays at 0 or
// above. The loop's `k >= 0` stands for NCBI's `ASSERT(k >= 0)`, which only debug builds check:
// a release build of NCBI would read before the row there; LOSAT stops instead (not reached).
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
// ```c
//    charOrder[0] =  1;  /*A*/
//    charOrder[1] =  16; /*R*/
//    charOrder[2] =  13; /*N*/
//    charOrder[3] =  4;  /*D*/
//    charOrder[4] =  3;  /*C*/
//    charOrder[5] =  15; /*Q*/
//    charOrder[6] =  5;  /*E*/
//    charOrder[7] =  7;  /*G*/
//    charOrder[8] =  8;  /*H*/
//    charOrder[9] =  9;  /*I*/
//    charOrder[10] = 11; /*L*/
//    charOrder[11] = 10; /*K*/
//    charOrder[12] = 12; /*M*/
//    charOrder[13] =  6; /*F*/
//    charOrder[14] = 14; /*P*/
//    charOrder[15] = 17; /*S*/
//    charOrder[16] = 18; /*T*/
//    charOrder[17] = 20; /*W*/
//    charOrder[18] = 22; /*Y*/
//    charOrder[19] = 19; /*V*/
//
//    for(c = 0; c < EFFECTIVE_ALPHABET; c++)
//      probabilities[c] = posSearch->match_weights[columnNumber][charOrder[c]];
// ```
const CHAR_ORDER: [usize; EFFECTIVE_ALPHABET] = [
    1, 16, 13, 4, 3, 15, 5, 7, 8, 9, 11, 10, 12, 6, 14, 17, 18, 20, 22, 19,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3086-3118
// ```c
//   const double kPseudoMult = 500.0;  /* Was PSEUDO_MULTIPLIER */
//   const double kPseudoNumerator = 0.0457;  /*numerator of entropy-based method, was PSEUDO_NUMERATOR */
//   const double kPseudoExponent = 0.8;  /*exponent of denominator, was PSEUDO_EXPONENT */
//   const double kPseudoSmallInitial = 5.5; /*small number of pseudocounts to
//                               avoid 0 probabilities in entropy-based method, was PSEUDO_SMALL_INITIAL */
//
//   s_fillColumnProbabilities(&(columnProbabilitiesInitial[0]), posSearch, columnNumber);
//   s_adjustColumnProbabilities(&(columnProbabilitiesInitial[0]),
//                             &(columnProbabilitiesAdjusted[0]),
//                               kPseudoSmallInitial,
//                               backgroundProbabilities, observations);
//   relativeEntropy = s_computeRelativeEntropy(&(columnProbabilitiesAdjusted[0]),
//                                            backgroundProbabilities);
//   pseudoDenominator = pow(relativeEntropy, kPseudoExponent);
//   alpha = kPseudoNumerator/pseudoDenominator;
//   if (alpha < (1.0 - kPosEpsilon))
//     returnValue = kPseudoMult * alpha/ (1- alpha);
//   else
//     returnValue = PSEUDO_MAX;
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:3014-3060
// ```c
//   overallSum = 0.0;
//   for(c = 0; c < EFFECTIVE_ALPHABET; c++) {
//     intermediateSums[c] =
//       (initialProbabilities[c] * observations) +
//       (standardProbabilities[c] * standardWeight);
//     overallSum += intermediateSums[c];
//   }
//   for(c = 0; c < EFFECTIVE_ALPHABET; c++)
//     probabilitiesToReturn[c] = intermediateSums[c]/overallSum;
// ...
// const double kPosEpsilon = 0.0001; /**< minimum return value of s_computeRelativeEntropy */
// ...
//    returnValue = 0.0;
//    for(c = 0; c < EFFECTIVE_ALPHABET; c++) {
//      if (newDistribution[c] > kPosEpsilon)
//        returnValue += (newDistribution[c] *
//                        log (newDistribution[c]/backgroundProbabilities[c]));
//    }
//    if (returnValue < kPosEpsilon)
//      returnValue = kPosEpsilon;
// ```
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
            relative_entropy += adjusted * (adjusted / background[c]).ln();
        }
    }
    if relative_entropy < K_POS_EPSILON {
        relative_entropy = K_POS_EPSILON;
    }
    let pseudo_denominator = relative_entropy.powf(K_PSEUDO_EXPONENT);
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
// ```c
//     ideal_lambda = sbp->kbp_ideal->Lambda;
//     freq_ratios = _PSIMatrixFrequencyRatiosNew(sbp->name);
//
//     /* Each column is a position in the query */
//     for (i = 0; i < internal_pssm->ncols; i++) {
//
//         /* True if all frequencies in column i are zero */
//         Boolean is_unaligned_column = TRUE;
//         const Uint4 kResidue = query[i];
//
//         for (j = 0; j < (Uint4) sbp->alphabet_size; j++) {
//
//             double qOverPEstimate = 0.0;
//
//             /* Division compensates for multiplication in
//              * _PSIComputeFreqRatios */
//             if (std_probs[j] > kEpsilon) {
//                 qOverPEstimate =
//                     internal_pssm->freq_ratios[i][j] / std_probs[j];
//             }
//
//             if (is_unaligned_column && qOverPEstimate != 0.0) {
//                 is_unaligned_column = FALSE;
//             }
//
//             /* Populate scaled matrix */
//             if (qOverPEstimate == 0.0 || std_probs[j] < kEpsilon) {
//                 internal_pssm->scaled_pssm[i][j] = BLAST_SCORE_MIN;
//             } else {
//                 double tmp = log(qOverPEstimate)/ideal_lambda;
//                 internal_pssm->scaled_pssm[i][j] = (int)
//                     BLAST_Nint(kPSIScaleFactor * tmp);
//             }
//
//             if ( (j == kXResidue || j == kStarResidue) &&
//                  (sbp->matrix->data[kResidue][kXResidue] != BLAST_SCORE_MIN) ) {
//                 internal_pssm->scaled_pssm[i][j] =
//                     sbp->matrix->data[kResidue][j] * kPSIScaleFactor;
//             }
//         }
//
//         if (is_unaligned_column) {
//             for (j = 0; j < (Uint4) sbp->alphabet_size; j++) {
//
//                 internal_pssm->pssm[i][j] = sbp->matrix->data[kResidue][j];
//
//                 if (freq_ratios->data[kResidue][j] != 0.0) {
//                     double tmp =
//                         kPSIScaleFactor * freq_ratios->bit_scale_factor *
//                         log(freq_ratios->data[kResidue][j])/NCBIMATH_LN2;
//
//                     internal_pssm->scaled_pssm[i][j] = (int)BLAST_Nint(tmp);
//                 } else {
//                     internal_pssm->scaled_pssm[i][j] = BLAST_SCORE_MIN;
//                 }
//             }
//         }
//     }
// ```
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
                let tmp = q_over_p_estimate.ln() / ideal_lambda;
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
                    let tmp = (K_PSI_SCALE_FACTOR * sbp.bit_scale_factor) as f64 * ratio.ln()
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
// ```c
//     factor = 1.0;
//     for ( ; ; ) {
//         ...
//         for (i = 0; i < internal_pssm->ncols; i++) {
//             for (j = 0; j < internal_pssm->nrows; j++) {
//                 if (scaled_pssm[i][j] != BLAST_SCORE_MIN) {
//                     pssm[i][j] =
//                         (int)BLAST_Nint(factor*scaled_pssm[i][j]/kPSIScaleFactor);
//                 } else {
//                     pssm[i][j] = BLAST_SCORE_MIN;
//                 }
//             }
//         }
//         _PSIUpdateLambdaK((const int**)pssm, query, query_length,
//                           std_probs, sbp);
//
//         new_lambda = sbp->kbp_psi[0]->Lambda;
//
//         if (new_lambda > ideal_lambda) {
//             if (first_time) {
//                 factor_high = 1.0 + kPositScalingPercent;
//                 factor = factor_high;
//                 factor_low = 1.0;
//                 too_high = TRUE;
//                 first_time = FALSE;
//             } else {
//                 if (too_high == FALSE) {
//                     break;
//                 }
//                 factor_high += (factor_high - 1.0);
//                 factor = factor_high;
//             }
//         } else if (new_lambda > 0) {
//             if (first_time) {
//                 factor_high = 1.0;
//                 factor_low = 1.0 - kPositScalingPercent;
//                 factor = factor_low;
//                 too_high = FALSE;
//                 first_time = FALSE;
//             } else {
//                 if (too_high == TRUE) {
//                     break;
//                 }
//                 factor_low += (factor_low - 1.0);
//                 factor = factor_low;
//             }
//         } else {
//             return PSIERR_POSITIVEAVGSCORE;
//         }
//     }
//
//     /* Binary search for kPositScalingNumIterations times */
//     for (index = 0; index < kPositScalingNumIterations; index++) {
//         ...
//         factor = (factor_high + factor_low)/2;
//         ... (the same rescaling and _PSIUpdateLambdaK)
//         if (new_lambda > ideal_lambda) {
//             factor_low = factor;
//         } else {
//             factor_high = factor;
//         }
//     }
// ```
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

/// The rescaling of `_PSIScaleMatrix`'s loops (quoted above `scale_matrix`).
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
// ```c
//     Blast_ScoreFreq* score_freqs =
//         _PSIComputeScoreProbabilities(pssm, query, query_length,
//                                       std_probs, sbp);
//
//     /* Calculate lambda and K */
//     Blast_KarlinBlkUngappedCalc(sbp->kbp_psi[0], score_freqs);
//
//     ASSERT(sbp->kbp_ideal);
//     ASSERT(sbp->kbp_psi[0]);
//     ASSERT(sbp->kbp_gap_std[0]);
//     ASSERT(sbp->kbp_gap_psi[0]);
//
//     sbp->kbp_gap_psi[0]->K =
//         sbp->kbp_psi[0]->K * sbp->kbp_gap_std[0]->K / sbp->kbp_ideal->K;
//     sbp->kbp_gap_psi[0]->logK = log(sbp->kbp_gap_psi[0]->K);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2699-2734
// ```c
//    kbp->Lambda = Blast_KarlinLambdaNR(sfp, BLAST_KARLIN_LAMBDA0_DEFAULT);
//    if (kbp->Lambda < 0.)
//       goto ErrExit;
//    ...
//    kbp->H = BlastKarlinLtoH(sfp, kbp->Lambda);
//    if (kbp->H < 0.)
//       goto ErrExit;
//    ...
//    kbp->K = BlastKarlinLHtoK(sfp, kbp->Lambda, kbp->H);
//    if (kbp->K < 0.)
//       goto ErrExit;
//    ...
// ErrExit:
//    kbp->Lambda = kbp->H = kbp->K = -1.;
// ```
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

/// `sbp->kbp_gap_psi[0]->K` after `_PSIUpdateLambdaK` (quoted above `update_lambda_k`).
pub(super) fn gapped_psi_k(kbp_psi: KarlinBlk, sbp: &PsiScoreBlk) -> f64 {
    kbp_psi.k * sbp.gap_std_k / sbp.ideal_k
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1795-1815
// ```c
// static BLAST_LetterProb Robinson_prob[] = {
//       { 'A', 78.05 },
//       { 'C', 19.25 },
//       { 'D', 53.64 },
//       ...
//       { 'W', 13.30 },
//       { 'Y', 32.16 }
//    };
// ```
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
// ```c
//     effective_length = _PSISequenceLengthWithoutX(query, query_length);
//
//     /* Get the minimum and maximum scores */
//     for (p = 0; p < query_length; p++) {
//         if (query[p] == kXResidue) {
//             continue;
//         }
//         for (r = 0; r < alphabet_size; r++) {
//             const int kScore = pssm[p][aa_alphabet[r]];
//
//             if (kScore <= BLAST_SCORE_MIN || kScore >= BLAST_SCORE_MAX) {
//                 continue;
//             }
//             max_score = MAX(kScore, max_score);
//             min_score = MIN(kScore, min_score);
//         }
//     }
//     ...
//     score_freqs = Blast_ScoreFreqNew(min_score, max_score);
//     ...
//     score_freqs->obs_min = min_score;
//     score_freqs->obs_max = max_score;
//     for (p = 0; p < query_length; p++) {
//         if (query[p] == kXResidue) {
//             continue;
//         }
//
//         for (r = 0; r < alphabet_size; r++) {
//             const int kScore = pssm[p][aa_alphabet[r]];
//
//             if (kScore <= BLAST_SCORE_MIN || kScore >= BLAST_SCORE_MAX) {
//                 continue;
//             }
//
//             /* Increment the weight for the score in position p, residue r */
//             score_freqs->sprob[kScore] +=
//                 (std_probs[aa_alphabet[r]]/effective_length);
//         }
//     }
//
//     ASSERT(score_freqs->score_avg == 0.0);
//     for (s = min_score; s <= max_score; s++) {
//         score_freqs->score_avg += (s * score_freqs->sprob[s]);
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2102-2104, 2118-2119 and
// 2703-2704
// ```c
//    if (lo >= 0 || hi <= 0 ||
//          lo < BLAST_SCORE_MIN || hi > BLAST_SCORE_MAX)
//       return 1;
//    ...
//    if (BlastScoreChk(score_min, score_max) != 0)
//       return NULL;
//    ...
//    if (kbp == NULL || sfp == NULL)
//       return 1;
// ```
// NCBI's `Blast_ScoreFreqNew` returns NULL when the PSSM's smallest score is 0 or more or its
// largest 0 or less; `Blast_KarlinBlkUngappedCalc` then leaves `kbp_psi[0]` as it was (the
// previous Lambda). This port does not model that branch: the scaled scores of a BLOSUM62 round
// (the only matrix LOSAT's psiblast accepts) keep negative and positive values at every factor
// the scaling tries, so it is not reached.
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
// ```c
//     if (diagnostics->information_content) {
//         double* info = _PSICalculateInformationContentFromFreqRatios(
//                 internal_pssm->freq_ratios, seq_weights->std_prob,
//                 diagnostics->query_length,
//                 diagnostics->alphabet_size);
//         ...
//     if (diagnostics->weighted_residue_freqs) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             for (r = 0; r < diagnostics->alphabet_size; r++) {
//                 diagnostics->weighted_residue_freqs[p][r] =
//                     seq_weights->match_weights[p][r];
//             }
//         }
//     }
//
//     if (diagnostics->frequency_ratios) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             for (r = 0; r < diagnostics->alphabet_size; r++) {
//                 diagnostics->frequency_ratios[p][r] =
//                     internal_pssm->freq_ratios[p][r];
//             }
//         }
//     }
//
//     if (diagnostics->gapless_column_weights) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             if (msa->num_matching_seqs[p] > 1
//                 && msa->cell[0][p].letter != kXResidue) {
//
//                 diagnostics->gapless_column_weights[p] =
//                     seq_weights->gapless_column_weights[p]
//                     / internal_pssm->pseudocounts[p];
//
//                 diagnostics->gapless_column_weights[p] *=
//                     (seq_weights->sigma[p] / aligned_block->size[p] - 1);
//             }
//             else {
//                 diagnostics->gapless_column_weights[p] = 0.0;
//             }
//         }
//     }
//
//     if (diagnostics->sigma) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             diagnostics->sigma[p] = seq_weights->sigma[p];
//         }
//     }
//
//     if (diagnostics->interval_sizes) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             diagnostics->interval_sizes[p] = aligned_block->size[p];
//         }
//     }
//
//     if (diagnostics->num_matching_seqs) {
//         for (p = 0; p < diagnostics->query_length; p++) {
//             diagnostics->num_matching_seqs[p] = msa->num_matching_seqs[p];
//         }
//     }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_psi_priv.c:2360-2378
// ```c
//     for (p = 0; p < query_length; p++) {
//
//         double info_sum = 0.0;
//
//         for (r = 0; r < alphabet_sz; r++) {
//
//             if (std_prob[r] > kEpsilon) {
//                 /* Division compensates for multiplication in
//                  * _PSIComputeFreqRatios */
//                 double qOverPEstimate = freq_ratios[p][r] / std_prob[r];
//                 if (qOverPEstimate > kEpsilon) {
//                     info_sum +=
//                         freq_ratios[p][r] * log(qOverPEstimate) / NCBIMATH_LN2;
//                 }
//             }
//
//         }
//         retval[p] = info_sum;
//     }
// ```
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
                            info_sum += ratios[r] * q_over_p_estimate.ln() / NCBIMATH_LN2;
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
