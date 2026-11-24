//! Tests of `delta_hooks` (series PR, item PS13a, fix round 1).

use crate::algorithm::pssm_search::in_pssm::query_ungapped_block;
use crate::config::ScoringMatrix;

/// `DeltaSide::cdd_pssm` refuses a query whose letters give a Karlin-Altschul block with the domain
/// search's matrix and none with `-matrix`, the case where NCBI's engine would read a null block.
/// With the matrices LOSAT knows (NCBI's), whether a query has the block depends on its set of
/// letters only: each letter gives a block under every matrix or under none, and every letter's
/// mean score against the standard composition is negative under every matrix, so a mix of letters
/// gives a block wherever one of its letters does. So the refusal is not reached; this test (every
/// letter, and every pair of letters) keeps it so if a matrix is added.
#[test]
fn every_letter_set_gives_a_karlin_block_under_every_matrix_or_none() {
    let matrices: Vec<ScoringMatrix> = ScoringMatrix::ALL
        .into_iter()
        .filter(|&matrix| matrix != ScoringMatrix::Identity)
        .collect();
    for first in 1u8..28 {
        for second in first..28 {
            let mut query = vec![first; 20];
            query.extend(std::iter::repeat_n(second, 20));
            let blocks: Vec<bool> = matrices
                .iter()
                .map(|&matrix| query_ungapped_block(matrix, &query).is_some())
                .collect();
            assert!(
                blocks.iter().all(|&block| block == blocks[0]),
                "letters {first} and {second}: {blocks:?} under {matrices:?}"
            );
        }
    }
}
