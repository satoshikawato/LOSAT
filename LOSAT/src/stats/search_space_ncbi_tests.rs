//! NCBI BLAST comparison of the effective search space (was tests/unit_tests.rs, mod ncbi_comparison),
//! through the one owner of `BLAST_CalcEffLengths` (`stats::setup::calc_eff_lengths`).

use crate::stats::tables::KarlinParams;

/// Test case from actual NCBI BLAST run
///
/// Query: test_sequences/query.fasta (6188 bp)
/// Subject: test_sequences/subject.fasta (6188 bp)
/// Parameters: reward=1, penalty=-2, gapopen=0, gapextend=0
/// NCBI BLAST reported: "Effective search space used: 38081241"
#[test]
fn test_effective_search_space_ncbi_comparison() {
    // BLASTN parameters for reward=1, penalty=-2, gapopen=0, gapextend=0
    let params = KarlinParams {
        lambda: 1.28,
        k: 0.46,
        h: 0.85,
        alpha: 1.5,
        beta: -2.0,
    };
    let lengths =
        crate::stats::setup::calc_eff_lengths(6188, 6188, 1, &params, 0, "BLASTN").unwrap();
    assert_eq!(lengths.eff_searchsp, 38_081_241);
    assert!(lengths.length_adjustment >= 0 && lengths.length_adjustment < 6188);
}
