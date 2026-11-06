//! Ports of NCBI's PSSM engine unit tests (`pssmcreate_unit_test.cpp`) for the stages of
//! `stages.rs` and the engine (`create_pssm_with_diagnostics`).

use super::*;
use crate::config::ScoringMatrix;
use crate::core::psi::{
    create_pssm_with_diagnostics, KarlinBlk, PsiDiagnosticsRequest, PsiMsa, PsiMsaCell, PsiScoreBlk,
};
use crate::utils::matrix::{protein_score, BLAST_SCORE_MIN};

const PURGE_TEST_QUERY_LENGTH: usize = 232;

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:549-565
// (`CPssmInputTestData::kQuery`, Gi 129295)
const PURGE_TEST_QUERY: [u8; PURGE_TEST_QUERY_LENGTH] = [
    15, 9, 10, 4, 11, 11, 19, 17, 17, 17, 18, 4, 11, 4, 18, 18, 11, 19, 11, 19, 13, 1, 9, 22, 6,
    10, 7, 12, 20, 10, 18, 1, 6, 13, 1, 5, 4, 18, 16, 5, 12, 14, 6, 8, 19, 18, 10, 15, 5, 17, 10,
    14, 19, 15, 12, 12, 3, 12, 13, 13, 17, 6, 13, 19, 1, 18, 11, 14, 1, 5, 10, 12, 10, 9, 11, 5,
    11, 14, 6, 1, 17, 7, 4, 11, 17, 12, 11, 19, 11, 11, 14, 4, 5, 19, 17, 4, 11, 5, 16, 9, 5, 10,
    18, 9, 13, 6, 5, 10, 11, 18, 5, 20, 18, 13, 14, 13, 18, 12, 5, 10, 16, 16, 19, 10, 19, 22, 11,
    14, 15, 12, 10, 9, 5, 5, 10, 22, 13, 11, 18, 17, 19, 11, 12, 1, 11, 7, 12, 18, 4, 11, 6, 9, 14,
    17, 1, 13, 11, 18, 7, 9, 17, 17, 1, 5, 17, 11, 10, 9, 17, 15, 1, 19, 8, 7, 1, 6, 12, 5, 11, 17,
    5, 4, 7, 9, 5, 12, 1, 7, 17, 18, 7, 19, 9, 5, 4, 9, 10, 8, 17, 14, 5, 17, 5, 15, 6, 16, 1, 4,
    8, 14, 6, 11, 6, 11, 9, 10, 8, 13, 14, 18, 13, 18, 9, 19, 22, 6, 7, 16, 22, 20, 17, 14,
];

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:487-513
// (`SetupDuplicateHit`'s `kGi_129296_`)
const GI_129296: [u8; 388] = [
    12, 4, 17, 9, 17, 19, 18, 13, 1, 10, 6, 3, 6, 4, 19, 6, 13, 5, 12, 10, 19, 8, 8, 19, 13, 5, 13,
    9, 11, 22, 3, 14, 11, 17, 9, 11, 18, 1, 11, 1, 12, 19, 22, 11, 7, 1, 16, 7, 13, 18, 5, 17, 15,
    12, 10, 10, 19, 11, 8, 6, 4, 17, 9, 18, 7, 1, 7, 17, 18, 18, 4, 17, 15, 3, 7, 17, 17, 5, 22,
    19, 8, 13, 11, 6, 10, 5, 11, 11, 17, 5, 9, 18, 16, 14, 13, 1, 18, 22, 17, 11, 5, 9, 1, 4, 10,
    11, 22, 19, 4, 10, 18, 6, 17, 19, 11, 14, 5, 22, 11, 17, 3, 1, 16, 10, 6, 22, 18, 7, 7, 19, 5,
    5, 19, 13, 6, 10, 18, 1, 1, 5, 5, 1, 16, 15, 11, 9, 13, 17, 20, 19, 5, 10, 5, 18, 13, 7, 15, 9,
    10, 4, 11, 11, 19, 17, 17, 17, 9, 4, 6, 7, 18, 18, 12, 19, 6, 9, 13, 18, 9, 22, 6, 10, 7, 9,
    20, 10, 9, 1, 6, 13, 18, 5, 4, 18, 16, 5, 12, 14, 6, 17, 12, 18, 10, 5, 5, 17, 10, 14, 19, 15,
    12, 12, 3, 12, 13, 13, 17, 6, 13, 19, 1, 18, 11, 14, 1, 5, 10, 12, 10, 9, 11, 5, 11, 14, 22, 1,
    17, 7, 4, 11, 17, 12, 11, 19, 11, 11, 14, 4, 5, 19, 17, 7, 11, 5, 16, 9, 5, 10, 18, 9, 13, 6,
    4, 10, 11, 16, 5, 20, 18, 17, 18, 13, 1, 12, 1, 10, 10, 17, 12, 10, 19, 22, 11, 14, 16, 12, 10,
    9, 5, 5, 10, 22, 13, 11, 18, 17, 9, 11, 12, 1, 11, 7, 12, 18, 4, 11, 6, 17, 16, 17, 1, 13, 11,
    18, 7, 9, 17, 17, 19, 4, 13, 11, 12, 9, 17, 4, 1, 19, 8, 7, 19, 6, 12, 5, 19, 13, 5, 5, 7, 18,
    5, 1, 18, 7, 17, 18, 7, 1, 9, 7, 13, 9, 10, 8, 17, 11, 5, 11, 5, 5, 6, 16, 1, 4, 8, 14, 6, 11,
    6, 6, 9, 16, 22, 13, 14, 18, 13, 1, 9, 11, 6, 6, 7, 16, 22, 20, 17, 14,
];

/// BLOSUM62 11/1's gapped block (`Blast_KarlinBlkGappedCalc`).
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:268
// ```c
//     {11, 1, (double) INT2_MAX, 0.267, 0.041, 0.14, 1.9, -30, 0.669720, 42.602800, 43.636200},
// ```
const BLOSUM62_11_1: KarlinBlk = KarlinBlk {
    lambda: 0.267,
    k: 0.041,
    h: 0.14,
};

fn blosum62() -> PsiScoreBlk {
    crate::algorithm::psiblast::pssm::score_block(ScoringMatrix::Blosum62, BLOSUM62_11_1).unwrap()
}

/// psiblast's request without `-out_ascii_pssm` (`PSIDiagnosticsRequestNewEx(FALSE)`).
const REQUEST: PsiDiagnosticsRequest = PsiDiagnosticsRequest {
    information_content: false,
    weighted_residue_frequencies: false,
    frequency_ratios: true,
    gapless_column_weights: false,
    sigma: false,
    interval_sizes: false,
    num_matching_seqs: false,
};

/// An MSA whose rows hold `rows` (every cell aligned), the first of them the query.
fn aligned_msa(rows: &[&[u8]]) -> PsiMsa {
    let mut msa = PsiMsa::new(rows[0].len(), rows.len() - 1);
    for (cells, letters) in msa.rows.iter_mut().zip(rows) {
        for (cell, &letter) in cells.iter_mut().zip(letters.iter()) {
            *cell = PsiMsaCell {
                letter,
                is_aligned: true,
            };
        }
    }
    msa
}

/// `SetupDuplicateHit`: the query and Gi 129296's first letters twice.
fn duplicate_hit_msa() -> PsiMsa {
    let hit = &GI_129296[..PURGE_TEST_QUERY_LENGTH];
    aligned_msa(&[&PURGE_TEST_QUERY, hit, hit])
}

/// `FindNonIdenticalHighScoringResidue`: the residue other than `res` with the highest BLOSUM62
/// score (the first of equal ones).
fn non_identical_high_scoring_residue(res: u8) -> u8 {
    let mut retval = ncbistdaa::GAP;
    let mut max_score = BLAST_SCORE_MIN;
    for i in 0..BLASTAA_SIZE as u8 {
        if i == res {
            continue;
        }
        let score = protein_score(ScoringMatrix::Blosum62, res, i);
        if score > max_score {
            max_score = score;
            retval = i;
        }
    }
    assert_ne!(retval, ncbistdaa::GAP);
    retval
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:762-772 testPurgeSelfHit
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:762-772
// ```c++
//         int rv = _PSIPurgeBiasedSegments(msa.get());
//         BOOST_REQUIRE_EQUAL(PSI_SUCCESS, rv);
//         const Uint4 kSelfHitIndex = 1;
// 		BOOST_REQUIRE_EQUAL(true, !!msa->use_sequence[kQueryIndex]);
// 		BOOST_REQUIRE_EQUAL(false, !!msa->use_sequence[kSelfHitIndex]);
// ```
#[test]
fn purge_removes_a_self_hit() {
    let msa = aligned_msa(&[&PURGE_TEST_QUERY, &PURGE_TEST_QUERY]);
    let mut packed = PackedMsa::new(&msa);
    purge_biased_segments(&mut packed);
    assert!(packed.use_sequence[K_QUERY_INDEX]);
    assert!(!packed.use_sequence[1]);
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:774-785 testPurgeDuplicateHit
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:781-784
// ```c++
//         const Uint4 kDuplicateHitIndex = 2;
//         BOOST_REQUIRE_EQUAL(false, !!msa->use_sequence[kDuplicateHitIndex]);
//         BOOST_REQUIRE_EQUAL(true, !!msa->use_sequence[kQueryIndex]);
//         BOOST_REQUIRE_EQUAL(true, !!msa->use_sequence[kQueryIndex + 1]);
// ```
#[test]
fn purge_removes_a_duplicate_hit() {
    let mut packed = PackedMsa::new(&duplicate_hit_msa());
    purge_biased_segments(&mut packed);
    assert!(!packed.use_sequence[2]);
    assert!(packed.use_sequence[K_QUERY_INDEX]);
    assert!(packed.use_sequence[K_QUERY_INDEX + 1]);
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:787-799 testPurgeNearIdenticalHits
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:533-545
// ```c++
//         SetupDuplicateHit();
//
//         const Uint4 kHitIndex = 2;  // index of the near identical hit
//         const Uint4 kNumIdenticalResidues = (Uint4) (GetQueryLength() *
//             (kPSINearIdentical + 0.01));
//
//         for (Uint4 i = kNumIdenticalResidues; i < GetQueryLength(); i++) {
//             Uint1& residue = m_msa->data[kHitIndex][i].letter;
//             residue = (residue + 1) % BLASTAA_SIZE;
// ```
#[test]
fn purge_removes_a_near_identical_hit() {
    let mut msa = duplicate_hit_msa();
    let identical = (PURGE_TEST_QUERY_LENGTH as f64 * (K_PSI_NEAR_IDENTICAL + 0.01)) as usize;
    for cell in &mut msa.rows[2][identical..] {
        cell.letter = (cell.letter + 1) % BLASTAA_SIZE as u8;
        assert!(cell.letter > 0);
    }
    let mut packed = PackedMsa::new(&msa);
    purge_biased_segments(&mut packed);
    assert!(!packed.use_sequence[2]);
    assert!(packed.use_sequence[K_QUERY_INDEX]);
    assert!(packed.use_sequence[K_QUERY_INDEX + 1]);
}

/// Whether the residue `j` of a column is one that the PSSM gives `BLAST_SCORE_MIN` (the gap,
/// B, Z, U and the codes from O on).
fn exceptional_residue(j: u8) -> bool {
    j == ncbistdaa::GAP
        || j == ncbistdaa::B
        || j == ncbistdaa::Z
        || j == ncbistdaa::U
        || j >= ncbistdaa::O
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:801-865 testQueryAlignedWithInternalGaps
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:401-436
// ```c++
//         string query_seq("MFKVYGYDSNIHKCGPCDNAKRLLTVKKQPFEFINIM");
//         query_seq += string("PEKGVFDDEKIAELLTKLGRDTQIGLTMPQVFAPDGSHIGGFD");
//         query_seq += string("QLREYFK");
//         ...
//         TAlignedSequence region(make_pair(make_pair(0U, 8U),
//                                           string("KVVVFIKP")));
//         ...
//             // align the second sequence to gaps
//             m_msa->data[1][i].letter = kGapResidue;
//             m_msa->data[1][i].is_aligned = true;
// ```
#[test]
fn a_query_aligned_with_internal_gaps_keeps_the_matrix_scores_at_the_gaps() {
    let letters = |text: &str| -> Vec<u8> {
        text.bytes()
            .map(crate::utils::matrix::aa_char_to_ncbistdaa)
            .collect()
    };
    let query = letters(concat!(
        "MFKVYGYDSNIHKCGPCDNAKRLLTVKKQPFEFINIM",
        "PEKGVFDDEKIAELLTKLGRDTQIGLTMPQVFAPDGSHIGGFD",
        "QLREYFK"
    ));
    assert_eq!(query.len(), 87);
    let mut aligned = vec![ncbistdaa::GAP; query.len()];
    for (start, text) in [
        (0, "KVVVFIKP"),
        (12, "TCPFCRKTQELLSQLPFLLEFVDITAT"),
        (41, "SDTNEIQDYLQQLTGA"),
        (62, "RTVPRVFIG"),
        (72, "KECIGGCTDLESMHK"),
    ] {
        let region = letters(text);
        aligned[start..start + region.len()].copy_from_slice(&region);
    }
    let msa = aligned_msa(&[&query, &aligned]);
    let pssm = create_pssm_with_diagnostics(&msa, 0, &blosum62(), REQUEST).unwrap();
    assert_eq!(pssm.scores.len(), query.len());
    for (i, scores) in pssm.scores.iter().enumerate() {
        if aligned[i] != ncbistdaa::GAP {
            continue;
        }
        for (j, &score) in scores.iter().enumerate() {
            if exceptional_residue(j as u8) {
                assert_eq!(score, BLAST_SCORE_MIN, "position {i} residue {j}");
            } else {
                let expected = protein_score(ScoringMatrix::Blosum62, query[i], j as u8);
                assert!(
                    expected - score <= 3,
                    "position {i} residue {j}: expected={expected} actual={score}"
                );
            }
        }
    }
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:867-1061 testMultiSeqAlignmentHasRegionsUnalignedToQuery
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:367-385
// ```c++
//         const TAlignedSegment kFirstAlignment(0, 100);
//         for (unsigned int i = kFirstAlignment.first;
//              i < kFirstAlignment.second; i++) {
//             m_msa->data[1][i].letter =
//                 FindNonIdenticalHighScoringResidue(kQuery[i], score_matrix);
//             m_msa->data[1][i].is_aligned = true;
//         }
//         ...
//         const TAlignedSegment kSecondAlignment(200, kQueryLength);
// ```
#[test]
fn a_region_unaligned_to_the_query_keeps_the_matrix_scores() {
    let length = PURGE_TEST_QUERY_LENGTH;
    let mut input = PsiMsa::new(length, 2);
    for (p, &letter) in PURGE_TEST_QUERY.iter().enumerate() {
        input.rows[0][p] = PsiMsaCell {
            letter,
            is_aligned: true,
        };
    }
    for (row, range) in [(1, 0..100), (2, 200..length)] {
        for p in range {
            input.rows[row][p] = PsiMsaCell {
                letter: non_identical_high_scoring_residue(PURGE_TEST_QUERY[p]),
                is_aligned: true,
            };
        }
    }
    let sbp = blosum62();

    let mut packed = PackedMsa::new(&input);
    purge_biased_segments(&mut packed);
    assert!(packed.use_sequence[K_QUERY_INDEX]);
    assert!(packed.use_sequence[1]);
    assert!(packed.use_sequence[2]);

    let mut msa = Msa::new(&packed);
    let mut aligned_blocks = AlignedBlock::new(length);
    compute_alignment_blocks(&mut msa, &mut aligned_blocks);
    for (first, second) in [(0, 99), (200, length - 1)] {
        for pos in first..second {
            assert_eq!(
                aligned_blocks.pos_extnt[pos].0, first as i32,
                "position {pos}"
            );
            assert_eq!(
                aligned_blocks.pos_extnt[pos].1, second as i32,
                "position {pos}"
            );
            assert_eq!(
                aligned_blocks.size[pos],
                (second - first + 1) as u32,
                "position {pos}"
            );
        }
    }
    for pos in 100..200 {
        let (left, right) = aligned_blocks.pos_extnt[pos];
        assert_eq!(left, -1, "position {pos}");
        assert_eq!(right, length as i32, "position {pos}");
        assert_eq!(
            (right - left + 1) as u32,
            aligned_blocks.size[pos],
            "position {pos}"
        );
    }

    let mut seq_weights = SequenceWeights::new(msa.num_seqs, length);
    compute_sequence_weights(&msa, &aligned_blocks, &mut seq_weights).unwrap();
    let column_weight = |pos: usize| -> f64 {
        (0..BLASTAA_SIZE)
            .filter(|&res| res != ncbistdaa::X as usize)
            .map(|res| seq_weights.match_weights[pos][res])
            .sum()
    };
    for (first, second) in [(0, 99), (200, length - 1)] {
        for pos in first..second {
            let total = column_weight(pos);
            assert!(total > 0.99 && total < 1.01, "position {pos}: {total}");
        }
    }
    for pos in 100..200 {
        assert_eq!(column_weight(pos), 0.0, "position {pos}");
    }

    let mut internal = InternalPssm::new(length);
    compute_freq_ratios(
        &msa,
        &mut seq_weights,
        &sbp,
        &aligned_blocks,
        0,
        &mut internal,
    );
    convert_freq_ratios_to_pssm(&mut internal, &msa.query, &sbp, &seq_weights.std_prob);
    scale_matrix(&msa.query, &seq_weights.std_prob, &mut internal, &sbp).unwrap();
    assert_eq!(msa.num_seqs, 3);
    for i in 0..length {
        if (1..=3).any(|s| msa.cell[s][i].is_aligned) {
            continue;
        }
        for j in 0..BLASTAA_SIZE as u8 {
            if j != ncbistdaa::GAP && exceptional_residue(j) {
                continue;
            }
            let actual = internal.pssm[i][j as usize];
            if j == ncbistdaa::GAP {
                assert_eq!(actual, BLAST_SCORE_MIN, "position {i} residue {j}");
            } else {
                let score = protein_score(ScoringMatrix::Blosum62, msa.query[i], j);
                assert!(
                    score - 1 <= actual && actual <= score + 1,
                    "position {i} residue {j}: expected={score} actual={actual}"
                );
            }
        }
    }
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1066-1110 testQueryIsOnlyAlignedSequenceInMsa
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1097-1104
// ```c++
//         for (size_t i = 0; i < pssm_input->GetQueryLength(); i++) {
//             BOOST_REQUIRE_EQUAL((int)-1,
//                                  (int)aligned_blocks->pos_extnt[i].left);
//             BOOST_REQUIRE_EQUAL((int)pssm_input->GetQueryLength(),
//                                  (int)aligned_blocks->pos_extnt[i].right);
//             BOOST_REQUIRE_EQUAL((int)pssm_input->GetQueryLength() + 2,
//                                  (int)aligned_blocks->size[i]);
//         }
// ```
// The rest of the case computes the weights in the structure group compatibility mode, which
// psiblast does not use and `stages.rs` does not port.
#[test]
fn the_query_alone_has_no_aligned_blocks() {
    let msa = aligned_msa(&[&PURGE_TEST_QUERY, &PURGE_TEST_QUERY]);
    let mut packed = PackedMsa::new(&msa);
    purge_biased_segments(&mut packed);
    assert!(packed.use_sequence[K_QUERY_INDEX]);
    assert!(!packed.use_sequence[1]);
    let mut msa = Msa::new(&packed);
    let mut aligned_blocks = AlignedBlock::new(PURGE_TEST_QUERY_LENGTH);
    compute_alignment_blocks(&mut msa, &mut aligned_blocks);
    for i in 0..PURGE_TEST_QUERY_LENGTH {
        assert_eq!(aligned_blocks.pos_extnt[i].0, -1);
        assert_eq!(
            aligned_blocks.pos_extnt[i].1,
            PURGE_TEST_QUERY_LENGTH as i32
        );
        assert_eq!(aligned_blocks.size[i], PURGE_TEST_QUERY_LENGTH as u32 + 2);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:108-131
// ```c++
//         const unsigned int  kQuerySize = 10;
//         const unsigned int  kNumSeqs = 2;
//         const unsigned char kQuery[] = { 3, 9, 14, 20, 6, 23, 1, 7, 16, 5 };
//         ...
//         // Add the flanking gaps
//         m_msa->data[1][0].letter =
//             m_msa->data[2][0].letter =
//             m_msa->data[2][m_dim.query_length-1].letter =
//             AMINOACID_TO_NCBISTDAA[(int)'-'];
// ```
const FLANKING_GAPS_QUERY: [u8; 10] = [3, 9, 14, 20, 6, 23, 1, 7, 16, 5];

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1196-1200 testRejectFlankingGaps
#[test]
fn flanking_gaps_are_rejected() {
    let mut msa = aligned_msa(&[&FLANKING_GAPS_QUERY[..]; 3]);
    msa.rows[1][0].letter = ncbistdaa::GAP;
    msa.rows[2][0].letter = ncbistdaa::GAP;
    msa.rows[2][9].letter = ncbistdaa::GAP;
    assert!(create_pssm_with_diagnostics(&msa, 0, &blosum62(), REQUEST).is_err());
}

// NCBI unit test (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:1202-1206 testRejectGapInQuery
// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/pssmcreate_unit_test.cpp:178-182
// ```c++
//         // Randomly assign a position in the query to contain a gap
//         CRandom r((CRandom::TValue)time(0));
//         int gap_position = r.GetRand(0, GetQueryLength() - 1);
//         m_query[gap_position] = AMINOACID_TO_NCBISTDAA[(int)'-'];
//         m_msa->data[0][gap_position].letter = m_query[gap_position];
// ```
// Every position in turn instead of a random one.
#[test]
fn a_gap_in_the_query_is_rejected() {
    for gap_position in 0..FLANKING_GAPS_QUERY.len() {
        let mut msa = aligned_msa(&[&FLANKING_GAPS_QUERY[..]; 3]);
        msa.rows[0][gap_position].letter = ncbistdaa::GAP;
        assert!(
            create_pssm_with_diagnostics(&msa, 0, &blosum62(), REQUEST).is_err(),
            "gap at {gap_position}"
        );
    }
}
