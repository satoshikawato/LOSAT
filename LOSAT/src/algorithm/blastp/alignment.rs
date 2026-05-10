//! BLASTP alignment reconstruction helpers.
//!
//! Reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c

use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, protein_score};

use super::encoding::ncbistdaa_to_ascii;

#[inline(always)]
fn blastp_positive_score<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query_residue: u8,
    subject_residue: u8,
) -> bool {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:745-818
    // ```c
    // if (*q == *s)
    //    num_ident++;
    // else if (matrix[*q][*s] > 0)
    //    num_pos ++;
    // ```
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue) > 0
    } else {
        protein_score(matrix, query_residue, subject_residue) > 0
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
// ```c
// static Int2
// s_Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query,
//                                         const Uint1* subject,
//                                         const BlastHSP* hsp,
//                                         Int4* num_ident_ptr,
//                                         Int4* align_length_ptr,
//                                         const BlastScoreBlk* sbp,
//                                         Int4* num_pos_ptr)
// {
//    ...
//    if (*q == *s)
//        num_ident++;
//    else if (matrix[*q][*s] > 0)
//        num_pos ++;
//    ...
//    *num_pos_ptr = num_pos + num_ident;
// }
// ```
#[derive(Debug, Clone)]
pub struct ProteinAlignmentView {
    pub query_aln: String,
    pub subject_aln: String,
    pub num_ident: usize,
    pub num_positives: usize,
    pub gap_letters: usize,
    pub gap_opens: usize,
    pub alignment_len: usize,
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
// ```c
// if (!hsp->gap_info) {
//    /* Ungapped case. Check that lengths are the same in query and subject,
//       then count number of matches. */
//    if (q_length != s_length)
//       return -1;
//    align_length = q_length;
//    for (i=0; i<align_length; i++) {
//       if (*q == *s)
//          num_ident++;
//       else if (matrix[*q][*s] > 0)
//          num_pos ++;
//       q++;
//       s++;
//    }
// } else {
//    ...
// }
// ```
pub fn build_alignment_view_with_matrix(
    query: &[u8],
    subject: &[u8],
    hit: &Hit,
    matrix: ScoringMatrix,
) -> Option<ProteinAlignmentView> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
    // ```c
    // if (*q == *s)
    //    num_ident++;
    // else if (matrix[*q][*s] > 0)
    //    num_pos ++;
    // ```
    if matrix == ScoringMatrix::Blosum62 {
        build_alignment_view_with_matrix_impl::<true>(query, subject, hit, matrix)
    } else {
        build_alignment_view_with_matrix_impl::<false>(query, subject, hit, matrix)
    }
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
// ```c
// static Int2
// s_Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query,
//                                         const Uint1* subject,
//                                         const BlastHSP* hsp,
//                                         Int4* num_ident_ptr,
//                                         Int4* align_length_ptr,
//                                         const BlastScoreBlk* sbp,
//                                         Int4* num_pos_ptr)
// {
//    ...
// }
// ```
fn build_alignment_view_with_matrix_impl<const BLOSUM62: bool>(
    query: &[u8],
    subject: &[u8],
    hit: &Hit,
    matrix: ScoringMatrix,
) -> Option<ProteinAlignmentView> {
    let q_offset = hit.q_start.checked_sub(1)?;
    let s_offset = hit.s_start.min(hit.s_end).checked_sub(1)?;

    let implicit_ungapped;
    let ops: &[GapEditOp] = if let Some(ref gap_info) = hit.gap_info {
        gap_info.as_slice()
    } else {
        let q_span = hit.q_end.checked_sub(hit.q_start)?.saturating_add(1);
        let s_span = hit
            .s_end
            .max(hit.s_start)
            .checked_sub(hit.s_start.min(hit.s_end))?
            .saturating_add(1);
        if q_span != hit.length || s_span != hit.length {
            return None;
        }
        implicit_ungapped = [GapEditOp::Sub(hit.length as u32)];
        &implicit_ungapped
    };

    let mut q_index = q_offset;
    let mut s_index = s_offset;
    let mut query_aln = String::new();
    let mut subject_aln = String::new();
    let mut num_ident = 0usize;
    let mut num_positives = 0usize;
    let mut gap_letters = 0usize;
    let mut gap_opens = 0usize;
    let mut alignment_len = 0usize;

    for op in ops {
        match *op {
            GapEditOp::Sub(n) => {
                for _ in 0..n as usize {
                    let q = *query.get(q_index)?;
                    let s = *subject.get(s_index)?;
                    query_aln.push(ncbistdaa_to_ascii(q));
                    subject_aln.push(ncbistdaa_to_ascii(s));
                    if q == s {
                        num_ident += 1;
                        num_positives += 1;
                    } else if blastp_positive_score::<BLOSUM62>(matrix, q, s) {
                        num_positives += 1;
                    }
                    q_index += 1;
                    s_index += 1;
                    alignment_len += 1;
                }
            }
            GapEditOp::Del(n) => {
                gap_opens += 1;
                for _ in 0..n as usize {
                    let s = *subject.get(s_index)?;
                    query_aln.push('-');
                    subject_aln.push(ncbistdaa_to_ascii(s));
                    s_index += 1;
                    gap_letters += 1;
                    alignment_len += 1;
                }
            }
            GapEditOp::Ins(n) => {
                gap_opens += 1;
                for _ in 0..n as usize {
                    let q = *query.get(q_index)?;
                    query_aln.push(ncbistdaa_to_ascii(q));
                    subject_aln.push('-');
                    q_index += 1;
                    gap_letters += 1;
                    alignment_len += 1;
                }
            }
        }
    }

    Some(ProteinAlignmentView {
        query_aln,
        subject_aln,
        num_ident,
        num_positives,
        gap_letters,
        gap_opens,
        alignment_len,
    })
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:746-824
// ```c
// static Int2
// s_Blast_HSPGetNumIdentitiesAndPositives(..., const BlastScoreBlk* sbp, ...)
// {
//     ...
//     else if (matrix[*q][*s] > 0)
//         num_pos ++;
// }
// ```
pub fn build_alignment_view(
    query: &[u8],
    subject: &[u8],
    hit: &Hit,
) -> Option<ProteinAlignmentView> {
    build_alignment_view_with_matrix(query, subject, hit, ScoringMatrix::Blosum62)
}

#[cfg(test)]
mod tests {
    use crate::utils::matrix::ncbistdaa;

    use super::*;

    #[test]
    fn test_build_alignment_view_ungapped() {
        let hit = Hit {
            identity: 75.0,
            length: 4,
            mismatch: 1,
            gapopen: 0,
            q_start: 1,
            q_end: 4,
            s_start: 1,
            s_end: 4,
            e_value: 1e-10,
            bit_score: 50.0,
            num_ident: 3,
            query_frame: 0,
            query_length: 4,
            q_idx: 0,
            s_idx: 0,
            raw_score: 20,
            gap_info: None,
            num_positives: 4,
        };

        let query = [1u8, 2, 3, 4];
        let subject = [1u8, 2, 5, 4];
        let view = build_alignment_view(&query, &subject, &hit).expect("alignment");
        assert_eq!(view.alignment_len, 4);
        assert_eq!(view.num_ident, 3);
        assert_eq!(view.gap_letters, 0);
    }

    #[test]
    fn test_build_alignment_view_gapped() {
        let hit = Hit {
            identity: 75.0,
            length: 4,
            mismatch: 0,
            gapopen: 1,
            q_start: 1,
            q_end: 3,
            s_start: 1,
            s_end: 4,
            e_value: 1e-10,
            bit_score: 50.0,
            num_ident: 3,
            query_frame: 0,
            query_length: 3,
            q_idx: 0,
            s_idx: 0,
            raw_score: 20,
            gap_info: Some(vec![
                GapEditOp::Sub(2),
                GapEditOp::Del(1),
                GapEditOp::Sub(1),
            ]),
            num_positives: 3,
        };

        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_encoding.c:115-121
        // ```c
        // const char NCBISTDAA_TO_AMINOACID[BLASTAA_SIZE] = {
        // '-','A','B','C','D','E','F','G','H','I','K','L','M',
        // 'N','P','Q','R','S','T','V','W','X','Y','Z','U','*',
        // 'O', 'J'};
        // ```
        let query = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::N];
        let subject = [ncbistdaa::A, ncbistdaa::R, ncbistdaa::D, ncbistdaa::C];
        let view = build_alignment_view(&query, &subject, &hit).expect("alignment");
        // NCBI reference: ncbi-blast/c++/include/algo/blast/core/gapinfo.h:45-51
        // ```c
        // eGapAlignDel = 0, /**< Deletion: a gap in query */
        // ...
        // eGapAlignIns = 6, /**< Insertion: a gap in subject */
        // ```
        assert_eq!(view.query_aln, "AR-N");
        assert_eq!(view.subject_aln, "ARDC");
        assert_eq!(view.gap_letters, 1);
        assert_eq!(view.gap_opens, 1);
    }
}
