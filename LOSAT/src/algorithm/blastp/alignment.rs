//! BLASTP alignment reconstruction helpers.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c

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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-818
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue) > 0
    } else {
        protein_score(matrix, query_residue, subject_residue) > 0
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-824
#[derive(Debug, Clone)]
#[cfg_attr(
    not(test),
    expect(
        dead_code,
        reason = "num_ident, num_positives, gap_letters, gap_opens and alignment_len are read only by unit tests"
    )
)]
pub struct ProteinAlignmentView {
    pub query_aln: String,
    pub subject_aln: String,
    pub num_ident: usize,
    pub num_positives: usize,
    pub gap_letters: usize,
    pub gap_opens: usize,
    pub alignment_len: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
pub fn build_alignment_view_with_matrix(
    query: &[u8],
    subject: &[u8],
    hit: &Hit,
    matrix: ScoringMatrix,
) -> Option<ProteinAlignmentView> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
    if matrix == ScoringMatrix::Blosum62 {
        build_alignment_view_with_matrix_impl::<true>(query, subject, hit, matrix)
    } else {
        build_alignment_view_with_matrix_impl::<false>(query, subject, hit, matrix)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
fn build_alignment_view_with_matrix_impl<const BLOSUM62: bool>(
    query: &[u8],
    subject: &[u8],
    hit: &Hit,
    matrix: ScoringMatrix,
) -> Option<ProteinAlignmentView> {
    // An ungapped HSP of length 0 (`-ungapped` with a low `-threshold`, which `Blast_HSPTest` keeps)
    // has no column: empty rows, whose range of the subject may start at its first letter.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:773-780
    if hit.length == 0 {
        return Some(ProteinAlignmentView {
            query_aln: String::new(),
            subject_aln: String::new(),
            num_ident: 0,
            num_positives: 0,
            gap_letters: 0,
            gap_opens: 0,
            alignment_len: 0,
        });
    }
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
    // The rows live in a `PairwiseHit` while the reports are written: their exact size.
    let columns = ops.iter().map(|op| op.num() as usize).sum();
    let mut query_aln = String::with_capacity(columns);
    let mut subject_aln = String::with_capacity(columns);
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

#[cfg(test)]
#[path = "alignment_tests.rs"]
mod tests;
