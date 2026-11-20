//! The letters of an alignment's rows as NCBI's reports read and compare them: the query's
//! letters (`QueryLetters`: a FASTA record's, or the Bioseq's of a PSSM query), how a query row's
//! lower case reads (`QueryCase`), the positives of two rows (`display_positives`), and the
//! columns whose `-` is a residue (`stored_gap_columns`).
//!
//! NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-121
//! ```c++
//!         CBioseq_Handle h = GetBioseqHandle(row);
//!         CSeqVector vec = h.GetSeqVector
//!             (CBioseq_Handle::eCoding_Iupac,
//!              IsPositiveStrand(row) ?
//!              CBioseq_Handle::eStrand_Plus :
//!              CBioseq_Handle::eStrand_Minus);
//! ```
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
//! ```c++
//!                     } else if (m_SeqLocChar==eLowerCase){
//!                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
//! ```
//! NCBI's reports read the letters of every row from its Bioseq and compare them as they are; a
//! mask lowers a letter only where outfmt 0 prints the row. A FASTA record's Bioseq holds its
//! letters in upper case (NCBI's reader keeps the lower case as a mask, fasta.cpp:945-949), and
//! LOSAT's outfmt 0 rows carry those masks in lower case, so the reports compare the capitals of
//! such a row (`QueryCase::Masked`). The Bioseq of a PSSM query (`-in_pssm`) holds iupacaa letters
//! as the file has them, in lower case too (`pssm_search::PsiSearch::bioseq_letters`), which the
//! search reads as their capitals; its rows show no mask (NCBI's PSSM search sets no filtered
//! regions: `CPsiBlastImpl::Run`, psiblast_impl.cpp:150), so the reports compare its letters as
//! they are (`QueryCase::Bioseq`). A subject's Bioseq always holds upper case (a database's
//! NCBIstdaa letters, or a FASTA record's), and its rows' lower case is a mask.

use crate::common::GapEditOp;
use crate::config::ScoringMatrix;

/// How the reports compare the letters of a query row (`QueryLetters::case`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) enum QueryCase {
    /// A record's row: a lower-case letter is a mask on the capital of the Bioseq, which NCBI
    /// compares.
    #[default]
    Masked,
    /// The letters of a PSSM query's Bioseq as they are (`PsiSearch::bioseq_letters`), without
    /// masks: NCBI compares a lower-case letter as it is, so it is identical to no subject letter
    /// and scores -1000 in the display matrix (`protein_display_score`).
    Bioseq,
}

impl QueryCase {
    /// The letter of the query's Bioseq that NCBI compares for each row letter, by byte: a loop
    /// over a row's letters takes the table once and looks each letter up, with no test of the
    /// case per letter.
    #[inline(always)]
    pub(crate) fn compared_letters(self) -> &'static [u8; 256] {
        match self {
            Self::Masked => &MASKED_COMPARED,
            Self::Bioseq => &BIOSEQ_COMPARED,
        }
    }
}

/// The table of [`QueryCase::compared_letters`]: each byte, in upper case where `upper`.
const fn compared_table(upper: bool) -> [u8; 256] {
    let mut table = [0u8; 256];
    let mut byte = 0;
    while byte < 256 {
        table[byte] = if upper {
            (byte as u8).to_ascii_uppercase()
        } else {
            byte as u8
        };
        byte += 1;
    }
    table
}

/// `QueryCase::Masked`'s letters: the capital of each lower-case letter.
static MASKED_COMPARED: [u8; 256] = compared_table(true);
/// `QueryCase::Bioseq`'s letters: each as it is.
static BIOSEQ_COMPARED: [u8; 256] = compared_table(false);

/// A query's letters as the reports read them: its Bioseq's (`CSeqVector` in Iupac coding), which
/// are a FASTA record's letters in upper case, a PSSM query's record letters (the capitals, `*`,
/// and `-` for the gap residue), or that PSSM query's Bioseq letters where they are not its
/// record's (`bioseq`: lower-case iupacaa, `PsiSearch::bioseq_letters`).
#[derive(Clone, Copy, Debug)]
pub(crate) struct QueryLetters<'a> {
    record: &'a [u8],
    bioseq: Option<&'a [u8]>,
}

impl<'a> QueryLetters<'a> {
    /// The letters of a query whose record holds `record` and whose Bioseq, where its letters are
    /// not the record's, `bioseq` (of the same length).
    pub(crate) fn new(record: &'a [u8], bioseq: Option<&'a [u8]>) -> Self {
        Self { record, bioseq }
    }

    /// The Bioseq's letter at `position` (0-based), `None` past the end.
    #[inline]
    pub(crate) fn get(self, position: usize) -> Option<u8> {
        match self.bioseq {
            Some(bioseq) => bioseq.get(position).copied(),
            None => self.record.get(position).map(u8::to_ascii_uppercase),
        }
    }

    /// Whether the reports count an alignment of these letters on other letters than the search
    /// read: the Bioseq's (`bioseq`), or a record's gap residue (`-` of a PSSM query, which the
    /// search scores by its column and the rows show as the letter of a gap, identical to a gap
    /// of the other row). A FASTA record holds no `-` (NCBI's reader drops it).
    pub(crate) fn shows_other_letters(self) -> bool {
        self.bioseq.is_some() || self.record.contains(&b'-')
    }

    /// How the reports compare the letters of a row made of these letters.
    pub(crate) fn case(self) -> QueryCase {
        if self.bioseq.is_some() {
            QueryCase::Bioseq
        } else {
            QueryCase::Masked
        }
    }
}

/// The sequence of a protein alignment that a row shows.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum RowSide {
    Query,
    Subject,
}

impl RowSide {
    /// Whether a column of the operation `op` holds a residue of this side's sequence: a
    /// substitution's does, a deletion is a gap in the query, an insertion a gap in the subject.
    #[inline]
    pub(crate) fn reads(self, op: &GapEditOp) -> bool {
        matches!(
            (self, op),
            (_, GapEditOp::Sub(_))
                | (Self::Query, GapEditOp::Ins(_))
                | (Self::Subject, GapEditOp::Del(_))
        )
    }
}

/// The columns (0-based) of the rows of an alignment, whose operations are `ops`, where the
/// letter of `side` is a gap character: a residue `gap` of `sequence`, whose residues the
/// alignment reads from `start` (0-based) on (`pairwise::PairwiseHit::query_stored_gaps`,
/// `subject_stored_gaps`; `pairwise::end_gap_span`). Empty without such a residue.
pub(crate) fn stored_gap_columns(
    sequence: &[u8],
    gap: u8,
    start: usize,
    ops: &[GapEditOp],
    side: RowSide,
) -> Vec<u32> {
    let residues: usize = ops
        .iter()
        .filter(|op| side.reads(op))
        .map(|op| op.num() as usize)
        .sum();
    let end = (start + residues).min(sequence.len());
    if !sequence
        .get(start..end)
        .is_some_and(|range| range.contains(&gap))
    {
        return Vec::new();
    }
    let mut columns = Vec::new();
    let (mut column, mut index) = (0u32, start);
    for op in ops {
        let reads = side.reads(op);
        for _ in 0..op.num() {
            if reads {
                if sequence.get(index) == Some(&gap) {
                    columns.push(column);
                }
                index += 1;
            }
            column += 1;
        }
    }
    columns
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2120-2154
// ```c++
// void CDisplaySeqalign::x_FillIdentityInfo(const string& sequence_standard,
//                                           const string& sequence ,
//                                           int& match, int& positive,
//                                           string& middle_line)
// {
//     match = 0;
//     positive = 0;
//     int min_length=min<int>((int)sequence_standard.size(), (int)sequence.size());
//     if(m_AlignOption & eShowMiddleLine){
//         middle_line = sequence;
//     }
//     for(int i=0; i<min_length; i++){
//         if(sequence_standard[i]==sequence[i]){
//             if(m_AlignOption & eShowMiddleLine){
//                 if(m_MidLineStyle == eBar ) {
//                     middle_line[i] = '|';
//                 } else if (m_MidLineStyle == eChar){
//                     middle_line[i] = sequence[i];
//                 }
//             }
//             match ++;
//         } else {
//             if ((m_AlignType&eProt)
//                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
//                 positive ++;
//                 if(m_AlignOption & eShowMiddleLine){
//                     if (m_MidLineStyle == eChar){
//                         middle_line[i] = '+';
//                     }
//                 }
//             } else {
//                 if (m_AlignOption & eShowMiddleLine){
//                     middle_line[i] = ' ';
//                 }
//             }
// ```
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1794-1798
// ```c++
//     x_FillIdentityInfo(aln_vec_info->alnRowInfo->sequence[0],
//                        aln_vec_info->alnRowInfo->sequence[1],
//                        aln_vec_info->match,
//                        aln_vec_info->positive,
//                        aln_vec_info->alnRowInfo->middleLine);
// ```
/// The positives of the rows `query` and `subject` (identities included), as NCBI counts them
/// on the Bioseqs' letters before a mask lowers any (`QueryCase`: the query's letters as `case`
/// reads them, the subject's capitals).
pub(crate) fn display_positives(
    query: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    case: QueryCase,
) -> usize {
    let compared = case.compared_letters();
    query
        .iter()
        .zip(subject)
        .filter(|&(&q, &s)| {
            let q = compared[usize::from(q)];
            let s = s.to_ascii_uppercase();
            q == s || crate::utils::matrix::protein_display_score(matrix, q, s) > 0
        })
        .count()
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2131-2155 and c++/src/objtools/align_format/align_format_util.cpp:3563-3581 (quoted above and at
    // `protein_display_score`): a masked capital compares as its capital; a PSSM query's
    // lower-case letter is identical to no subject letter and scores -1000.
    #[test]
    fn a_masked_letter_compares_as_its_capital_and_a_bioseq_letter_as_it_is() {
        let matrix = ScoringMatrix::Blosum62;
        assert_eq!(
            display_positives(b"aKv", b"AKI", matrix, QueryCase::Masked),
            3
        );
        assert_eq!(
            display_positives(b"aKv", b"AKI", matrix, QueryCase::Bioseq),
            1
        );
        assert_eq!(
            display_positives(b"-A", b"-a", matrix, QueryCase::Bioseq),
            2
        );
        let pssm = QueryLetters::new(b"AKV", Some(b"aKv"));
        assert_eq!(
            (pssm.get(0), pssm.get(1), pssm.get(3)),
            (Some(b'a'), Some(b'K'), None)
        );
        assert_eq!(pssm.case(), QueryCase::Bioseq);
        let record = QueryLetters::new(b"aKv", None);
        assert_eq!((record.get(0), record.get(2)), (Some(b'A'), Some(b'V')));
        assert_eq!(record.case(), QueryCase::Masked);
        assert_eq!(
            QueryCase::Masked.compared_letters()[usize::from(b'k')],
            b'K'
        );
        assert_eq!(
            QueryCase::Bioseq.compared_letters()[usize::from(b'k')],
            b'k'
        );
        assert_eq!(
            QueryCase::Masked.compared_letters()[usize::from(b'-')],
            b'-'
        );
        assert!(pssm.shows_other_letters());
        assert!(QueryLetters::new(b"A-K", None).shows_other_letters());
        assert!(!record.shows_other_letters());
    }

    // A gap character is a residue of its own side's sequence: the query reads one at a
    // substitution and an insertion, the subject at a substitution and a deletion.
    #[test]
    fn the_gap_residues_lie_at_their_sides_columns() {
        let ops = [
            GapEditOp::Sub(2),
            GapEditOp::Ins(1),
            GapEditOp::Del(1),
            GapEditOp::Sub(1),
        ];
        // Query `XA--K` read from 1: its `-` at positions 2 and 3 lie at columns 1 and 2 (the
        // insertion).
        assert_eq!(
            stored_gap_columns(b"XA--K", b'-', 1, &ops, RowSide::Query),
            [1, 2]
        );
        // Subject residues 0, 9, 0, 5 at columns 0, 1, 3, 4: code 0 at 0 and 3.
        assert_eq!(
            stored_gap_columns(&[0, 9, 0, 5], 0, 0, &ops, RowSide::Subject),
            [0, 3]
        );
        assert!(stored_gap_columns(b"MKVLA", b'-', 0, &ops, RowSide::Query).is_empty());
    }
}
