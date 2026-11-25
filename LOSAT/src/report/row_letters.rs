//! The letters of an alignment's rows as NCBI's reports read and compare them: the query's
//! letters (`QueryLetters`: a FASTA record's, or the Bioseq's of a PSSM query), how a query row's
//! lower case reads (`QueryCase`), the positives of two rows (`display_positives`), and the
//! columns whose `-` is a residue (`stored_gap_columns`).
//!
//! NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-121
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
//! NCBI's reports read the letters of every row from its Bioseq and compare them as they are; a
//! mask lowers a letter only where outfmt 0 prints the row. A FASTA record's Bioseq holds its
//! letters in upper case (NCBI's reader keeps the lower case as a mask, fasta.cpp:945-949), and
//! the rows of LOSAT's hit records may carry those masks in lower case (`HspMasks::lower_row`,
//! as the library API hands them out), so the reports compare the capitals of such a row
//! (`QueryCase::Masked`). The Bioseq of a PSSM query (`-in_pssm`) holds iupacaa letters
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
/// `subject_stored_gaps`). Empty without such a residue. A row's end gaps are read from the edit
/// script (`serial::seqalign::row_end_gaps`), where such a residue is a letter.
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
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1794-1798
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
#[path = "row_letters_tests.rs"]
mod tests;
