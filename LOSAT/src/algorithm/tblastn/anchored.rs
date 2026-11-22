//! TBLASTN's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them
//! (`stage_e_report::pairwise_hits`).
//!
//! - Letters: the searched protein letters of the queries, upper case (a PSSM query's Bioseq
//!   letters as they are, `row_letters::QueryLetters`), and the nucleotides of the
//!   subjects (the records' letters, or a packed subject's), in their records' positions (the
//!   searched letters start at the `-query_loc` or `-subject_loc` interval); the section translates
//!   a subject's rows with the database's genetic code (`subject_genetic_code`, NCBI's
//!   `m_SlaveGeneticCode`, `-db_gencode`).
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::record_label`,
//!   `query_anchored::database_label`).
//! - Masks: the query's SEG and `-lcase_masking` masks (`QueryFrame::seg_masks`, as outfmt 0's
//!   rows show them; none for a PSSM query). The subjects' masks have no frame
//!   (`CSeqLocInfo::eFrameNotSet`), and a subject row of width 3 has one, so NCBI shows none of them
//!   (`query_anchored::prepare_rows`): none are read.

use super::stage_e_report::ShownSubjects;
use super::x_subject_store::SubjectLetters;
use crate::algorithm::tblastx::report::DisplayCode;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::query_anchored::{
    database_label, query_record_letters, record_label, Mask, ReportSequences,
};
use crate::report::row_letters::QueryLetters;

/// What the sections of a TBLASTN report read (module doc).
pub(super) struct AnchoredSequences<'a> {
    pub(super) queries: &'a [FastaRecord],
    /// The PSSM query's Bioseq letters where they are not its record's
    /// (`PsiSearch::bioseq_letters`), which the sections show as they are (`QueryLetters`); a run
    /// of such a query has no other.
    pub(super) query_bioseq: Option<&'a [u8]>,
    pub(super) query_placements: &'a Placements,
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) subject_letters: SubjectLetters<'a>,
    pub(super) subject_placements: &'a Placements,
    pub(super) shown: &'a ShownSubjects<'a>,
    pub(super) show_gis: bool,
    /// The masks of each query, inclusive positions in its record.
    pub(super) query_masks: Vec<Vec<(u32, u32)>>,
    /// The database's genetic code.
    pub(super) code: DisplayCode,
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A protein query is never a nucleotide subject's Bioseq.
    fn subject_is_query(&self, _q_idx: usize, _s_idx: u32) -> bool {
        false
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        query_record_letters(
            QueryLetters::new(self.queries[q_idx].seq(), self.query_bioseq),
            self.query_placements.offset(q_idx),
            from,
            to,
            out,
        );
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        let shift = self.subject_placements.offset(s_idx);
        let subject = self.subject_letters.nt(s_idx);
        let end = (to as usize + 1 - shift).min(subject.len());
        let start = (from as usize - shift).min(end);
        let mut buffer = Vec::new();
        out.extend(
            subject
                .fetch(start, end, &mut buffer)
                .iter()
                .map(u8::to_ascii_uppercase),
        );
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        record_label(&self.queries[q_idx], self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        let rank = s_idx as usize;
        match self.shown.oid(rank) {
            Some(oid) => database_label(self.shown.parsed_subject(rank), oid, out),
            None => record_label(&self.subject_records[rank], self.show_gis, out),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        if let Some(masks) = self.query_masks.get(q_idx) {
            out.extend(masks.iter().map(|&(from, to)| (from, to, 0)));
        }
    }

    fn subject_masks(&self, _q_idx: usize, _s_idx: u32, _out: &mut Vec<Mask>) {}

    fn subject_genetic_code(&self) -> Option<DisplayCode> {
        Some(self.code)
    }
}
