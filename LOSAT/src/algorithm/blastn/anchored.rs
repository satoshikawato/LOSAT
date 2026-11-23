//! BLASTN's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them
//! (`pairwise::Rows`).
//!
//! - Letters: the letters of the queries' records (with `-query_loc` the whole record, as outfmt
//!   0's rows read them, `pairwise::QueryRecords`) and the subjects' searched letters, upper case,
//!   in their records' positions (the searched letters start at the `-subject_loc` interval); a
//!   database's reported sequences are the letters that outfmt 0 reads (`SubjectLetters`).
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::record_label`,
//!   `query_anchored::database_label`).
//! - Masks: the query's `TMaskedQueryRegions` (`pairwise::Rows::query_masks`) and each subject's
//!   masks of its HSP list with the query (`GetSubjectMasks`, `pairwise::Rows::subject_masks`).

use std::sync::Arc;

use super::hsp::BlastnHitList;
use super::pairwise::{QueryRecords, Rows};
use super::x_subject_store::SubjectLetters;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::db_subject::ParsedSubject;
use crate::report::query_anchored::{
    database_label, query_record_letters, record_label, subject_is_query_record, Mask,
    ReportSequences,
};

/// What the sections of a BLASTN report read (module doc).
pub(crate) struct AnchoredSequences<'a> {
    pub(crate) queries: &'a [FastaRecord],
    /// The records whose letters the sections read (`pairwise::QueryRecords`).
    pub(crate) query_records: QueryRecords<'a>,
    /// The subjects' letters by subject index (a database's by rank).
    pub(crate) subjects: SubjectLetters<'a>,
    pub(crate) subject_placements: &'a Placements,
    /// The database sequences with parsed Seq-ids by rank (a `-db` search).
    pub(crate) parsed: &'a [Option<Arc<ParsedSubject>>],
    /// The ordinal id of each rank and the database's OID of each ordinal id (a `-db` search;
    /// `None` for `-subject` records).
    pub(crate) listed_ordinals: Option<(&'a [u32], crate::seqdb::OidMap<'a>)>,
    pub(crate) show_gis: bool,
    /// The masks that the report shows (outfmt 0's `Rows`).
    pub(crate) rows: &'a Rows<'a>,
    /// The final hit lists of the report's queries (the subjects' masks of each list).
    pub(crate) hit_lists: &'a [Option<BlastnHitList>],
}

/// Appends the half-open masks `masks` as inclusive positions (`shown_masks`).
fn push_masks(masks: &[crate::utils::dust::MaskedInterval], out: &mut Vec<Mask>) {
    out.extend(shown_masks(masks));
}

/// The half-open masks `masks` as inclusive positions (`Mask`), those of no position left out:
/// the masks of the query-anchored reports and of outfmt 0's rows (`pairwise::Rows`).
pub(crate) fn shown_masks(
    masks: &[crate::utils::dust::MaskedInterval],
) -> impl Iterator<Item = Mask> + '_ {
    masks
        .iter()
        .filter(|mask| mask.start < mask.end)
        .map(|mask| (mask.start as u32, mask.end as u32 - 1, 0))
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A `-subject` record read with `-parse_deflines` that is the query's record
    /// (`subject_is_query_record`; a database sequence with the query's id is refused before the
    /// reports).
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool {
        self.listed_ordinals.is_none()
            && match (
                self.queries.get(q_idx),
                self.subjects.records().get(s_idx as usize),
            ) {
                (Some(query), Some(subject)) => subject_is_query_record(query, subject),
                _ => false,
            }
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        // NCBI reference (598d8ae6): c++/include/objtools/alnmgr/alnvec.hpp:304-305
        // The letters of the query's Bioseq at record positions (`QueryRecord::letters`, from
        // position 0), whose `CSeqVector` stops at the record's end.
        if let Some(record) = self.query_records.record(q_idx) {
            query_record_letters(
                crate::report::row_letters::QueryLetters::new(record.letters, None),
                0,
                from,
                to,
                out,
            );
        }
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        let shift = self.subject_placements.offset(s_idx);
        let end = (to as usize + 1 - shift).min(self.subjects.len(s_idx));
        let start = (from as usize - shift).min(end);
        out.extend(self.subjects.upper_letters(s_idx, start, end));
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        record_label(&self.queries[q_idx], self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        match self.listed_ordinals {
            Some((ordinals, oids)) => database_label(
                self.parsed.get(s_idx).and_then(Option::as_deref),
                oids.oid(ordinals[s_idx]),
                out,
            ),
            None => record_label(&self.subjects.records()[s_idx], self.show_gis, out),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        push_masks(self.rows.query_masks(q_idx), out);
    }

    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<Mask>) {
        let Some(hit_list) = self.hit_lists.get(q_idx).and_then(Option::as_ref) else {
            return;
        };
        if let Some(hsp_list) = hit_list
            .hsplist_array
            .iter()
            .find(|list| list.hsps.first().is_some_and(|hsp| hsp.s_idx == s_idx))
        {
            push_masks(&self.rows.subject_masks(&hsp_list.hsps), out);
        }
    }
}
