//! BLASTX's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them
//! (`report::alignment`, `report::pairwise_hits`).
//!
//! - Letters: the searched nucleotides of the queries and the residues of the subjects (the
//!   records' letters, a database's reported sequences), upper case, in their records' positions
//!   (the searched letters start at the `-query_loc` or `-subject_loc` interval); the section
//!   translates a query row with the query's genetic code (`query_genetic_code`, NCBI's
//!   `m_MasterGeneticCode`, `-query_gencode`).
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::seq_id_label`,
//!   `query_anchored::database_label`).
//! - Masks: the query's masks of every frame (`blast_xml::query_masks`, NCBI's
//!   `TMaskedQueryRegions` with each context's frame, which a group's query row shows where it is
//!   that group's frame), and each subject's `-lcase_masking` runs of its HSP list with the query
//!   (`x_soft_subject::x_subject_masks`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1547-1548,1579-1584
//! (`CBlastFormat::PrintOneResultSet`, formats 0 to 4)
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2100-2101
//! NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1244-1251
//! (`x_FillSeqid`, the labels: showalign.cpp:4545-4572, cited at `query_anchored::fill_seqid`.)

use std::collections::HashMap;

use std::sync::Arc;

use super::input::FastaRecord;
use super::runtime::BatchResults;
use crate::algorithm::tblastx::report::DisplayCode;
use crate::report::db_subject::ParsedSubject;
use crate::report::query_anchored::{database_label, seq_id_label, Mask, ReportSequences};

/// What the sections of a BLASTX report read (module doc).
pub(super) struct AnchoredSequences<'a> {
    pub(super) queries: &'a [FastaRecord],
    pub(super) subjects: &'a [FastaRecord],
    /// A database's reported sequences by rank: their parsed Seq-ids and OIDs (`None` for
    /// `-subject` records).
    pub(super) database: Option<(&'a [Option<Arc<ParsedSubject>>], &'a [crate::seqdb::Oid])>,
    pub(super) show_gis: bool,
    /// The masks of each query (`Mask`), by query index.
    pub(super) query_masks: Vec<Vec<Mask>>,
    /// The subjects' masks by (query, subject) (`x_soft_subject::x_subject_masks`).
    pub(super) subject_masks: HashMap<(u32, u32), Vec<(i32, i32)>>,
    /// The query's genetic code.
    pub(super) code: DisplayCode,
}

impl<'a> AnchoredSequences<'a> {
    /// The masks of the queries of `batches` (by query index of `queries`): the masks of each
    /// query's contexts in their order (`blast_xml::query_masks`).
    pub(super) fn query_masks(queries: &[FastaRecord], batches: &[BatchResults]) -> Vec<Vec<Mask>> {
        let mut masks = vec![Vec::new(); queries.len()];
        for batch in batches {
            for qidx in 0..batch.queries.len() {
                let absolute = batch.query_ordinal + qidx;
                masks[absolute] = query_record_masks(batch, qidx, &queries[absolute]);
            }
        }
        masks
    }
}

/// The masks of query `qidx` of `batch` (its record `query`) as `Mask`s: its contexts' masks in
/// their order, each with its context's frame, in record positions (`blast_xml::query_masks`,
/// NCBI's `TMaskedQueryRegions`). The query-anchored reports (`AnchoredSequences::query_masks`)
/// and outfmt 0's query rows (`report::pairwise_hits`) place them (`x_FillLocList`).
pub(super) fn query_record_masks(
    batch: &BatchResults,
    qidx: usize,
    query: &FastaRecord,
) -> Vec<Mask> {
    super::blast_xml::query_masks(batch, qidx, query)
        .into_iter()
        .map(|mask| (mask.from as u32, mask.to as u32, mask.frame as i8))
        .collect()
}

/// Appends the letters `from..=to` (positions in `record`'s whole record) of `record`, whose
/// searched letters start at its interval's start, upper case, as far as the record has them.
fn record_letters(record: &FastaRecord, from: u32, to: u32, out: &mut Vec<u8>) {
    let shift = record.placement().offset;
    for position in from as usize - shift..=to as usize - shift {
        let Some(&letter) = record.sequence.get(position) else {
            break;
        };
        out.push(letter.to_ascii_uppercase());
    }
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A nucleotide query is never a protein subject's Bioseq.
    fn subject_is_query(&self, _q_idx: usize, _s_idx: u32) -> bool {
        false
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        record_letters(&self.queries[q_idx], from, to, out);
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        record_letters(&self.subjects[s_idx as usize], from, to, out);
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        let query = &self.queries[q_idx];
        seq_id_label(query.believed(), query.seq_id(), self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        match self.database {
            Some((parsed, oids)) => database_label(
                parsed.get(s_idx as usize).and_then(Option::as_deref),
                oids[s_idx as usize],
                out,
            ),
            None => {
                let subject = &self.subjects[s_idx as usize];
                seq_id_label(subject.believed(), subject.seq_id(), self.show_gis, out);
            }
        }
    }

    fn query_seq_id(&self, q_idx: usize) -> crate::report::seqid::SeqId<'_> {
        self.queries[q_idx].seq_id()
    }

    fn subject_seq_id(&self, s_idx: u32) -> crate::report::seqid::SeqId<'_> {
        match self.database {
            Some((parsed, oids)) => crate::report::query_anchored::database_seq_id(
                parsed.get(s_idx as usize).and_then(Option::as_deref),
                oids[s_idx as usize],
            ),
            None => self.subjects[s_idx as usize].seq_id(),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        if let Some(masks) = self.query_masks.get(q_idx) {
            out.extend_from_slice(masks);
        }
    }

    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<Mask>) {
        if let Some(masks) = self.subject_masks.get(&(q_idx as u32, s_idx)) {
            out.extend(
                masks
                    .iter()
                    .map(|&(from, to)| (from.max(0) as u32, to.max(0) as u32, 0)),
            );
        }
    }

    fn query_genetic_code(&self) -> Option<DisplayCode> {
        Some(self.code)
    }
}
