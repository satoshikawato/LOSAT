//! BLASTP's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them.
//!
//! - Letters: the searched letters as the rows decode them (`ncbistdaa_to_ascii` of the query's
//!   unmasked context, `query_nomask_sequence`, and of the subject's residues), with an O that the
//!   search read as X shown as O (`displayed_letter`); a section's positions are the records', the
//!   searched letters start at the `-query_loc` or `-subject_loc` interval.
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::fill_seqid`): a record read
//!   with `-parse_deflines` by its ids, else its local id (`Query_N`, `Subject_N`); a database
//!   sequence by its parsed Seq-ids (`ParsedSubject::anchored_label`), else by its OID.
//! - Masks: the query's `TMaskedQueryRegions` (`query_record_masks`) and each subject's
//!   `-lcase_masking` runs of its HSP list with the query (`x_soft_subject::x_subject_masks`).

use std::collections::HashMap;
use std::sync::Arc;

use super::{displayed_letter, query_nomask_sequence, ShownSubjects};
use crate::algorithm::blastp::encoding::{ncbistdaa_to_ascii, EncodedProtein};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::db_subject::ParsedSubject;
use crate::report::query_anchored::{fill_seqid, ReportSequences};

/// What the sections of a BLASTP report read (module doc).
pub(super) struct AnchoredSequences<'a> {
    pub(super) contexts: &'a [QueryContext],
    pub(super) query_records: &'a [FastaRecord],
    pub(super) query_placements: &'a Placements,
    /// The encoded subjects by subject index (a database's by rank), with their sentinels.
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) shown_subjects: ShownSubjects<'a>,
    /// The records of a `-subject` set (their labels).
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) subject_placements: &'a Placements,
    /// The database sequences with parsed Seq-ids by rank, and the ordinal id of each rank with
    /// the database's OID of each ordinal id (a `-db` search; `None` for `-subject` records).
    pub(super) parsed: &'a [Option<Arc<ParsedSubject>>],
    pub(super) listed_ordinals: Option<(&'a [u32], crate::seqdb::OidMap<'a>)>,
    pub(super) show_gis: bool,
    pub(super) query_masks: Option<&'a [Vec<(i32, i32)>]>,
    pub(super) subject_masks: HashMap<(u32, u32), Vec<(i32, i32)>>,
}

/// Whether the `-subject` record `subject` is query `query`'s own record: the query is one record
/// with a subject (the same definition line and letters, `BelievedIds::shared`), and the subject
/// has its id (any other subject with the query's id is refused before the reports).
pub(super) fn subject_is_query_record(query: &FastaRecord, subject: &FastaRecord) -> bool {
    match (query.believed(), subject.believed()) {
        (Some(query), Some(subject)) => query.shared.is_some() && subject.best == query.best,
        _ => false,
    }
}

/// The label of a record's rows: its ids with `-parse_deflines`, else its local id.
fn record_label(record: &FastaRecord, show_gis: bool, out: &mut Vec<u8>) {
    match record.believed() {
        Some(ids) => fill_seqid(&ids.ids, show_gis, out),
        None => record.seq_id().sam_name(out),
    }
}

/// Appends the masks `masks` as positions.
fn push_masks(masks: &[(i32, i32)], out: &mut Vec<(u32, u32)>) {
    out.extend(
        masks
            .iter()
            .map(|&(from, to)| (from.max(0) as u32, to.max(0) as u32)),
    );
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A `-subject` record read with `-parse_deflines` that is the query's record
    /// (`subject_is_query_record`; a database sequence with the query's id is refused before the
    /// reports).
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool {
        self.listed_ordinals.is_none()
            && match (
                self.query_records.get(q_idx),
                self.subject_records.get(s_idx as usize),
            ) {
                (Some(query), Some(subject)) => subject_is_query_record(query, subject),
                _ => false,
            }
    }

    fn query_letters(&self, q_idx: usize, from: u32, to: u32, out: &mut Vec<u8>) {
        let shift = self.query_placements.offset(q_idx);
        let letters = query_nomask_sequence(&self.contexts[q_idx]);
        let input = self
            .query_records
            .get(q_idx)
            .map_or(&[][..], |record| record.seq());
        for position in from as usize - shift..=to as usize - shift {
            let Some(&residue) = letters.get(position) else {
                break;
            };
            let letter =
                displayed_letter(ncbistdaa_to_ascii(residue), input.get(position).copied());
            out.push(letter as u8);
        }
    }

    fn subject_letters(&self, s_idx: u32, from: u32, to: u32, out: &mut Vec<u8>) {
        let shift = self.subject_placements.offset(s_idx as usize);
        let encoded = &self.subjects[s_idx as usize].aa_seq;
        let letters = encoded
            .get(1..encoded.len().saturating_sub(1))
            .unwrap_or_default();
        let input = self.shown_subjects.letters(s_idx as usize);
        for position in from as usize - shift..=to as usize - shift {
            let Some(&residue) = letters.get(position) else {
                break;
            };
            let shown = input.and_then(|input| input.get(position));
            out.push(displayed_letter(ncbistdaa_to_ascii(residue), shown) as u8);
        }
    }

    fn query_label(&self, q_idx: usize, out: &mut Vec<u8>) {
        record_label(&self.query_records[q_idx], self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        match self.listed_ordinals {
            // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbvol.cpp:2231-2236
            // ```c++
            //                 if (seqid.Which() == CSeq_id::e_General) {
            //                     CDbtag & dbt = seqid.SetGeneral();
            //
            //                     if (dbt.GetDb() == "BL_ORD_ID") {
            //                         int vol_oid = dbt.GetTag().GetId();
            //                         dbt.SetTag().SetId(m_VolStart + vol_oid);
            // ```
            // A database sequence: its parsed Seq-ids, else `gnl|BL_ORD_ID|<oid>`, whose label is
            // the tag (`GetLabel` of a `general` id, `query_anchored::fill_seqid`): its OID in the
            // opened database, not its ordinal id where a list or an alias range leaves sequences
            // out (`seqdb::Oid`, as `report::record::SubjectNames::Database`; review 3 of T3a, D2).
            Some((ordinals, oids)) => match self.parsed.get(s_idx).and_then(Option::as_deref) {
                Some(parsed) => parsed.anchored_label(out),
                None => {
                    out.extend_from_slice(oids.oid(ordinals[s_idx]).get().to_string().as_bytes())
                }
            },
            None => record_label(&self.subject_records[s_idx], self.show_gis, out),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<(u32, u32)>) {
        if let Some(masks) = self.query_masks.and_then(|masks| masks.get(q_idx)) {
            push_masks(masks, out);
        }
    }

    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<(u32, u32)>) {
        if let Some(masks) = self.subject_masks.get(&(q_idx as u32, s_idx)) {
            push_masks(masks, out);
        }
    }
}

/// The refusal of a run whose query `record` (one record as a query and as a subject under
/// `-parse_deflines`) has a query-anchored section that NCBI BLAST+ aborts on
/// (`pairwise::blastp_section_aborts`), named as the refusal of its masks in outfmt 0 names it
/// (`believed::refuse_shared_masks`).
pub(super) fn aborted_section(record: &FastaRecord, program: &str) -> anyhow::Error {
    let mut fasta = Vec::new();
    if let Some(ids) = record.believed() {
        crate::seqid::label::fasta(&ids.best, &mut fasta);
    }
    anyhow::anyhow!(
        "the record {} is both a query and a subject under -parse_deflines and NCBI BLAST+ aborts \
         on its query-anchored alignments (outfmt 1 to 4), whose merge leaves out the subject's \
         rows; this is not supported by LOSAT's {program}",
        String::from_utf8_lossy(&fasta)
    )
}
