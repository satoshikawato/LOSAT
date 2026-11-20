//! BLASTP's sequences for the query-anchored reports (outfmt 1 to 4, `report::query_anchored`):
//! the letters, labels and masks of a query's section, read where outfmt 0's rows read them.
//!
//! - Letters: the searched letters as the rows decode them (`ncbistdaa_to_ascii` of the query's
//!   unmasked context, `query_nomask_sequence`, and of the subject's residues), with an O that the
//!   search read as X shown as O (`displayed_letter`) and a PSSM query's Bioseq letters as they are
//!   (lower-case iupacaa: `ShownQueries`, `ShownLetters::shown`); a section's positions are the
//!   records', the searched letters start at the `-query_loc` or `-subject_loc` interval.
//! - Labels: `x_FillSeqid` of the Bioseq's Seq-ids (`query_anchored::fill_seqid`): a record read
//!   with `-parse_deflines` by its ids, else its local id (`Query_N`, `Subject_N`); a database
//!   sequence by its parsed Seq-ids (`ParsedSubject::anchored_label`), else by its OID.
//! - Masks: the query's `TMaskedQueryRegions` (`query_record_masks`) and each subject's
//!   `-lcase_masking` runs of its HSP list with the query (`x_soft_subject::x_subject_masks`).

use std::collections::HashMap;
use std::sync::Arc;

use super::{displayed_letter, query_nomask_sequence, ShownLetters, ShownQueries, ShownSubjects};
use crate::algorithm::blastp::encoding::{ncbistdaa_to_ascii, EncodedProtein};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::db_subject::ParsedSubject;
use crate::report::query_anchored::{
    database_label, record_label, subject_is_query_record, Mask, ReportSequences,
};

/// What the sections of a BLASTP report read (module doc).
pub(super) struct AnchoredSequences<'a> {
    pub(super) contexts: &'a [QueryContext],
    /// The query records (their letters, with a PSSM query's Bioseq letters, and their labels).
    pub(super) queries: ShownQueries<'a>,
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

/// Appends the masks `masks` as positions.
fn push_masks(masks: &[(i32, i32)], out: &mut Vec<Mask>) {
    out.extend(
        masks
            .iter()
            .map(|&(from, to)| (from.max(0) as u32, to.max(0) as u32, 0)),
    );
}

impl ReportSequences for AnchoredSequences<'_> {
    /// A `-subject` record read with `-parse_deflines` that is the query's record
    /// (`subject_is_query_record`; a database sequence with the query's id is refused before the
    /// reports).
    fn subject_is_query(&self, q_idx: usize, s_idx: u32) -> bool {
        self.listed_ordinals.is_none()
            && match (
                self.queries.records.get(q_idx),
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
            .queries
            .letters(q_idx)
            .unwrap_or(ShownLetters::Ascii(&[]));
        for position in from as usize - shift..=to as usize - shift {
            let Some(&residue) = letters.get(position) else {
                break;
            };
            out.push(input.shown(ncbistdaa_to_ascii(residue), position) as u8);
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
        record_label(&self.queries.records[q_idx], self.show_gis, out);
    }

    fn subject_label(&self, s_idx: u32, out: &mut Vec<u8>) {
        let s_idx = s_idx as usize;
        match self.listed_ordinals {
            // A database sequence: its parsed Seq-ids, else its OID (`database_label`).
            Some((ordinals, oids)) => database_label(
                self.parsed.get(s_idx).and_then(Option::as_deref),
                oids.oid(ordinals[s_idx]),
                out,
            ),
            None => record_label(&self.subject_records[s_idx], self.show_gis, out),
        }
    }

    fn query_masks(&self, q_idx: usize, out: &mut Vec<Mask>) {
        if let Some(masks) = self.query_masks.and_then(|masks| masks.get(q_idx)) {
            push_masks(masks, out);
        }
    }

    fn subject_masks(&self, q_idx: usize, s_idx: u32, out: &mut Vec<Mask>) {
        if let Some(masks) = self.subject_masks.get(&(q_idx as u32, s_idx)) {
            push_masks(masks, out);
        }
    }
}
