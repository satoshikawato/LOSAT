//! BLASTN's SAM records (outfmt 17): the report records of each query's final hit list
//! (`report::record`) added to the run's [`SamDocument`], one query at a time. The records are
//! built as for the Seq-annot of outfmt 8 and 12 (`algorithm::blastn::seqalign`); the SAM
//! formatter reads the Dense-seg segments, the scores, the sequences' lengths and, with `SQ`,
//! their letters.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:750-755
//! ```c++
//!     if (m_FormatType == CFormattingArgs::eSAM) {
//!         if (results.HasAlignments()) {
//!             m_SamFormatter->Print(*results.GetSeqAlign());
//!         }
//!         return;
//!     }
//! ```
//! (a query without alignments writes nothing: no unmapped record). The Seq-aligns are those of
//! `BlastHitList2SeqAlign_OMF`, without `PruneSeqalign` (`algorithm::blastn::seqalign`).

use anyhow::{Context, Result};

use super::hsp::BlastnHitList;
use super::x_subject_store::SubjectLetters;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::sam::{SamDocument, SamId, SamLetters, SamRow, SamSequence};
use crate::report::seqid::SeqId;

/// The lengths of the subjects the hit lists name, by index (rank for a database).
pub(crate) struct SubjectLengths<'a> {
    /// The letters that were searched: `SubjectRecord::searched_length`.
    pub searched: &'a dyn Fn(usize) -> usize,
    /// The whole Bioseq (`GetBioseqLength`): also where `-subject_loc` searched an interval.
    pub whole: &'a dyn Fn(usize) -> usize,
}

/// Adds the records of the queries of `hit_lists` (indexed by query) to `document`.
/// `queries_before` is the number of queries earlier calls added (a query's reference key).
#[allow(clippy::too_many_arguments)]
pub(crate) fn add_blastn_sam(
    document: &mut SamDocument,
    hit_lists: &[Option<BlastnHitList>],
    queries: &QueryNames<'_>,
    queries_before: usize,
    subjects: SubjectNames<'_>,
    subject_placements: &Placements,
    lengths: &SubjectLengths<'_>,
    letters: &mut dyn SamLetters,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastn,
        db_title: subjects.db_title(),
        gapped: true,
    };
    let mut sam_subjects: Vec<SamSequence<'_>> = Vec::new();
    let mut hsp_records: Vec<HspRecord<'_>> = Vec::new();
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        sam_subjects.clear();
        hsp_records.clear();
        for hsp_list in &hit_list.hsplist_array {
            if hsp_list.hsps.is_empty() {
                continue;
            }
            let s_idx = hsp_list.hsps[0].s_idx as usize;
            let subject = sam_subjects.len() as u32;
            let (id, offset, key) = match subjects {
                SubjectNames::Records(records) => (
                    SeqId::Local(&records[s_idx].local_id),
                    subject_placements.offset(s_idx) as i32,
                    s_idx as u64,
                ),
                // A sequence with parsed Seq-ids is named by an id of its Bioseq
                // (`ParsedSubject::sam_id`).
                SubjectNames::Database {
                    ordinals, parsed, ..
                } => (
                    match parsed.get(s_idx).and_then(Option::as_deref) {
                        Some(subject) => SeqId::Db(subject.sam_id()),
                        None => SeqId::BlOrdId(ordinals[s_idx]),
                    },
                    0,
                    u64::from(ordinals[s_idx]),
                ),
            };
            sam_subjects.push(SamSequence {
                record: SamId::Subject(SubjectRecord {
                    id,
                    searched_length: (lengths.searched)(s_idx) as i32,
                    offset,
                }),
                key,
                index: s_idx,
                length: (lengths.whole)(s_idx) as i32,
            });
            for hsp in &hsp_list.hsps {
                debug_assert_eq!(hsp.s_idx as usize, s_idx, "one subject per HSP list");
                hsp_records.push(HspRecord {
                    subject,
                    q_offset: hsp.internal_q_offset_0 as i32,
                    q_end: hsp.internal_q_end_0 as i32,
                    s_offset: hsp.internal_s_offset_0 as i32,
                    s_end: hsp.internal_s_end_0 as i32,
                    q_frame: if hsp.query_frame < 0 { -1 } else { 1 },
                    s_frame: 1,
                    score: hsp.raw_score,
                    num: 0,
                    evalue: hsp.e_value,
                    bit_score: hsp.bit_score,
                    num_ident: hsp.num_ident as i32,
                    num_positives: 0,
                    comp_adjust: 0,
                    ops: hsp.gap_info.as_deref().unwrap_or(&[]),
                });
            }
        }
        let query_record = queries.record(q_idx);
        let query = SamSequence {
            record: SamId::Query(query_record),
            key: (queries_before + q_idx) as u64,
            index: q_idx,
            length: queries
                .placements
                .length(q_idx, query_record.searched_length as usize) as i32,
        };
        document.query(&run, &query, &sam_subjects, &hsp_records, letters)?;
    }
    Ok(())
}

/// Where the letters of the subjects come from.
pub(crate) enum SubjectSource<'a> {
    /// A `-subject` set: the records (`SubjectLetters`, in record coordinates).
    Records(SubjectLetters<'a>),
    /// A database: letters `[from, to)` of the sequence of a rank (the range decoder of
    /// `seqdb`).
    Database(&'a dyn Fn(usize, usize, usize) -> Result<Vec<u8>>),
}

/// The letters of the queries and the subjects for SEQ.
pub(crate) struct BlastnSamLetters<'a> {
    /// The searched query records (with `-query_loc`, the interval's letters).
    pub queries: &'a [FastaRecord],
    pub query_placements: &'a Placements,
    pub subjects: SubjectSource<'a>,
}

impl SamLetters for BlastnSamLetters<'_> {
    fn letters(
        &mut self,
        row: SamRow,
        index: usize,
        from: usize,
        to: usize,
        out: &mut Vec<u8>,
    ) -> Result<()> {
        match row {
            SamRow::Query => {
                let offset = self.query_placements.offset(index);
                let record = self.queries[index].seq();
                let letters = from
                    .checked_sub(offset)
                    .zip(to.checked_sub(offset))
                    .and_then(|(from, to)| record.get(from..to))
                    .context("a SAM record lies outside the searched query letters")?;
                out.extend(letters.iter().map(u8::to_ascii_uppercase));
            }
            SamRow::Subject => match &self.subjects {
                SubjectSource::Records(letters) => {
                    out.extend_from_slice(&letters.upper_letters_uncounted(index, from, to))
                }
                SubjectSource::Database(read) => out.extend_from_slice(&read(index, from, to)?),
            },
        }
        Ok(())
    }
}
