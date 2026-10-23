//! TBLASTN's Seq-annot of Seq-aligns (outfmt 8, 12): the report records of its final hit lists
//! (`report::record`), one query at a time, written by `report::serial::write_query`.
//!
//! NCBI makes the Seq-aligns of each query from its final hit list
//! (`BlastHitList2SeqAlign_OMF`): subject by subject in hit-list order, each subject's HSPs
//! sorted by e-value, one Seq-align per HSP of score other than 0 with one Std-seg per
//! edit-script operation, the subject row translated. The report walks the lists after
//! `sort_hsps_for_report` and `drop_zero_score_hsps`, as the tabular rows do.
//!
//! What NCBI's `BlastHSP` holds, from LOSAT's `LinkedHsp` and its `KappaHspPayload`:
//! - `query.offset/end`: `q_start`, `q_end` (protein letters of the searched query, end
//!   exclusive); `subject.offset/end`: `s_start`, `s_end` (letters of the subject's frame);
//!   `query.frame` 0, `subject.frame`: `frame`;
//! - `score`, `num`, `evalue`: `hsp.score`, `num` (0 with `-sum_stats F`, which links nothing:
//!   NCBI's `num` is then 0 too), `evalue` (the sum e-value of a linked set); `bit_score`,
//!   `num_ident`, `num_positives`: the payload's engine values (`num_ident`, `num_positives`;
//!   not `report_num_*`, the tabular formatter's recount from the displayed rows);
//! - `comp_adjustment_method`: the payload's `matrix_adjust_rule` (0, 1 or 2);
//! - `gap_info`: the payload's `edit_script`;
//! - the lengths of the construction: the searched query's length, the subject's searched
//!   length (`SubjectLetters::length`, with `-subject_loc` the interval's).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:331-342
//! ```c
//!             switch (align->matrix_adjust_rule) {
//!             case eDontAdjustMatrix:
//!                 new_hsp->comp_adjustment_method = eNoCompositionBasedStats;
//!                 break;
//!             case eCompoScaleOldMatrix:
//!                 new_hsp->comp_adjustment_method = eCompositionBasedStats;
//!                 break;
//!             default:
//!                 new_hsp->comp_adjustment_method = eCompositionMatrixAdjust;
//!                 break;
//!             }
//! ```

use std::io::Write;

use anyhow::{Context, Result};

use super::stage_d_results::KappaResultHitList;
use super::stage_e_report::{ReportRanges, ShownSubjects};
use super::x_subject_store::SubjectLetters;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::{write_query, SerialFormat};

/// Writes the documents of every query of `hitlists` (indexed by query) to `writer`, each
/// query's warnings before it. `subject_records` and `subject_letters` are indexed by the hit
/// lists' subject index (a `-db` search's reported sequences by rank, `shown`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_tblastn_seqalign(
    format: SerialFormat,
    writer: &mut dyn Write,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    hitlists: &[KappaResultHitList],
    ranges: &ReportRanges<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastn,
        db_title: shown.db_title(),
        gapped: true,
    };
    let mut subjects: Vec<SubjectRecord<'_>> = Vec::new();
    let mut hsps: Vec<HspRecord<'_>> = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    for (q_idx, (query, hitlist)) in query_records.iter().zip(hitlists).enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        query_records_of(
            hitlist,
            subject_records,
            subject_letters,
            shown,
            ranges,
            &mut subjects,
            &mut hsps,
        )?;
        let query = QueryRecord {
            id: SeqId::Local(&query.local_id),
            searched_length: i32::try_from(query.seq().len())?,
            offset: i32::try_from(ranges.queries.offset(q_idx))?,
        };
        buf.clear();
        if write_query(format, &run, &query, &subjects, &hsps, &mut buf)? {
            writer.write_all(&buf)?;
        }
    }
    Ok(())
}

/// The subject and HSP records of one query's final hit list (the elements of its Seq-annot in
/// outfmt 8, and of the `alignments` of the archive of outfmt 11), appended to the cleared
/// `subjects` and `hsps`.
pub(super) fn query_records_of<'a>(
    hitlist: &'a KappaResultHitList,
    subject_records: &'a [FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    ranges: &ReportRanges<'_>,
    subjects: &mut Vec<SubjectRecord<'a>>,
    hsps: &mut Vec<HspRecord<'a>>,
) -> Result<()> {
    subjects.clear();
    hsps.clear();
    for list in hitlist.lists() {
        if list.hsps.hsps.is_empty() {
            continue;
        }
        let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
        let subject = subjects.len() as u32;
        subjects.push(SubjectRecord {
            id: match shown.ordinal(oid) {
                Some(ordinal) => SeqId::BlOrdId(ordinal),
                None => SeqId::Local(&subject_records[oid].local_id),
            },
            searched_length: i32::try_from(subject_letters.length(oid))?,
            // A database's sequences have no subject remap (`RunRecord::remaps_subjects`).
            offset: i32::try_from(ranges.subjects.offset(oid))?,
        });
        for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
            let hsp = &linked.hsp;
            hsps.push(HspRecord {
                subject,
                q_offset: hsp.q_start,
                q_end: hsp.q_end,
                s_offset: hsp.s_start,
                s_end: hsp.s_end,
                q_frame: 0,
                s_frame: hsp.frame,
                score: hsp.score,
                num: linked.num,
                evalue: linked.evalue,
                bit_score: payload.bit_score,
                num_ident: payload.num_ident,
                num_positives: i32::try_from(payload.num_positives)?,
                comp_adjust: match payload.matrix_adjust_rule {
                    EMatrixAdjustRule::DontAdjustMatrix => 0,
                    EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                    _ => 2,
                },
                ops: &payload.edit_script,
            });
        }
    }
    Ok(())
}
