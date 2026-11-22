//! BLASTN's Seq-annot of Seq-aligns (outfmt 8, 12): the report records of its final hit lists
//! (`report::record`), one query at a time, written by `report::serial::write_query`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1549-1649
//! The final hit lists are in that order already: the tabular formats walk the same
//! Seq-aligns (blast_format.cpp:772-800), and BLASTN writes its rows from these lists in this
//! order (`write_output_blastn_hitlists_to_writer`).
//!
//! What NCBI's `BlastHSP` holds, from LOSAT's `BlastnHsp`:
//! - `query.offset/end`: `internal_q_offset_0/internal_q_end_0` (in the query's context: the
//!   reverse complement for a minus-strand HSP), `subject.offset/end`:
//!   `internal_s_offset_0/internal_s_end_0`, all in the searched letters (the record
//!   coordinates that `shift_to_records` makes are not read);
//! - `query.frame`: `query_frame`; `subject.frame`: +1 (BLASTN's subject is the plus strand);
//! - `gap_info`: `gap_info` (`None` only for an empty script, `HspRecord::ops`);
//! - `num_ident`: `num_ident`; `num_positives`: 0, because NCBI counts positives only with a
//!   protein matrix (core/blast_hits.c:770-788), so BLASTN never writes them (LOSAT's
//!   `num_positives` holds the matches); `num` 0 and `comp_adjustment_method` 0 (no linking,
//!   no composition adjustment).
//! - the query length of the construction: the searched length (`query_data.GetSeqLength`,
//!   blast_seqalign.cpp:1972, with `-query_loc` the interval's).

use std::io::Write;

use anyhow::Result;

use super::hsp::BlastnHitList;
use crate::blastinput::app::Program;
use crate::blastinput::seq_range::Placements;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryNames, RunRecord, SubjectNames, SubjectRecord};
use crate::report::serial::{write_hit_list_seqaligns, SerialFormat};

/// Writes the documents of every query of `hit_lists` (indexed by query) to `writer`, each
/// query's warnings before it (`write_hit_list_seqaligns`). `gapped`: the search is not
/// `-ungapped` (an ungapped search's Seq-aligns are Dense-diags, one Seq-align per subject).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_blastn_seqalign(
    format: SerialFormat,
    gapped: bool,
    hit_lists: &[Option<BlastnHitList>],
    queries: &QueryNames<'_>,
    subjects: SubjectNames<'_>,
    subject_placements: &Placements,
    writer: &mut dyn Write,
    warnings: &mut QueryWarnings<'_>,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastn,
        rps: false,
        db_title: subjects.db_title(),
        gapped,
    };
    write_hit_list_seqaligns(
        format,
        &run,
        hit_lists,
        queries,
        writer,
        warnings,
        |_, hit_list, subject_records, hsp_records| {
            query_records(
                hit_list,
                &subjects,
                subject_placements,
                subject_records,
                hsp_records,
            );
        },
    )
}

/// The subject and HSP records of one query's hit list (what `write_query` and the archive
/// write), replacing the contents of `subject_records` and `hsp_records`.
pub(crate) fn query_records<'a>(
    hit_list: &'a BlastnHitList,
    subjects: &SubjectNames<'a>,
    subject_placements: &Placements,
    subject_records: &mut Vec<SubjectRecord<'a>>,
    hsp_records: &mut Vec<HspRecord<'a>>,
) {
    subject_records.clear();
    hsp_records.clear();
    for hsp_list in &hit_list.hsplist_array {
        if hsp_list.hsps.is_empty() {
            continue;
        }
        let s_idx = hsp_list.hsps[0].s_idx as usize;
        let subject = subject_records.len() as u32;
        // A database sequence's length is not read: BLASTN's subject row is the plus
        // strand (`s_GetAlignmentStart`).
        subject_records.push(subjects.record(s_idx, subject_placements.offset(s_idx), 0));
        for hsp in &hsp_list.hsps {
            debug_assert_eq!(hsp.s_idx as usize, s_idx, "one subject per HSP list");
            hsp_records.push(HspRecord {
                subject,
                q_offset: hsp.internal_q_offset_0,
                q_end: hsp.internal_q_end_0,
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
}
