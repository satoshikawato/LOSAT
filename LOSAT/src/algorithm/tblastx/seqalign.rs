//! TBLASTX's Seq-annot of Seq-aligns (outfmt 8, 12): the report records of the final HSP list
//! (`report::record`), one query at a time, written by `report::serial::write_query`.
//!
//! TBLASTX is an ungapped search, so NCBI makes one `diags` Seq-align per subject of a query's
//! final hit list (`BLASTUngappedHspListToSeqAlign`), one Std-seg per HSP with the HSP's scores,
//! subject by subject in hit-list order, each subject's HSPs sorted by e-value. The final list
//! (`hit_collector.rs`) holds the HSPs in that order, query by query, one run of HSPs per subject
//! (the order of the tabular rows).
//!
//! What NCBI's `BlastHSP` holds, from LOSAT's `TblastxHsp`:
//! - `query.offset/end`, `query.frame`: `sort_query_offset`, `sort_query_end` (letters of the
//!   query's frame in the searched letters, end exclusive), `query_frame`;
//!   `subject.offset/end`, `subject.frame`: `sort_subject_offset`, `sort_subject_end`,
//!   `subject_frame` (these are not moved by `shift_to_records`, which moves the tabular
//!   coordinates only);
//! - `score`, `num`, `evalue`, `bit_score`: `raw_score`, `num`, `e_value` (the sum e-value of
//!   the linked set), `bit_score`;
//! - `num_ident`, `num_positives`: the engine's counts (`report::engine_counts`, made before
//!   `set_displayed_identities` replaces `num_ident` by the display count);
//! - no composition adjustment, no edit script;
//! - the lengths of the construction: the searched lengths of the query and the subject
//!   (with `-query_loc` / `-subject_loc`, of the intervals), whose starts the remaps add.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1606-1626
//! ```c++
//!         GetFilteredRedundantSeqids(*seqinfo_src, hsp_list->oid, seqid_list, subject_id->IsGi());
//!         // stores a CSeq_align for each matching sequence
//!         vector<CRef<CSeq_align > > hit_align;
//!         if (is_gapped) {
//!                 BLASTHspListToSeqAlign(prog,
//!                                        hsp_list,
//!                                        query_id,
//!                                        subject_id,
//!                                        query_length,
//!                                        subj_length,
//!                                        is_ooframe,
//!                                        seqid_list,
//!                                        hit_align);
//!         } else {
//!                 BLASTUngappedHspListToSeqAlign(prog,
//!                                                hsp_list,
//!                                                query_id,
//!                                                subject_id,
//!                                                query_length,
//!                                                subj_length,
//!                                                seqid_list,
//! ```

use std::io::Write;

use anyhow::Result;

use super::blast_engine::TblastxHsp;
use super::report::SubjectLetters;
use crate::blastinput::app::Program;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::Placements;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::{write_query_to, SerialFormat};

/// The subjects of a TBLASTX search as the Seq-aligns name them: the records of `-subject`
/// (`Subject_N`) with their `-subject_loc` placements, or a database's reported sequences by
/// rank (`ordinals`) with the database's title, named by their Seq-align ids where they have
/// parsed Seq-ids (`parsed`, by rank; empty without them).
pub(crate) enum SeqalignSubjects<'a> {
    Records {
        records: &'a [FastaRecord],
        placements: &'a Placements,
    },
    Database {
        ordinals: &'a [u32],
        title: &'a [u8],
        parsed: &'a [Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>],
    },
}

/// The report records of one query's HSPs (`hits`, in report order, with their engine counts):
/// one subject record per run of HSPs of a subject and one HSP record per HSP, replacing the
/// contents of `subject_records` and `hsp_records`. The Seq-annot formats (outfmt 8, 9, 12) and
/// the BLAST archive (outfmt 11) write the Seq-aligns from them.
pub(crate) fn query_records_of<'a>(
    hits: &[TblastxHsp],
    engine: &[(u32, u32)],
    subjects: &'a SeqalignSubjects<'a>,
    subject_letters: SubjectLetters<'_>,
    subject_records: &mut Vec<SubjectRecord<'a>>,
    hsp_records: &mut Vec<HspRecord<'a>>,
) -> Result<()> {
    subject_records.clear();
    hsp_records.clear();
    let mut last_subject: Option<u32> = None;
    for (hsp, &(num_ident, num_positives)) in hits.iter().zip(engine) {
        let hit = &hsp.hit;
        if last_subject != Some(hit.s_idx) {
            last_subject = Some(hit.s_idx);
            let s_idx = hit.s_idx as usize;
            subject_records.push(match subjects {
                SeqalignSubjects::Records {
                    records,
                    placements,
                } => SubjectRecord {
                    id: records[s_idx].seq_id(),
                    searched_length: i32::try_from(subject_letters.length(s_idx))?,
                    offset: i32::try_from(placements.offset(s_idx))?,
                    use_this_seqid: &[],
                },
                SeqalignSubjects::Database {
                    ordinals, parsed, ..
                } => {
                    let subject = parsed.get(s_idx).and_then(Option::as_deref);
                    SubjectRecord {
                        id: match subject {
                            Some(subject) => SeqId::Db(subject.alignment_id()),
                            None => SeqId::BlOrdId(ordinals[s_idx]),
                        },
                        searched_length: i32::try_from(subject_letters.length(s_idx))?,
                        offset: 0,
                        use_this_seqid: crate::report::db_subject::use_this_seqid(subject),
                    }
                }
            });
        }
        hsp_records.push(HspRecord {
            subject: subject_records.len() as u32 - 1,
            q_offset: i32::try_from(hit.sort_query_offset)?,
            q_end: i32::try_from(hit.sort_query_end)?,
            s_offset: i32::try_from(hit.sort_subject_offset)?,
            s_end: i32::try_from(hit.sort_subject_end)?,
            q_frame: i8::try_from(hit.query_frame)?,
            s_frame: hsp.subject_frame,
            score: hit.raw_score,
            num: hsp.num,
            evalue: hit.e_value,
            bit_score: hit.bit_score,
            num_ident: i32::try_from(num_ident)?,
            num_positives: i32::try_from(num_positives)?,
            comp_adjust: 0,
            ops: &[],
        });
    }
    Ok(())
}

/// Writes the documents of every query to `writer`, each query's warnings before it. `hits`
/// is the final list in report order and `engine` its engine counts by HSP index
/// (`report::engine_counts`).
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_tblastx_seqalign(
    format: SerialFormat,
    writer: &mut dyn Write,
    hits: &[TblastxHsp],
    engine: &[(u32, u32)],
    queries: &[FastaRecord],
    query_placements: &Placements,
    subjects: &SeqalignSubjects<'_>,
    subject_letters: SubjectLetters<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let run = RunRecord {
        program: Program::Tblastx,
        db_title: match subjects {
            SeqalignSubjects::Records { .. } => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped: false,
    };
    let mut subject_records: Vec<SubjectRecord<'_>> = Vec::new();
    let mut hsp_records: Vec<HspRecord<'_>> = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    let mut next = 0usize;
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let start = next;
        while next < hits.len() && hits[next].hit.q_idx as usize == q_idx {
            next += 1;
        }
        query_records_of(
            &hits[start..next],
            &engine[start..next],
            subjects,
            subject_letters,
            &mut subject_records,
            &mut hsp_records,
        )?;
        let query_record = QueryRecord {
            id: query.seq_id(),
            searched_length: i32::try_from(query.seq().len())?,
            offset: i32::try_from(query_placements.offset(q_idx))?,
        };
        write_query_to(
            format,
            &run,
            &query_record,
            &subject_records,
            &hsp_records,
            &mut buf,
            &mut *writer,
        )?;
    }
    debug_assert_eq!(next, hits.len(), "the final list is in query order");
    Ok(())
}
