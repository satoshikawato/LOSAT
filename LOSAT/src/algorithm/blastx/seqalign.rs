//! BLASTX's Seq-annot of Seq-aligns (outfmt 8, 12): the report records of each query batch's
//! final hit lists (`report::record`), one query at a time, written by
//! `report::serial::write_query`.
//!
//! NCBI makes the Seq-aligns of each query from its final hit list
//! (`BlastHitList2SeqAlign_OMF`): subject by subject in hit-list order, each subject's HSPs
//! sorted by e-value (`render_queries` sorts them, `HspList::sort_evalue`). A gapped search
//! makes one Seq-align per HSP of score other than 0, with one Std-seg per edit-script
//! operation and the query row translated; an ungapped search (`-ungapped`) one `diags`
//! Seq-align per subject, one Std-seg per HSP.
//!
//! What NCBI's `BlastHSP` holds, from LOSAT's `Hsp`:
//! - `query.offset/end`, `query.frame`: `hsp.q_start`, `hsp.q_end` (letters of the query's
//!   frame, end exclusive), `hsp.frame`; `subject.offset/end`: `hsp.s_start`, `hsp.s_end`,
//!   `subject.frame` 0;
//! - `score`, `num`, `evalue`, `bit_score`: `hsp.score`, `num`, `evalue` (the sum e-value of a
//!   linked set), `bit_score`;
//! - `num_ident`, `num_positives`: `identity`, `positive` (`statistics::identity_and_positive`:
//!   NCBI's `Blast_HSPGetNumIdentitiesAndPositives` on the unmasked query frame, positives by
//!   the match's composition-adjusted matrix, `RedoneHsp::identity_matrix`), not the
//!   tabular formatter's recount from the displayed rows;
//! - `comp_adjustment_method`: `composition_method`; `gap_info`: `edit_script` (empty for an
//!   ungapped HSP);
//! - the lengths of the construction: the query's nucleotide length (LOSAT's BLASTX has no
//!   `-query_loc`) and the subject's length; no subject range (`-subject_loc` is rejected).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:518-531
//! ```c
//!             query = query_blk->sequence +
//!                 query_info->contexts[hsp->context].query_offset;
//!             query_nomask = query_blk->sequence_nomask +
//!                 query_info->contexts[hsp->context].query_offset;
//!         }
//!
//!         /* Translate subject if needed. */
//!         if (program_number == eBlastTypeTblastn) {
//!             ...
//!         }
//!         else
//!             status = Blast_HSPGetNumIdentitiesAndPositives(query_nomask, subject, hsp, scoring_options, 0, sbp);
//! ```

use std::io::Write;

use anyhow::Result;

use super::database::ReportedSubjects;
use super::input::FastaRecord;
use super::runtime::BatchResults;
use crate::blastinput::app::Program;
use crate::report::record::{HspRecord, QueryRecord, RunRecord, SubjectRecord};
use crate::report::seqid::SeqId;
use crate::report::serial::{write_query_to, SerialFormat};

/// The subjects of a BLASTX search as the Seq-aligns name them: the records of `-subject`
/// (`Subject_N`), or a database's reported sequences by rank with the database's title.
pub(super) enum SeqalignSubjects<'a> {
    Records,
    Database {
        reported: &'a ReportedSubjects,
        title: &'a [u8],
    },
}

/// Writes the documents of every query of `batches` to `writer`. A structured format writes
/// none of a query's own warnings (`ReportFormat::writes_query_warnings`), so nothing goes
/// to standard error here.
pub(super) fn write_blastx_seqalign(
    format: SerialFormat,
    writer: &mut impl Write,
    gapped: bool,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    names: &SeqalignSubjects<'_>,
    batches: &[BatchResults],
) -> Result<()> {
    let run = RunRecord {
        program: Program::Blastx,
        db_title: match names {
            SeqalignSubjects::Records => None,
            SeqalignSubjects::Database { title, .. } => Some(title),
        },
        gapped,
    };
    let mut subject_records: Vec<SubjectRecord<'_>> = Vec::new();
    let mut hsp_records: Vec<HspRecord<'_>> = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            query_records(
                lists,
                subjects,
                names,
                &mut subject_records,
                &mut hsp_records,
            )?;
            let query_record = QueryRecord {
                id: query.seq_id(),
                searched_length: i32::try_from(query.sequence.len())?,
                offset: 0,
            };
            write_query_to(
                format,
                &run,
                &query_record,
                &subject_records,
                &hsp_records,
                &mut buf,
                writer,
            )?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// The subject and HSP records of one query's hit lists (what `write_query` and the archive
/// write), replacing the contents of `subject_records` and `hsp_records`.
pub(super) fn query_records<'a>(
    lists: &'a [super::results::HspList],
    subjects: &'a [FastaRecord],
    names: &SeqalignSubjects<'a>,
    subject_records: &mut Vec<SubjectRecord<'a>>,
    hsp_records: &mut Vec<HspRecord<'a>>,
) -> Result<()> {
    subject_records.clear();
    hsp_records.clear();
    for list in lists {
        if list.hsps.is_empty() {
            continue;
        }
        let subject_index = subject_records.len() as u32;
        let subject = &subjects[list.oid];
        subject_records.push(SubjectRecord {
            id: match names {
                SeqalignSubjects::Records => subject.seq_id(),
                SeqalignSubjects::Database { reported, .. } => {
                    match reported.parsed_subject(list.oid) {
                        Some(subject) => SeqId::Db(subject.alignment_id()),
                        None => SeqId::BlOrdId(reported.ordinal(list.oid)),
                    }
                }
            },
            searched_length: i32::try_from(subject.sequence.len())?,
            offset: 0,
        });
        for hsp in &list.hsps {
            hsp_records.push(HspRecord {
                subject: subject_index,
                q_offset: hsp.hsp.q_start,
                q_end: hsp.hsp.q_end,
                s_offset: hsp.hsp.s_start,
                s_end: hsp.hsp.s_end,
                q_frame: hsp.hsp.frame,
                s_frame: 0,
                score: hsp.hsp.score,
                num: hsp.num,
                evalue: hsp.evalue,
                bit_score: hsp.bit_score,
                num_ident: i32::try_from(hsp.identity)?,
                num_positives: i32::try_from(hsp.positive)?,
                comp_adjust: u8::try_from(hsp.composition_method)?,
                ops: &hsp.edit_script,
            });
        }
    }
    Ok(())
}
