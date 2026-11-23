//! Parts of `post_process_hits_and_write` (the report of one call): the refusals of the query-anchored
//! formats (outfmt 1 to 4), the SAM writer (outfmt 17), the structured documents (outfmt 5, 13 to 16)
//! and the subjects' names that the Seq-align, archive and SAM writers take.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-755 (the XML, XML2 and
//! SAM branches of `CBlastFormat::x_PrintStructuredReport`, cited in `write_structured_xml` and
//! `write_sam`). The refusals have no NCBI counterpart (see `refuse_anchored_reports`).
//!
//! This module is a child of `run` for what `run` keeps private (the types `BlastnReportInputs` and
//! `SamTarget`, the module `database`). Each function holds the statements of one part of
//! `post_process_hits_and_write`, in their order, and reads what that part read.

use super::{
    database, BlastnHitList, BlastnOutputFormat, BlastnReportInputs, SamTarget, SubjectLetters,
};
use crate::algorithm::blastn::pairwise::shown_query_masks;
use crate::algorithm::blastn::sam;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::query_batch::NO_DATA_MESSAGE;
use crate::report::pairwise::PairwiseHit;
use crate::report::record::{QueryNames, SubjectNames};
use anyhow::Result;
use std::sync::Arc;

/// Refuses the call, before any format of it is written, when one of its query-anchored reports
/// (outfmt 1 to 4) is a report that NCBI BLAST+ aborts: a query whose own record is one of the
/// subjects (`-subject`) and a reported database sequence that NCBI's scope finds as a query of its
/// batch. `rows_of` renders the hits of a range of queries; `pairwise_hits` holds those of every
/// query when the call made them already.
///
/// No NCBI counterpart: the refusal is LOSAT's (PD-LOSAT-NCBI-DEFECTS rule 3); NCBI writes bytes before
/// the abort that depend on its buffers. It runs once per call, after `subject_titles` is made.
pub(super) fn refuse_anchored_reports(
    output_formats: &[BlastnOutputFormat],
    report: &BlastnReportInputs<'_>,
    hit_lists: &[Option<BlastnHitList>],
    pairwise_hits: Option<&[PairwiseHit]>,
    subject_titles: &[Arc<[u8]>],
    rows_of: &dyn Fn(std::ops::Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    // The query-anchored formats (outfmt 1 to 4, line RG T3b), before any format of the call is
    // written, as BLASTP's (`blast_engine::run_call`): a query whose own record is one of the
    // subjects (`-subject`, `-parse_deflines`) and whose section NCBI BLAST+ aborts on
    // (`pairwise::section_aborts`, PD-LOSAT-NCBI-DEFECTS rule 3), and a reported database sequence
    // that NCBI's scope finds as a query of its batch, by its Seq-align id (a generated
    // `lcl|Query_N`) or, in a description line, by its display id
    // (`common::database::refuse_query_ids_of_anchored_reports`): NCBI shows the query in its
    // place. NCBI's bytes before an abort depend on its buffers, so LOSAT writes none of the call's.
    let anchored_view = output_formats.iter().find_map(|format| match format {
        BlastnOutputFormat::Pairwise(view)
            if *view != crate::report::pairwise::PairwiseView::Pairwise =>
        {
            Some(*view)
        }
        _ => None,
    });
    if let Some(view) = anchored_view {
        let sort = report.formatting.report_sort();
        let rows_of_query = |q_idx: usize| -> Result<std::borrow::Cow<'_, [PairwiseHit]>> {
            Ok(match pairwise_hits {
                Some(hits) => std::borrow::Cow::Borrowed(hits),
                None => std::borrow::Cow::Owned(rows_of(q_idx..q_idx + 1)?),
            })
        };
        if report.database_reader.is_none() {
            for (q_idx, record) in report.queries.iter().enumerate() {
                if record.believed().is_none_or(|ids| ids.shared.is_none()) {
                    continue;
                }
                let is_query = |s_idx: u32| {
                    report
                        .subjects
                        .records()
                        .get(s_idx as usize)
                        .is_some_and(|subject| {
                            crate::report::query_anchored::subject_is_query_record(record, subject)
                        })
                };
                let rows = rows_of_query(q_idx)?;
                let aborts = crate::report::pairwise::section_aborts(
                    view,
                    sort,
                    report.num_alignments,
                    q_idx,
                    &rows,
                    subject_titles,
                    is_query,
                    |group| {
                        crate::report::pairwise::blastn_section_hsps(
                            report.num_alignments,
                            group,
                            &[],
                            is_query,
                        )
                    },
                )?;
                if aborts {
                    return Err(crate::report::query_anchored::aborted_section(
                        record, "BLASTN",
                    ));
                }
            }
        }
        if let (Some(ordinals), Some(reader)) = (report.database_ordinals, report.database_reader) {
            let batches = crate::blastinput::query_batch::AnchoredBatchScopes::new(
                report.queries,
                if report.batch_spans.is_empty() {
                    vec![0..report.queries.len()]
                } else {
                    report.batch_spans.to_vec()
                },
            );
            let oids = reader.oid_map();
            crate::algorithm::common::database::refuse_query_ids_of_anchored_reports(
                hit_lists.len(),
                |q_idx| {
                    hit_lists[q_idx]
                        .iter()
                        .flat_map(|hit_list| &hit_list.hsplist_array)
                        .filter_map(|list| list.hsps.first().map(|hsp| hsp.s_idx))
                        .collect()
                },
                |q_idx| {
                    let rows = rows_of_query(q_idx)?;
                    Ok(crate::report::pairwise::blastp_described_subjects(
                        view.report_sort(sort),
                        None,
                        report.num_descriptions,
                        q_idx,
                        &rows,
                        |previous, s_idx| {
                            crate::report::db_subject::rank_joins_as_copy(
                                report.database_parsed,
                                previous as usize,
                                s_idx as usize,
                            )
                        },
                    ))
                },
                report.database_parsed,
                |rank| oids.oid(ordinals[rank as usize]),
                |q_idx| batches.of_query(q_idx),
                batches.generated(),
                "BLASTN",
            )?;
        }
    }
    Ok(())
}

/// Adds the queries of the call to the run's SAM document (outfmt 17) and, in the call that ends the
/// run, writes it. The posts of the queries (`before_query`) come first, then the document is filled,
/// then `finish`. Its NCBI counterpart, the SAM branch of
/// `CBlastFormat::x_PrintStructuredReport`, is cited at the top of the body.
pub(super) fn write_sam(
    sam_target: &mut Option<SamTarget<'_>>,
    report: &BlastnReportInputs<'_>,
    hit_lists: &[Option<BlastnHitList>],
    query_lengths: &[usize],
    format_warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    writer: &mut crate::cli::NcbiOstream<'_>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:750-755
    // The queries' records are added to the run's document, which NCBI writes when its
    // formatter is destroyed (`SamDocument::finish`): after the last call, or in a call
    // that is the run's only one.
    let target = sam_target
        .as_mut()
        .expect("a SAM document is made for outfmt 17");
    if let Some(warnings) = format_warnings {
        for q_idx in 0..hit_lists.len() {
            warnings.before_query(q_idx, writer)?;
        }
    }
    let names = QueryNames {
        records: report.queries,
        searched_lengths: query_lengths,
        placements: &report.ranges.queries.placements,
    };
    match (
        report.database_ordinals,
        report.database,
        report.database_reader,
    ) {
        (Some(ordinals), Some(_), Some(reader)) => {
            // A database's sequences are searched whole; the lengths come from its
            // index (no letters).
            let lengths = reader.sam_lengths(ordinals)?;
            let cache = database::SamLetterCache::default();
            let letters_of = |rank: usize, from: usize, to: usize| {
                reader.sam_letters(&cache, ordinals[rank], from, to)
            };
            sam::add_blastn_sam(
                target.document,
                report.ungapped,
                hit_lists,
                &names,
                report.queries_before,
                subject_names_of(report),
                report.ranges.subjects.placements,
                &sam::SubjectLengths {
                    searched: &|rank| lengths[rank],
                    whole: &|rank| lengths[rank],
                },
                &mut sam::BlastnSamLetters {
                    queries: report.queries,
                    query_placements: &report.ranges.queries.placements,
                    subjects: sam::SubjectSource::Database(&letters_of),
                },
            )?;
        }
        _ => {
            // The records as read with `-subject_loc` (never packed), else the
            // searched ones: the letters are asked for in record coordinates.
            let whole_records = match report.ranges.subjects.records {
                Some(records) => SubjectLetters::new(records, None),
                None => report.subjects,
            };
            sam::add_blastn_sam(
                target.document,
                report.ungapped,
                hit_lists,
                &names,
                report.queries_before,
                SubjectNames::Records(report.subjects.records()),
                report.ranges.subjects.placements,
                &sam::SubjectLengths {
                    searched: &|s_idx| report.subjects.len(s_idx),
                    whole: &|s_idx| {
                        report
                            .ranges
                            .subjects
                            .placements
                            .length(s_idx, report.subjects.len(s_idx))
                    },
                },
                &mut sam::BlastnSamLetters {
                    queries: report.queries,
                    query_placements: &report.ranges.queries.placements,
                    subjects: sam::SubjectSource::Records(whole_records),
                },
            )?;
        }
    }
    if target.finish {
        target.document.finish(writer)?;
    }
    Ok(())
}

/// What one structured document of the call (outfmt 5, 13 to 16) is: its kind and what the call's
/// other formats and the caller give it.
pub(super) struct XmlDocument<'a> {
    /// The XML2 / JSON kind of outfmt 13 to 16; `None` for outfmt 5 (BLAST XML).
    pub(super) kind: Option<crate::report::blast_xml2::Xml2Kind>,
    /// Every structured format of the call is an XML2 one: a `-subject` document holds back
    /// (`blast_xml::write_blastn_xml`).
    pub(super) only_xml2: bool,
    /// The caller's receiver of the side files of outfmt 13 and 14, if it gives one.
    pub(super) side_sink: Option<&'a crate::report::blast_xml2::SideFileSink>,
}

/// Outfmt 5 (BLAST XML) and 13 to 16 (the XML2 and JSON documents): the structured report of the
/// queries of the call (`blast_xml::write_blastn_xml`), the side files of outfmt 13 and 14 into the
/// caller's receiver where it gives one. `unsearched` marks the queries of the batches that were
/// not searched; `subject_title` is the database's title.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-710
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:745-749
/// The BLAST XML, XML2 and JSON branches of `CBlastFormat::x_PrintStructuredReport`.
pub(super) fn write_structured_xml(
    structured: XmlDocument<'_>,
    report: &BlastnReportInputs<'_>,
    hit_lists: &[Option<BlastnHitList>],
    query_lengths: &[usize],
    unsearched: &[bool],
    subject_title: &str,
    writer: &mut crate::cli::NcbiOstream<'_>,
    format_warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> Result<()> {
    // The XML2 / JSON documents of outfmt 13 to 16 take the path of outfmt 5: the same
    // per-query records, written by another writer.
    let XmlDocument {
        kind: xml2_kind,
        only_xml2: xml_only_xml2,
        side_sink,
    } = structured;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1444
    // The CBlastOutput document of the queries (`blast_xml`).
    let is_database = report.database.is_some();
    let xml_stats = super::blastn_xml_stats(report, query_lengths)?;
    let xml_messages: Vec<Vec<String>> = report
        .queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| {
            let mut messages = Vec::new();
            if report.context_karlin[2 * q_idx + report.first_context].is_none() {
                messages.push(if query.seq().is_empty() {
                    NO_DATA_MESSAGE.to_string()
                } else {
                    crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()
                });
            }
            if !unsearched[q_idx]
                && report
                    .database_reader
                    .is_some_and(|database| database.filtered_empty())
            {
                messages.push(crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string());
            }
            messages
        })
        .collect();
    let version = format!("BLASTN {}", crate::wording::REPORT_VERSION);
    let filter = crate::report::blast_xml::filter::filter_string(
        report.dust,
        None,
        report.repeat_filter,
        report.window_masker,
        report.soft_masking,
    );
    let subject_records = report.subjects.records();
    let database_length = |rank: usize| report.subjects.len(rank);
    let record_length = |s_idx: usize| {
        report
            .ranges
            .subjects
            .placements
            .length(s_idx, report.subjects.len(s_idx))
    };
    // The taxids of the reported sequences and the taxonomy database's names, read for
    // the XML2 and JSON formats only.
    let (xml_taxids, xml_taxonomy, xml_parsed_taxonomy, xml_taxdb_notes): (
        Vec<i32>,
        Vec<crate::report::blast_xml::HitTaxonomy>,
        Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
        crate::algorithm::common::database::TaxdbNotes,
    ) = match (xml2_kind, report.database_ordinals, report.database_reader) {
        (Some(_), Some(ordinals), Some(database)) if is_database => database.read_taxids(
            ordinals,
            report.database_parsed,
            hit_lists,
            report.database_loads,
        )?,
        _ => (Vec::new(), Vec::new(), Vec::new(), Default::default()),
    };
    // The masked regions of the queries in record coordinates and the queries of the
    // batches that were not searched, which only XML2 shows (outfmt 5 shows neither).
    let (xml_masks, xml_searched): (Vec<Vec<crate::report::blast_xml::QueryMask>>, Vec<bool>) =
        if xml2_kind.is_some() {
            (
                (0..report.queries.len())
                    .map(|q_idx| {
                        shown_query_masks(
                            &report.query_masks[q_idx],
                            query_lengths[q_idx],
                            report.ranges.queries.placements.offset(q_idx),
                        )
                        .iter()
                        .map(|mask| crate::report::blast_xml::QueryMask {
                            from: mask.start as i32,
                            to: mask.end as i32 - 1,
                            frame: 0,
                        })
                        .collect()
                    })
                    .collect(),
                unsearched.iter().map(|&skipped| !skipped).collect(),
            )
        } else {
            (Vec::new(), Vec::new())
        };
    let xml_subjects = match (report.database_ordinals, report.database_reader) {
        (Some(ordinals), Some(reader)) if is_database => {
            crate::report::blast_xml::XmlSubjects::Database {
                ordinals,
                oids: reader.oid_map(),
                titles: crate::report::blast_xml::XmlTitles::Records(subject_records),
                // Only XML2 shows the taxid.
                taxids: &xml_taxids,
                taxonomy: &xml_taxonomy,
                length: &database_length,
                parsed: report.database_parsed,
                parsed_taxonomy: &xml_parsed_taxonomy,
                loads: report.database_loads,
            }
        }
        _ => crate::report::blast_xml::XmlSubjects::Records {
            records: subject_records,
            length: &record_length,
            protein: false,
        },
    };
    let xml = crate::algorithm::blastn::blast_xml::BlastnXml {
        run: crate::report::blast_xml::XmlRun {
            program: "blastn",
            version: &version,
            reference: if report.megablast {
                crate::report::blast_xml::REFERENCE_MEGABLAST
            } else {
                crate::report::blast_xml::REFERENCE_GAPPED_BLAST
            },
            db_name: if is_database { subject_title } else { "" },
            matrix: None,
            expect: report.expect,
            gap_open: report.gap_open,
            gap_extend: report.gap_extend,
            reward_penalty: Some((report.reward, report.penalty)),
            filter: &filter,
            db_num: report.database.map_or(0, |database| database.num_sequences),
            db_len: report.database.map_or(0, |database| database.total_letters),
            db_len_overflow: report
                .database
                .is_some_and(|database| database.total_letters_overflow),
            positives: None,
            mask_char: b'N',
        },
        queries: report.queries,
        searched_lengths: query_lengths,
        stats: xml_stats,
        messages: xml_messages,
        masks: xml_masks,
        searched: xml_searched,
        subjects: xml_subjects,
        taxdb_notes: xml_taxdb_notes,
        letters: report.subjects,
        epilog: report.epilog,
        only_xml2: xml_only_xml2,
    };
    match xml2_kind {
        None => crate::algorithm::blastn::blast_xml::write_blastn_xml(
            &xml,
            &mut crate::report::blast_xml::BlastXmlWriter::new(report.queries_before),
            hit_lists,
            writer,
            format_warnings,
        )?,
        Some(kind) => {
            // The gap costs are written for a gapped search (NCBI reads them from
            // the options when `GetGappedMode()`); `cbs` and the genetic codes
            // are not those of a nucleotide search (blastxml2_format.cpp:531-546).
            let params = crate::report::blast_xml2::Xml2Run {
                target: if is_database {
                    crate::report::blast_xml2::Xml2Target::Database
                } else {
                    crate::report::blast_xml2::Xml2Target::Subjects(subject_records.len())
                },
                subject_names: if is_database {
                    None
                } else {
                    crate::report::blast_xml2::subject_names(
                        subject_records.iter().map(FastaRecord::believed),
                    )
                },
                gapped: !report.ungapped,
                cbs: 0,
                query_gencode: None,
                db_gencode: None,
            };
            let mut document =
                crate::report::blast_xml2::Xml2Writer::for_run(kind, params, report.out, side_sink);
            crate::algorithm::blastn::blast_xml::write_blastn_xml(
                &xml,
                &mut document,
                hit_lists,
                writer,
                format_warnings,
            )?;
        }
    }
    Ok(())
}

/// The names of the reported subjects as the Seq-align, archive and SAM writers take them: a
/// database's ordinals with its OID map, definition and parsed ids, else the subject records.
///
/// No NCBI counterpart: NCBI reads the subjects' names from the Seq-align ids and the database alone;
/// this is where the Seq-align, archive and SAM database writers pick which of LOSAT's two
/// sources a call has (the SAM writer's subject-set arm names the records directly).
pub(super) fn subject_names_of<'a>(report: &'a BlastnReportInputs<'_>) -> SubjectNames<'a> {
    match (
        report.database_ordinals,
        report.database,
        report.database_reader,
    ) {
        (Some(ordinals), Some(database), Some(reader)) => SubjectNames::Database {
            ordinals,
            oids: reader.oid_map(),
            title: &database.definition,
            parsed: report.database_parsed,
        },
        _ => SubjectNames::Records(report.subjects.records()),
    }
}
