//! TBLASTN's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit lists
//! (`report::blast_xml`), one query at a time, once per run.
//!
//! NCBI makes the Iteration of each query from its Seq-aligns (`BlastHitList2SeqAlign_OMF`,
//! see `seqalign`); the XML walks the lists as the Seq-annot formats do: subject by subject in
//! hit-list order, each subject's HSPs sorted by e-value (`render` sorts them), the HSPs of
//! score 0 dropped. What NCBI's `s_SeqAlignSetToXMLHsps` reads of an HSP is its Seq-align
//! (the coordinates in the records, the score, the bit score, the e-value or the sum e-value
//! of a linked set) and the two rows of the alignment, from which it counts identities and
//! positives again: nothing of the engine's counts and nothing of outfmt 0's midline. The rows
//! are the ones outfmt 0 shows (`aligned_sequences`, the subject translated again with B, Z
//! and J for the ambiguous codons) without the lower case of the query's masks; the masks are
//! applied to the query row as NCBI's `s_MaskQuerySeq` does (`report::blast_xml::hsp`).
//!
//! The XML2 and JSON documents of outfmt 15 and 16 (`report::blast_xml2`) are written from the
//! same per-query records by `Xml2Writer`: `write` is generic over the writer
//! (`StructuredXml`). They add the masks of every query (also one without hits), the
//! `searched` flag of the query's batch, the taxid of a database sequence, and `cbs`,
//! `db-gencode` and the gap costs of `XmlValues`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
//! ```c++
//!     } else if (m_FormatType == CFormattingArgs::eXml) {
//!         CRef<CSearchResults> res(const_cast<CSearchResults*>(&results));
//!         res->TrimSeqAlign(m_HitlistSize);
//!         m_AccumulatedResults.push_back(res);
//!         ...
//!         	CCmdLineBlastXMLReportData report_data(m_AccumulatedQueries,
//!                                                    m_AccumulatedResults,
//!                                                    *m_Options, m_DbName, m_DbIsAA,
//!                                                    m_QueryGenCode, m_DbGenCode,
//!                                                    m_IsRemoteSearch);
//!         	BlastXML_FormatReport(xml_output, &report_data, &m_Outfile,
//!                                   m_BlastXMLIncremental.GetPointer());
//! ```
//! (`PrintOneResultSet` is called once per query's result set, the accumulated lists hold that
//! one query.)

use std::io::Write;

use anyhow::{ensure, Context, Result};

use super::database::XmlTaxa;
use super::search_gapped::TargetTranslation;
use super::search_seed::tblastn_report_mask_intervals;
use super::stage_d_pipeline::LocalStageDScoring;
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{
    aligned_sequences_into, record_coordinates, ReportRanges, ShownSubjects,
};
use super::x_subject_store::SubjectLetters;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ProteinScoringSpec;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names::{self, QueryText};
use crate::report::blast_xml::subjects::SubjectNames;
use crate::report::blast_xml::{
    filter, QueryMask, StructuredXml, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
    REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::{Xml2Run, Xml2Target};
use crate::report::query_warnings::{replaced_o_message, QueryWarnings, INVALID_QUERY_MESSAGE};
use crate::report::serial::seqalign::seqalign_evalue;
use crate::stats::tables::lookup_protein_params;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;

const PROGRAM: &str = "TBLASTN";

/// What the document says of the run besides the hit lists (`Parameters`), made once per run.
#[derive(Clone, Copy, Debug)]
pub(super) struct XmlValues {
    /// `Parameters_expect`: the e-value threshold.
    pub expect: f64,
    /// `-soft_masking`: the filter string ends with `m;`.
    pub soft_masking: bool,
    /// Whether the search is gapped (`Parameters_gap-open` and `-gap-extend` of XML2).
    pub gapped: bool,
    /// `-comp_based_stats` as NCBI's number (0 or 2; `Parameters_cbs` of XML2 when non-zero).
    pub cbs: i32,
    /// `-db_gencode` (`Parameters_db-gencode` of XML2).
    pub db_gencode: i32,
}

/// The document's texts, made (and refused) before the first byte of any format is written:
/// the definition line and message of every query, the names of the subjects with hits, the
/// statistics of every query.
pub(super) struct PreparedXml {
    values: XmlValues,
    /// `None` for a matrix name the formatting table does not know (with a query that is an
    /// error; `prepare`): only a run without a query reaches `write` with it.
    positives: Option<FormattingMatrix>,
    filter: Vec<u8>,
    matrix_name: String,
    gap_open: i32,
    gap_extend: i32,
    db_name: String,
    db_num: i32,
    db_len: i64,
    texts: Vec<QueryText>,
    stats: Vec<XmlStats>,
    /// By subject index: the names of the subjects with hits.
    subject_names: Vec<Option<SubjectNames>>,
    /// By subject index: `Hit_len`, the length of the subject's record.
    subject_lengths: Vec<i32>,
    /// By query: whether NCBI searched the query's batch (it has a valid context in it).
    /// Only XML2 shows it.
    searched: Vec<bool>,
    /// XML2 and JSON: the masks of every query are made, not those of the queries with hits
    /// alone (a query without hits lists its masks too).
    all_masks: bool,
    /// What XML2 needs of the run beyond outfmt 5: the target, the gap costs, `cbs` and
    /// `db-gencode`.
    xml2: Xml2Run,
    /// NCBI's messages about the taxonomy database (XML2 and JSON with `-db`).
    taxa_notes: crate::algorithm::common::database::TaxdbNotes,
}

/// NCBI's error for a matrix that its formatting table does not know, thrown by the report
/// data's constructor when a query is formatted.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
/// ```c++
///     else {
///         string prog_name = Blast_ProgramNameFromType(
///                                            m_Options.GetProgramType());
///         if (prog_name != "blastn" && prog_name != "megablast") {
///             NCBI_THROW(blast::CBlastException, eInvalidArgument,
///                         "unsupported score matrix");
/// ```
/// (the matrix table is made from the name as typed, as BLASTP's.)
pub(super) fn unsupported_matrix(matrix_name: &str, has_queries: bool) -> Result<()> {
    if FormattingMatrix::new(matrix_name).is_none() && has_queries {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: "BLAST engine error: unsupported score matrix\n".to_string(),
        }
        .into());
    }
    Ok(())
}

/// The messages of a query's result (the warnings the other formats write to standard error,
/// `report::query_warning`): the O characters read as X, a record without letters, a query
/// without a valid context.
fn query_messages(query: &FastaRecord, valid: bool) -> Vec<String> {
    let mut messages = Vec::new();
    messages.extend(replaced_o_message(query.seq()));
    if query.seq().is_empty() {
        messages.push(crate::algorithm::blastn::blast_engine::NO_DATA_MESSAGE.to_string());
    } else if !valid {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    messages
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    values: XmlValues,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    hitlists: &[KappaResultHitList],
    parameters: &LocalSubjectParameters,
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    seg: Option<&SegParams>,
    ranges: &ReportRanges<'_>,
    xml2: bool,
    taxa: XmlTaxa,
) -> Result<PreparedXml> {
    ensure!(
        query_validity.len() == query_records.len()
            && parameters.lengths.len() >= query_records.len(),
        "TBLASTN query statistics count mismatch"
    );
    unsupported_matrix(matrix_name, !query_records.is_empty())?;
    let positives = FormattingMatrix::new(matrix_name);
    let spec = ProteinScoringSpec {
        matrix: scoring.matrix,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
    };
    let gapped = lookup_protein_params(&spec);
    let stats = (0..query_records.len())
        .map(|q_idx| {
            if !query_validity[q_idx] {
                return XmlStats::NONE;
            }
            let length = &parameters.lengths[q_idx];
            XmlStats {
                hsp_len: length.length_adjustment as i32,
                eff_space: length.eff_searchsp,
                kappa: gapped.k,
                lambda: gapped.lambda,
                entropy: gapped.h,
            }
        })
        .collect();
    let has_hits = |hitlist: &KappaResultHitList| {
        hitlist
            .lists()
            .iter()
            .any(|list| !list.hsps.hsps.is_empty())
    };
    let texts = query_records
        .iter()
        .zip(hitlists)
        .zip(query_validity)
        .map(|((query, hitlist), &valid)| {
            names::query_text(
                PROGRAM,
                &query.local_id,
                query.title(),
                true,
                &query_messages(query, valid),
                has_hits(hitlist),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    // The names of the subjects with hits, made once each.
    let mut subject_names: Vec<Option<SubjectNames>> = Vec::new();
    let count = match shown {
        ShownSubjects::Fasta { .. } => subject_records.len(),
        ShownSubjects::Database { reported, .. } => reported.records.len(),
    };
    subject_names.resize_with(count, || None);
    let mut subject_lengths = vec![0i32; count];
    for list in hitlists.iter().flat_map(KappaResultHitList::lists) {
        if list.hsps.hsps.is_empty() {
            continue;
        }
        let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
        if subject_names[oid].is_some() {
            continue;
        }
        let record = subject_records
            .get(oid)
            .context("TBLASTN subject OID out of range")?;
        subject_names[oid] = Some(match shown.ordinal(oid) {
            None => SubjectNames::of_record(record.local_id.as_bytes(), record.title(), false)
                .map_err(|refusal| refusal.error(PROGRAM, "a subject title"))?,
            // The taxid is shown by XML2 alone (`taxa` is empty for the other formats).
            Some(ordinal) => match shown.parsed_subject(oid) {
                Some(parsed) => SubjectNames::of_parsed(
                    parsed,
                    taxa.parsed_taxonomy.get(oid).map_or(&[], Vec::as_slice),
                ),
                None => SubjectNames::of_database_sequence(
                    ordinal,
                    record.title(),
                    taxa.taxids.get(oid).copied().unwrap_or(0),
                    taxa.taxonomy.get(oid).cloned().unwrap_or_default(),
                ),
            }
            .map_err(|refusal| refusal.error(PROGRAM, "a database title"))?,
        });
        // The whole record (`-subject_loc` searches part of it).
        subject_lengths[oid] =
            i32::try_from(ranges.subjects.length(oid, subject_letters.length(oid)))?;
    }
    let (db_name, db_num, db_len) = match shown {
        ShownSubjects::Fasta { .. } => (String::new(), 0, 0),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.to_string(),
            i32::try_from(pairwise.num_sequences)?,
            i64::try_from(pairwise.total_letters)?,
        ),
    };
    let xml2_run = Xml2Run {
        target: match shown {
            ShownSubjects::Fasta { .. } => Xml2Target::Subjects(subject_records.len()),
            ShownSubjects::Database { .. } => Xml2Target::Database,
        },
        gapped: values.gapped,
        cbs: values.cbs,
        query_gencode: None,
        db_gencode: Some(values.db_gencode),
    };
    Ok(PreparedXml {
        values,
        positives,
        filter: filter::filter_string(
            None,
            seg.map(|seg| (seg.window as i32, seg.locut, seg.hicut)),
            None,
            values.soft_masking,
        ),
        matrix_name: matrix_name.to_string(),
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        db_name,
        db_num,
        db_len,
        texts,
        stats,
        subject_names,
        subject_lengths,
        searched: query_batch_skipped
            .iter()
            .map(|&skipped| !skipped)
            .collect(),
        all_masks: xml2,
        xml2: xml2_run,
        taxa_notes: taxa.notes,
    })
}

/// The masked regions of a query in NCBI's form (`Blast_GetSeqLocInfoVector`): the SEG (and
/// lower-case) intervals of the searched letters (half-open, not yet joined), joined where
/// they overlap (`BlastSeqLocCombine(filter_out, 0)`, blast_filter.c:1256; adjacent intervals
/// stay apart), other than an empty one and the whole searched interval, as inclusive
/// positions in the record (`offset` is the interval's start), frame 0.
fn query_masks(
    intervals: &[(usize, usize)],
    length: usize,
    offset: usize,
    out: &mut Vec<QueryMask>,
) {
    out.clear();
    let ranges = intervals
        .iter()
        .map(|&(start, end)| (start as i32, end as i32 - 1))
        .collect();
    for (from, to) in crate::algorithm::blastx::query_setup::combine_masks(ranges) {
        if from > to || (from == 0 && to == length as i32 - 1) {
            continue;
        }
        out.push(QueryMask {
            from: from + offset as i32,
            to: to + offset as i32,
            frame: 0,
        });
    }
}

impl PreparedXml {
    /// What the XML2 and JSON documents need of the run beyond outfmt 5's (`Xml2Writer::new`).
    pub(super) fn xml2_run(&self) -> Xml2Run {
        self.xml2
    }

    /// Writes the document of every query of `hitlists` (indexed by query) through `document`
    /// (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before
    /// it, with the document's header before the first query and its end after the last when
    /// `ranges.epilog`.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn write<W: StructuredXml>(
        &self,
        document: &mut W,
        writer: &mut dyn Write,
        query_records: &[FastaRecord],
        subject_letters: SubjectLetters<'_>,
        hitlists: &[KappaResultHitList],
        genetic_code: u8,
        seg: Option<&SegParams>,
        mask_lowercase: bool,
        ranges: &ReportRanges<'_>,
        mut warnings: Option<&mut QueryWarnings<'_>>,
    ) -> Result<()> {
        let Some(positives) = &self.positives else {
            // A run without a query (with one, `prepare` refuses) and a matrix the formatting
            // table does not know (`-matrix blosum62`: TBLASTN runs it, the table does not):
            // NCBI never builds the table and writes the epilog alone (for outfmt 5 its end of
            // line).
            if ranges.epilog {
                let mut buf = Vec::new();
                document.finish(&mut buf)?;
                writer.write_all(&buf)?;
            }
            return Ok(());
        };
        let run = XmlRun {
            program: "tblastn",
            version: "TBLASTN 2.17.0+",
            reference: REFERENCE_GAPPED_BLAST,
            db_name: self.db_name.as_bytes(),
            matrix: Some(&self.matrix_name),
            expect: self.values.expect,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            reward_penalty: None,
            filter: &self.filter,
            db_num: self.db_num,
            db_len: self.db_len,
            positives: Some(positives),
            mask_char: b'X',
        };
        let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
        let mut taxdb_printer = self.taxa_notes.printer();
        let mut buf: Vec<u8> = Vec::new();
        let mut masks: Vec<QueryMask> = Vec::new();
        let mut x_window: Vec<u8> = Vec::new();
        let mut qrow = String::new();
        let mut srow = String::new();
        for (q_idx, (query, hitlist)) in query_records.iter().zip(hitlists).enumerate() {
            // The query's warnings come before its document (`QueryWarnings`).
            if let Some(warnings) = warnings.as_deref_mut() {
                warnings.before_query(q_idx, &mut *writer)?;
            }
            let text = &self.texts[q_idx];
            let query_letters = query.seq();
            masks.clear();
            // A query of a batch NCBI did not search has no masks (its result is made
            // without the search, probe `BATCH_SIZE=10` with `-lcase_masking`); XML2 and JSON
            // list the masks of every other query, also one without hits.
            if self.searched[q_idx]
                && (self.all_masks
                    || hitlist
                        .lists()
                        .iter()
                        .any(|list| !list.hsps.hsps.is_empty()))
            {
                // The masks of the report: the SEG and lower-case intervals, joined as NCBI
                // joins them (`query_masks`).
                let intervals = tblastn_report_mask_intervals(query_letters, seg, mask_lowercase);
                query_masks(
                    &intervals,
                    query_letters.len(),
                    ranges.queries.offset(q_idx),
                    &mut masks,
                );
            }
            let xml_query = XmlQuery {
                id: query.local_id.as_bytes(),
                title: &text.title,
                len: i32::try_from(query_letters.len())?,
                stats: self.stats[q_idx],
                message: text.message.as_deref(),
                // Read by XML2 alone.
                searched: self.searched[q_idx],
                masks: &masks,
            };
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                for list in hitlist.lists() {
                    if list.hsps.hsps.is_empty() {
                        continue;
                    }
                    let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
                    ensure!(
                        list.payloads.len() == list.hsps.hsps.len(),
                        "TBLASTN report payload count mismatch"
                    );
                    let names = self.subject_names[oid]
                        .as_ref()
                        .context("names of a subject with hits")?;
                    let record_length = self.subject_lengths[oid];
                    let hit = XmlHit::of(names, record_length);
                    let subject_seq = subject_letters.nt(oid);
                    let mut translation = TargetTranslation::of(subject_seq, &code);
                    // NCBI's messages about the taxonomy database: opening and the volume's
                    // lookup at the first `Hit` that needs them, the lookup of the `Hit` every
                    // time (nothing for the other formats).
                    if let Some(warnings) = warnings.as_deref_mut() {
                        taxdb_printer.at_hit(
                            oid,
                            names.taxid,
                            &names.taxonomy,
                            &mut *warnings.sink,
                        )?;
                    }
                    sink.hit(&hit, |hsps| {
                        for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                            let hsp = &linked.hsp;
                            let frame = i32::from(hsp.frame);
                            // NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105:
                            // translated coordinates use three nucleotide letters per
                            // residue and reverse strand endpoints run in descending order.
                            let subject_len = i32::try_from(subject_seq.len())?;
                            let (s_start, s_end) = if frame > 0 {
                                (3 * hsp.s_start + frame, 3 * hsp.s_end + frame - 1)
                            } else {
                                (
                                    subject_len - 3 * hsp.s_start + frame + 1,
                                    subject_len - 3 * hsp.s_end + frame + 2,
                                )
                            };
                            // The rows, without the lower case of the masks.
                            aligned_sequences_into(
                                query_letters,
                                subject_seq,
                                &mut x_window,
                                &code,
                                &mut translation,
                                hsp.frame,
                                hsp.q_start,
                                hsp.s_start,
                                hsp.s_end,
                                &payload.edit_script,
                                &[],
                                &mut qrow,
                                &mut srow,
                            )?;
                            // The coordinates in the searched letters, moved into the records.
                            let (q_start, q_end, s_start, s_end) = record_coordinates(
                                ranges,
                                q_idx,
                                oid,
                                (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                            )?;
                            // The Seq-align's endpoints, lower first: the minus frame's
                            // outfmt 6 order is not the XML's.
                            let (hit_from, hit_to) = (s_start.min(s_end), s_start.max(s_end));
                            hsps.hsp(&XmlHsp {
                                score: hsp.score,
                                bit_score: payload.bit_score,
                                evalue: seqalign_evalue(linked.evalue),
                                query_from: q_start,
                                query_to: q_end,
                                hit_from,
                                hit_to,
                                query_frame: 0,
                                hit_frame: translation_frame(
                                    frame > 0,
                                    hit_from,
                                    hit_to,
                                    record_length,
                                ),
                                qrow: qrow.as_bytes(),
                                srow: srow.as_bytes(),
                                ops: Some(&payload.edit_script),
                                query_first: q_start - 1,
                                query_step: 1,
                            })?;
                        }
                        Ok(())
                    })?;
                }
                Ok(())
            });
            // A failed query leaves its partial bytes in `buf`, which the error return drops.
            written?;
            writer.write_all(&buf)?;
            buf.clear();
        }
        if ranges.epilog {
            document.finish(&mut buf)?;
            writer.write_all(&buf)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn masks_of(intervals: &[(usize, usize)], length: usize, offset: usize) -> Vec<(i32, i32)> {
        let mut out = Vec::new();
        query_masks(intervals, length, offset, &mut out);
        out.iter()
            .map(|mask| {
                assert_eq!(mask.frame, 0);
                (mask.from, mask.to)
            })
            .collect()
    }

    #[test]
    fn adjacent_intervals_stay_two_masks_and_overlapping_ones_join() {
        // A lower-case run [96, 121) and a SEG interval [121, 153) of a `-query_loc 21-195`
        // query (fuzz seed 2, case_200159): NCBI keeps both (blast_filter.c:993), and its XML
        // mask loop then leaves the letter at 141 (0-based) unmasked when 140 is the
        // second-last column of a segment.
        assert_eq!(
            masks_of(&[(121, 153), (96, 121)], 175, 20),
            vec![(116, 140), (141, 172)]
        );
        // One shared position joins them (inclusive 120 > 120 is false: still apart).
        assert_eq!(
            masks_of(&[(96, 121), (120, 153)], 175, 20),
            vec![(116, 140), (140, 172)]
        );
        // Two shared positions join them.
        assert_eq!(
            masks_of(&[(96, 122), (120, 153)], 175, 20),
            vec![(116, 172)]
        );
        // The whole searched interval is not a mask.
        assert_eq!(
            masks_of(&[(0, 100), (100, 175)], 175, 20),
            vec![(20, 119), (120, 194)]
        );
        assert_eq!(masks_of(&[(0, 175)], 175, 20), Vec::<(i32, i32)>::new());
    }

    #[test]
    fn report_intervals_keep_the_lowercase_run_apart_from_seg() {
        // Without SEG, the lower-case runs alone, unjoined.
        let query = b"ACDEFghiklMNPQrsT";
        assert_eq!(
            super::super::search_seed::tblastn_report_mask_intervals(query, None, true),
            vec![(5, 10), (14, 16)]
        );
        assert!(
            super::super::search_seed::tblastn_report_mask_intervals(query, None, false).is_empty()
        );
    }
}
