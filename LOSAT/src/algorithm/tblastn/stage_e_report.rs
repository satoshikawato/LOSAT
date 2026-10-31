//! NCBI TBLASTN local-subject result formatting.

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, ensure, Context, Result};

use super::search_gapped::TargetTranslation;
use super::search_seed::encode_tblastn_lookup_query;
use super::stage_d_pipeline::{LocalStageDScoring, SubjectSetNumbers};
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::x_subject_store::SubjectLetters;
use crate::api::local_blast::{FormatProbe, HspIndex, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::common::{GapEditOp, Hit};
use crate::config::ProteinScoringSpec;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
};
use crate::report::pairwise::{
    write_tblastn_pairwise_report, BlastpPairwiseQuery, BlastpPairwiseReport, PairwiseConfig,
    PairwiseDatabase, PairwiseHit,
};
use crate::sequence::packed_subject::SubjectNt;
use crate::stats::spouge::lookup_protein_gumbel_params;
use crate::stats::tables::lookup_protein_params;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
// ```c++
//     info.name = seqdb->GetDBNameList();
//     info.definition = seqdb->GetTitle();
//     if (info.definition.empty())
//         info.definition = info.name;
// ```
/// The subject set as the reports name it.
pub(super) enum ShownSubjects<'a> {
    /// `-subject`: `User specified sequence set (Input: PATH)`; the records are the whole set.
    Fasta { path: &'a Path },
    /// `-db`: the `-db` value (outfmt 7's `# Database:` line), the database's `Database:`
    /// lines (outfmt 0) and the reported sequences by rank (their tabular ids).
    Database {
        name: &'a str,
        pairwise: &'a PairwiseDatabase,
        reported: &'a super::database::ReportedSubjects,
    },
}

impl ShownSubjects<'_> {
    /// The name after `# Database: ` (outfmt 7) and in the pairwise report of a subject set.
    pub(super) fn name(&self) -> String {
        match self {
            Self::Fasta { path } => {
                format!("User specified sequence set (Input: {})", path.display())
            }
            Self::Database { name, .. } => name.to_string(),
        }
    }

    /// The `Database:` lines of a database search's pairwise report.
    fn pairwise_database(&self) -> Option<&PairwiseDatabase> {
        match self {
            Self::Fasta { .. } => None,
            Self::Database { pairwise, .. } => Some(pairwise),
        }
    }

    /// The database title of the Seq-annots (`RunRecord::db_title`; NCBI's
    /// `m_DbInfo.front().definition`, the `Database:` line's title), `None` for `-subject`.
    pub(super) fn db_title(&self) -> Option<&[u8]> {
        self.pairwise_database()
            .map(|pairwise| pairwise.definition.as_slice())
    }

    /// The ordinal ID of the reported database sequence of rank `rank`, `None` for `-subject`.
    pub(super) fn ordinal(&self, rank: usize) -> Option<u32> {
        match self {
            Self::Fasta { .. } => None,
            Self::Database { reported, .. } => Some(reported.ordinals[rank]),
        }
    }

    /// The reported database sequences with parsed Seq-ids by rank (empty without them).
    fn parsed_subjects(&self) -> &[Option<Arc<crate::report::db_subject::ParsedSubject>>] {
        match self {
            Self::Fasta { .. } => &[],
            Self::Database { reported, .. } => &reported.parsed,
        }
    }
}

impl<'a> ShownSubjects<'a> {
    /// The reported database sequence of rank `rank` when it has parsed Seq-ids.
    pub(super) fn parsed_subject(
        &self,
        rank: usize,
    ) -> Option<&'a crate::report::db_subject::ParsedSubject> {
        match *self {
            Self::Fasta { .. } => None,
            Self::Database { reported, .. } => reported.parsed_subject(rank),
        }
    }
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:759-835,1411-1458:
// if (m_FormatType == eTabular || m_FormatType == eTabularWithComments) {
//     x_PrintTabularReport(results, itr_num); return;
// }
// The formatter receives already sorted query hitlists from Stage D.
// NCBI c++/src/app/blast/blast_formatter.cpp:429-467:
// ITERATE(CSearchResultSet, result, *results) {
//     ...
//         formatter.PrintOneResultSet(**result, queries);
//     ...
// }
// Each requested format prints the same final result without searching again.
/// The number of the BLAST archive among the formats (outfmt 11).
pub(super) const ARCHIVE_FORMAT: i32 = 11;

/// The format number of a checked `-outfmt` value, as NCBI's `ParseFormattingString` reads
/// it (a custom specification of outfmt 0 is ignored); -1 for a value it does not read.
pub(super) fn format_number(outfmt: &str) -> i32 {
    crate::blastinput::app::parse_formatting_string(outfmt).map_or(-1, |choice| choice.number)
}

/// The XML2 / JSON document of a checked `-outfmt` value (16, 15, 14 or 13; `None` for the
/// others). Outfmt 13 and 14 write side files next to `-out`: only the command line reaches them
/// (`run_local` stops at NCBI's XInclude check).
pub(super) fn xml2_kind(outfmt: &str) -> Option<crate::report::blast_xml2::Xml2Kind> {
    use crate::blastinput::app::ReportFormat;
    crate::report::blast_xml2::Xml2Kind::of(ReportFormat::from_number(format_number(outfmt))?)
}

/// Whether a checked `-outfmt` value is a format of the structured XML family that takes the
/// path of outfmt 5: BLAST XML (5) and the XML2 and JSON documents (16, 15, 14, 13).
pub(super) fn is_structured_xml(outfmt: &str) -> bool {
    format_number(outfmt) == 5 || xml2_kind(outfmt).is_some()
}

/// Whether the reports `outfmts` (with `hits`, the caller's hit records) read the letters of the
/// subjects: outfmt 0 and the hit records (`pairwise_hits`), the XML family
/// (`blast_xml::XmlReport::write`) and the tabular rows (`tabular::reads_rows`). The Seq-aligns
/// and the archive (8, 9, 11, 12) and the other tabular reports (`write_tabular`, the default
/// rows of 6 and 7, has no spec) read the lengths alone: their counts are the search's. A value
/// that this does not know reads them.
pub(super) fn reads_subject_letters<'a>(
    hits: bool,
    mut outfmts: impl Iterator<Item = &'a str>,
) -> bool {
    hits || outfmts.any(|outfmt| match format_number(outfmt) {
        8 | 9 | 11 | 12 => false,
        6 | 7 | 10 | 20 => super::tabular::tabular_spec(outfmt)
            .map_or(false, |spec| super::tabular::reads_rows(&spec)),
        _ => true,
    })
}

/// The ranges of a search as the reports use them: where the searched letters of each
/// query and subject lie in its record (`-query_loc`, `-subject_loc`), and whether the
/// reports end with the epilog (a search that stops at a query batch writes none).
pub(super) struct ReportRanges<'a> {
    pub queries: &'a crate::blastinput::seq_range::Placements,
    pub subjects: &'a crate::blastinput::seq_range::Placements,
    pub epilog: bool,
}

#[allow(clippy::too_many_arguments)]
pub(super) fn render(
    outputs: &mut ReportOutputs<'_>,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    subject_set: SubjectSetNumbers,
    hitlists: &mut [KappaResultHitList],
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    report_counts: crate::blastinput::report_counts::ReportCounts,
    genetic_code: u8,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    query_warning_lines: &[Vec<u8>],
    ranges: &ReportRanges<'_>,
    xml_values: super::blast_xml::XmlValues,
    xml_taxa: super::database::XmlTaxa,
    out: Option<&std::path::Path>,
    archive: Option<&super::archive::ArchiveRender<'_>>,
    tabular: &super::tabular::TabularInputs,
) -> Result<()> {
    ensure!(
        query_records.len() == hitlists.len(),
        "TBLASTN report query count mismatch"
    );
    // NCBI c++/src/algo/blast/api/blast_seqalign.cpp:1572-1577:
    // Blast_HSPListSortByEvalue(hsp_list); /* during Seq-align conversion */
    for hitlist in hitlists.iter_mut() {
        hitlist.sort_hsps_for_report();
        hitlist.drop_zero_score_hsps();
    }
    let hitlists: &[KappaResultHitList] = hitlists;
    // BLAST XML (outfmt 5, and the XML2 and JSON documents of 15 and 16): its texts are made,
    // and what NCBI cannot write is refused, before any format writes (`-out` is already
    // created, empty, as NCBI's).
    //
    // NCBI writes the XML2 and JSON documents' start before the search (`PrintProlog`) and
    // fails for a matrix of its formatting table's that it does not know when it makes the
    // first report's data: the start alone is written then (probe `tblastn_sub_lcmatrix`).
    if outputs
        .formats
        .iter()
        .any(|format| is_structured_xml(format.outfmt))
    {
        if let Err(error) =
            super::blast_xml::unsupported_matrix(matrix_name, !query_records.is_empty())
        {
            for format in outputs.formats.iter_mut() {
                if let Some(kind) = xml2_kind(format.outfmt) {
                    let mut writer = format.sink.open()?;
                    writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                    writer.flush()?;
                }
            }
            return Err(error);
        }
    }
    let xml_is_xml2 = outputs
        .formats
        .iter()
        .any(|format| xml2_kind(format.outfmt).is_some());
    let xml = if outputs
        .formats
        .iter()
        .any(|format| is_structured_xml(format.outfmt))
    {
        Some(super::blast_xml::prepare(
            xml_values,
            query_records,
            subject_records,
            subject_letters,
            shown,
            hitlists,
            parameters,
            query_validity,
            query_batch_skipped,
            scoring,
            matrix_name,
            seg,
            ranges,
            xml_is_xml2,
            xml_taxa,
        )?)
    } else {
        None
    };
    // The tabular subject IDs, made before any format is written: a subject that NCBI
    // writes part of a row for and then stops at (RP-20) is rejected first
    // (`tabular_subject_ids`). The tabular formats 6, 7, 10 and 20 show them.
    let subject_ids = if outputs
        .formats
        .iter()
        .any(|format| matches!(format_number(format.outfmt), 6 | 7 | 10 | 20))
    {
        match shown {
            ShownSubjects::Fasta { .. } => tabular_subject_ids(hitlists, subject_records)?,
            ShownSubjects::Database { reported, .. } => reported.tabular_subject_ids()?,
        }
    } else {
        Vec::new()
    };
    // NCBI c++/src/algo/blast/format/blast_format.cpp:1411:
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // The alignments are rendered once from the final result, for outfmt 0 and for
    // the caller's hit records, in the order of the tabular rows.
    let pairwise_hits = if outputs.hits.is_some()
        || outputs
            .formats
            .iter()
            .any(|format| format_number(format.outfmt) == 0)
    {
        Some(pairwise_hits(
            query_records,
            subject_records,
            subject_letters,
            hitlists,
            scoring,
            genetic_code,
            seg,
            mask_lowercase,
            ranges,
        )?)
    } else {
        None
    };
    if let (Some(hits_sink), Some(pairwise_hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref())
    {
        hits_sink(pairwise_hits);
    }
    // The warnings of each query are written before its report, once, with the first
    // format (`QueryWarnings`); the taxonomy database's messages of XML2 and JSON go to the
    // diagnostics at every format.
    let mut warning_lines = Some(query_warning_lines);
    let observer = &mut outputs.observer;
    for (format_index, format) in outputs.formats.iter_mut().enumerate() {
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a call with several formats, `run_local`).
        if let Some(archive) = archive {
            archive
                .handle
                .set_muted(format_number(format.outfmt) != ARCHIVE_FORMAT);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // A tabular report with a field list or a delimiter, or outfmt 10 or 20
        // (`tabular::tabular_spec`); `None` for the default 6 and 7.
        let custom_tabular = super::tabular::tabular_spec(format.outfmt);
        // So do the warnings of the formatters of a tabular field list that shows a taxonomy
        // name (`x_CheckTaxDB`, `tabular`).
        let mut warnings = (warning_lines.is_some()
            || xml2_kind(format.outfmt).is_some()
            || format_number(format.outfmt) == ARCHIVE_FORMAT
            || custom_tabular
                .as_ref()
                .is_some_and(crate::report::tabular::TabularSpec::checks_taxdb))
        .then(|| crate::report::query_warnings::QueryWarnings {
            before: warning_lines.take().unwrap_or(&[]),
            sink: &mut *outputs.diagnostics,
        });
        let mut writer = format.sink.open()?;
        match format_number(format.outfmt) {
            // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
            // (`tabular`); never the default 6 and 7 below.
            6 | 7 | 10 | 20 if custom_tabular.is_some() => super::tabular::write_tblastn_tabular(
                &super::tabular::TblastnTabular {
                    spec: custom_tabular.as_ref().expect("a custom tabular report"),
                    inputs: tabular,
                    queries: query_records,
                    letters: subject_letters,
                    subject_ids: &subject_ids,
                    shown,
                    query_batch_skipped,
                    ranges,
                    genetic_code,
                    matrix: scoring.matrix,
                },
                hitlists,
                &mut writer,
                warnings.as_mut(),
                probe.as_mut(),
            )?,
            6 | 7 => write_tabular(
                &mut writer,
                format_number(format.outfmt) == 7,
                query_records,
                subject_records,
                subject_letters,
                &subject_ids,
                shown,
                hitlists,
                query_batch_skipped,
                ranges,
                probe.as_mut(),
                warnings.as_mut(),
            )?,
            0 => write_pairwise(
                &mut writer,
                pairwise_hits
                    .as_deref()
                    .expect("pairwise hits prepared for TBLASTN outfmt 0"),
                query_records,
                subject_records,
                subject_letters,
                shown,
                subject_set,
                parameters,
                ungapped_karlin,
                query_validity,
                query_batch_skipped,
                scoring,
                matrix_name,
                report_counts,
                ranges,
                probe.as_mut(),
                warnings.as_mut(),
            )?,
            // BLAST XML: one document per run, a header before the first query and a
            // trailer after the last (`blast_xml`).
            5 => xml
                .as_ref()
                .expect("XML prepared for TBLASTN outfmt 5")
                .write(
                    &mut crate::report::blast_xml::BlastXmlWriter::new(0),
                    &mut writer,
                    query_records,
                    subject_letters,
                    hitlists,
                    genetic_code,
                    seg,
                    mask_lowercase,
                    ranges,
                    warnings.as_mut(),
                )?,
            // The XML2 and JSON documents (outfmt 16 and 15, and 14 and 13 as side files next to
            // `out` with this output as the master file): the same records, written by
            // `Xml2Writer`.
            13..=16 => {
                let kind = xml2_kind(format.outfmt).expect("a checked XML2 or JSON format");
                let xml = xml
                    .as_ref()
                    .expect("XML prepared for TBLASTN outfmt 13 to 16");
                xml.write(
                    &mut crate::report::blast_xml2::Xml2Writer::for_run(kind, xml.xml2_run(), out),
                    &mut writer,
                    query_records,
                    subject_letters,
                    hitlists,
                    genetic_code,
                    seg,
                    mask_lowercase,
                    ranges,
                    warnings.as_mut(),
                )?;
                // A run that stops before its first report (a database whose sequences cannot
                // be read) still writes the document's start, which NCBI writes before the
                // search.
                if query_records.is_empty() && !ranges.epilog {
                    writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                }
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
            //         if (results.HasAlignments()) {
            //     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
            //             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
            //             if(m_FormatType == CFormattingArgs::eAsnText)
            //             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //             else
            //             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //         }
            // ```
            // (`eAsnBinary` is the same pruning and `MSerial_AsnBinary`,
            // blast_format.cpp:699-705.) The Seq-annot of each query with alignments.
            number @ (8 | 9 | 12) => super::seqalign::write_tblastn_seqalign(
                match number {
                    8 => crate::report::serial::SerialFormat::AsnText,
                    9 => crate::report::serial::SerialFormat::AsnBinary,
                    _ => crate::report::serial::SerialFormat::Json,
                },
                &mut writer,
                query_records,
                subject_records,
                subject_letters,
                shown,
                hitlists,
                ranges,
                warnings.as_mut(),
            )?,
            // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:295
            // ```c++
            //                 if (fmt_args->ArchiveFormatRequested(args)) {
            //                     formatter.WriteArchive(*query_factory, *m_OptsHndl, *results, 0, m_Bah.GetMessages());
            //                     m_Bah.ResetMessages();
            // ```
            // One BLAST archive per query batch (`archive`).
            ARCHIVE_FORMAT => {
                let archive = archive.expect("a BLAST archive request has its handle");
                let gapped = lookup_protein_params(&ProteinScoringSpec {
                    matrix: scoring.matrix,
                    gap_open: scoring.gap_open,
                    gap_extend: scoring.gap_extend,
                });
                super::archive::write_tblastn_archives(
                    archive,
                    &super::archive::ArchiveData {
                        queries: query_records,
                        subject_records,
                        subject_letters,
                        shown,
                        hitlists,
                        parameters,
                        ungapped_karlin,
                        query_validity,
                        query_batch_skipped,
                        gapped,
                        seg,
                        mask_lowercase,
                        ranges,
                    },
                    warnings.as_mut(),
                    &mut writer,
                )?
            }
            _ => bail!("unsupported TBLASTN outfmt {}", format.outfmt),
        }
        writer.flush()?;
    }
    // The messages that `search` writes after the formats (the stop messages) are the archive's
    // again, whichever format came last (as in blastn and blastp).
    if let Some(archive) = archive {
        archive.handle.set_muted(false);
    }
    Ok(())
}

// NCBI c++/src/objtools/align_format/tabular.cpp:1266-1338:
// x_PrintQueryAndDbNames(...);
// if (align_set) { if (num_hits != 0) PrintFieldNames(...);
//                 m_Ostream << "# " << num_hits << " hits found\n"; }
// CBlastTabularInfo::PrintNumProcessed appends the final query count.
#[allow(clippy::too_many_arguments)]
fn write_tabular(
    writer: &mut impl Write,
    comments: bool,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    subject_ids: &[Arc<[u8]>],
    shown: &ShownSubjects<'_>,
    hitlists: &[KappaResultHitList],
    query_batch_skipped: &[bool],
    ranges: &ReportRanges<'_>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> Result<()> {
    // NCBI c++/src/algo/blast/api/local_blast.cpp:177-224:
    // a skipped batch has a null align set for each query it contains.
    ensure!(
        query_batch_skipped.len() == query_records.len(),
        "TBLASTN batch validity count mismatch"
    );
    // The rows are printed in the order of the final HSP list, so a running count is
    // each row's HSP index.
    let mut hsp_index: HspIndex = 0;
    for (q_idx, ((query, hitlist), &search_skipped)) in query_records
        .iter()
        .zip(hitlists)
        .zip(query_batch_skipped)
        .enumerate()
    {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let hit_count: usize = hitlist
            .lists()
            .iter()
            .map(|list| list.hsps.hsps.len())
            .sum();
        if comments {
            writeln!(writer, "# TBLASTN 2.17.0+")?;
            // NCBI reference: c++/src/objtools/align_format/tabular.cpp:1305-1308
            // ```c
            //     CAlignFormatUtil::AcknowledgeBlastQuery(bioseq, kLineLength, m_Ostream,
            //                                             m_ParseLocalIds, kHtmlFormat,
            //                                             kTabularFormat, rid);
            // ```
            // The title's bytes (`write_outfmt7_query_line`).
            crate::report::outfmt6::write_outfmt7_query_line(writer, &query.query_line())?;
            writeln!(writer, "# Database: {}", shown.name())?;
            if hit_count > 0 {
                writeln!(writer, "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")?;
            }
            // NCBI c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // if (status != 0) { CRef<CSeq_align_set> tmp_align;
            //                    sa_vec.push_back(tmp_align); }
            // NCBI c++/src/objtools/align_format/tabular.cpp:1266-1338:
            // if (align_set) { m_Ostream << "# " << num_hits << " hits found\n"; }
            // A skipped search has a null alignment set, so the count line
            // is omitted. An ordinary searched no-hit query retains it.
            if !search_skipped {
                writeln!(writer, "# {hit_count} hits found")?;
            }
        }
        for list in hitlist.lists() {
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            let subject = subject_records
                .get(oid)
                .context("TBLASTN subject OID out of range")?;
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                let hsp = &linked.hsp;
                let frame = i32::from(hsp.frame);
                // NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105:
                // if (segment->frame < 0) {
                //     *start = seq_length - 3*segment->offset + segment->frame;
                //     *end = seq_length - 3*segment->end + segment->frame + 1;
                // } else if (segment->frame > 0) {
                //     *start = 3*segment->offset + segment->frame - 1;
                //     *end = 3*segment->end + segment->frame - 2;
                // }
                // Stage D's local subject coordinate conversion is 1-offset.
                let subject_len = i32::try_from(subject_letters.length(oid))?;
                let (s_start, s_end) = if frame > 0 {
                    (3 * hsp.s_start + frame, 3 * hsp.s_end + frame - 1)
                } else {
                    (
                        subject_len - 3 * hsp.s_start + frame + 1,
                        subject_len - 3 * hsp.s_end + frame + 2,
                    )
                };
                // The coordinates in the searched letters, moved into the records
                // (`record_coordinates`).
                let (q_start, q_end, s_start, s_end) = record_coordinates(
                    ranges,
                    q_idx,
                    oid,
                    (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                )?;
                let percent =
                    format_percent_identity_ncbi(payload.report_num_ident, payload.align_length, 3);
                let bits = format_bitscore_ncbi(payload.bit_score);
                let evalue = format_evalue_ncbi_tabular(linked.evalue);
                // NCBI c++/src/objtools/align_format/tabular.cpp:1100-1108:
                // ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
                //     // Add tab in front of field, except for the first field.
                //     if (iter != m_FieldsToShow.begin())
                //         m_Ostream << m_FieldDelimiter;
                //     x_PrintField(*iter);
                // }
                // m_Ostream << "\n";
                // One printed row is one HSP; the probe marks it without changing it.
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-497
                // ```c++
                //         if (sid_in->IsLocal()) {
                //             string id_token;
                //             vector<string> title_tokens;
                //             title_tokens =
                //                 NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
                //             if(title_tokens.empty()){
                //                 id_token = NcbiEmptyString;
                //             } else {
                //                 id_token = title_tokens[0];
                //             }
                //
                //             if (id_token == NcbiEmptyString || parse_local) {
                //                 const CObject_id& obj_id = sid_in->GetLocal();
                //                 if (obj_id.IsStr())
                //                     id_token = obj_id.GetStr();
                // ```
                // The query ID is the first word of its title, or its local ID without a title
                // (`FastaRecord::shown_id`); the subject's is `tabular_subject_ids`'s. Both are
                // written as bytes.
                writer.write_all(query.shown_id())?;
                writer.write_all(b"\t")?;
                writer.write_all(&subject_ids[oid])?;
                writeln!(
                    writer,
                    "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    percent,
                    payload.align_length,
                    payload.report_mismatches,
                    payload.gap_opens,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    evalue,
                    // NCBI c++/src/objtools/align_format/align_format_util.cpp:986-993:
                    // snprintf(bit_score_buf, ..., kBitScoreFormat, bit_score);
                    // kBitScoreFormat is "%4.1lf", retaining width for scores below ten.
                    bits,
                )?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
                hsp_index += 1;
            }
        }
    }
    // NCBI reference: c++/src/app/blast/tblastn_app.cpp:342
    // ```c++
    //         formatter.PrintEpilog(opt);
    // ```
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog.
    if comments && ranges.epilog {
        writeln!(writer, "# BLAST processed {} queries", query_records.len())?;
    }
    Ok(())
}

/// The 1-based reported coordinates of an HSP (query start and end, subject start and end
/// in nucleotides) moved from the searched letters into the records: by the start of the
/// query's interval and of the subject's interval, on either frame. Without ranges nothing
/// moves.
///
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
/// ```c++
/// 	if (seqinfo_src->CanReturnPartialSequence() == true)
/// 	{
///         	CConstRef<CSeq_loc> subj_loc = seqinfo_src->GetSeqLoc(kOid);
///         	NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
///              	   RemapToQueryLoc(*iter, query_loc);
///              	   if ( !is_ooframe )
///                    	s_RemapToSubjectLoc(*iter, *subj_loc);
/// ```
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
/// ```c++
/// void RemapToQueryLoc(CRef<CSeq_align> sar, const CSeq_loc & query)
/// {
///     ...
///     if (query.IsInt()) {
///         q_shift = query.GetInt().GetFrom();
///     }
///
///     if (q_shift > 0) {
///         sar->OffsetRow(query_row, q_shift);
///     }
/// }
/// ```
/// A subject interval has the strand both, so `RemapAlignToLoc` (seq_align_util.cpp:71-78)
/// moves its row by the interval's start.
pub(super) fn record_coordinates(
    ranges: &ReportRanges<'_>,
    q_idx: usize,
    s_idx: usize,
    (q_start, q_end, s_start, s_end): (i32, i32, i32, i32),
) -> Result<(i32, i32, i32, i32)> {
    let q_shift = i32::try_from(ranges.queries.offset(q_idx))?;
    let s_shift = i32::try_from(ranges.subjects.offset(s_idx))?;
    Ok((
        q_start + q_shift,
        q_end + q_shift,
        s_start + s_shift,
        s_end + s_shift,
    ))
}

/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`), as BLASTN's and TBLASTX's.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
/// ```c++
///         } catch (const CException&) {
///             list<CRef<CSeq_id> > subject_ids;
///             CRef<CSeq_id> id(new CSeq_id());
///             id->Assign(align.GetSeq_id(1));
///             subject_ids.push_back(id);
///             SetSubjectId(subject_ids);
///             bioseqs_found = false;
///         }
/// ```
/// A subject with a row whose `GenerateDefline` throws (a `Subject_` first word and a title
/// in no encoding that `CUtf8::GuessEncoding` recognises) makes NCBI 2.17.0 write part of the
/// row and stop with a `CCoreException` that names its build's source files (exit 255;
/// `AUTHORITY.md` §J-6, §K-1, RP-20). No result can be defined for it, so LOSAT rejects such
/// a subject explicitly before it writes any format; a subject without rows is never
/// named.
fn tabular_subject_ids(
    hitlists: &[KappaResultHitList],
    subjects: &[FastaRecord],
) -> Result<Vec<Arc<[u8]>>> {
    let ids: Vec<Option<Arc<[u8]>>> = subjects
        .iter()
        .map(|subject| {
            crate::report::defline::tabular_subject_id(
                subject.local_id.as_bytes(),
                &subject.title,
                false,
            )
            .ok()
            .map(Arc::from)
        })
        .collect();
    for list in hitlists.iter().flat_map(KappaResultHitList::lists) {
        if list.hsps.hsps.is_empty() {
            continue;
        }
        let s_idx = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
        if ids.get(s_idx).is_some_and(Option::is_none) {
            bail!(
                "subject record {} ({}) has a title whose first word starts with 'Subject_' and whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises; NCBI BLAST+ writes part of a tabular row for such a subject and stops with an exception that names its build's source files (exit 255), which is not supported by LOSAT's TBLASTN",
                s_idx + 1,
                String::from_utf8_lossy(subjects[s_idx].shown_id())
            );
        }
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| Arc::from(&b""[..])))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::stage_d_pipeline::{
        run_local_for_report, LocalStageDProfile, LocalStageDScoring,
    };
    use crate::api::local_blast::{FormatOutput, OutputSink};
    use crate::utils::seg::SegParams;

    // The reports that read no letters of the reported database sequences: a value read wrongly
    // here reaches `SubjectLetters::Lengths::nt`, which panics.
    #[test]
    fn only_the_reports_that_show_letters_read_the_subjects() {
        let reads = |outfmt: &str| reads_subject_letters(false, std::iter::once(outfmt));
        for outfmt in [
            "0",
            "5",
            "13",
            "14",
            "15",
            "16",
            "6 sseq",
            "6 qseq",
            "6 btop",
            "6 positive",
            "6 ppos",
            "7 std sseq",
            "10 sseq",
            "20 btop",
            // Alone in the counts' block, `nident` reads the rows of some HSPs.
            "6 nident",
            "not a format",
        ] {
            assert!(reads(outfmt), "{outfmt}");
        }
        for outfmt in [
            "6",
            "7",
            "10",
            "20",
            "6 std",
            "6 qseqid sseqid slen evalue",
            "6 nident pident",
            "8",
            "9",
            "11",
            "12",
        ] {
            assert!(!reads(outfmt), "{outfmt}");
        }
        assert!(reads_subject_letters(true, std::iter::once("6")));
        assert!(!reads_subject_letters(false, std::iter::empty()));
        assert!(reads_subject_letters(false, ["6", "0"].into_iter()));
    }

    // Outfmt 5 and the XML2 and JSON documents (16 and 15) take the structured-XML path of
    // `render`; no other format does.
    #[test]
    fn the_xml_family_is_5_and_13_to_16() {
        use crate::report::blast_xml2::Xml2Kind;
        for (outfmt, kind, structured) in [
            ("5", None, true),
            ("13", Some(Xml2Kind::Files13), true),
            ("14", Some(Xml2Kind::Files14), true),
            ("15", Some(Xml2Kind::Single15), true),
            ("16", Some(Xml2Kind::Single16), true),
            (" 16 std", Some(Xml2Kind::Single16), true),
            ("0", None, false),
            ("6", None, false),
            ("12", None, false),
        ] {
            assert_eq!(xml2_kind(outfmt), kind, "{outfmt}");
            assert_eq!(is_structured_xml(outfmt), structured, "{outfmt}");
        }
    }

    // NCBI c++/src/app/blast/tblastn_app.cpp:288-301:
    // results = lcl_blast.Run(); formatter.PrintOneResultSet(**result, query);
    // Compare complete saved formatter bytes after the Stage D local search.
    // NCBI c++/src/app/blast/blast_formatter.cpp:429-467:
    // ITERATE(CSearchResultSet, result, *results) {
    //     ...
    //         formatter.PrintOneResultSet(**result, queries);
    //     ...
    // }
    // All three formats are printed from one result.
    #[test]
    fn initial_tabular_oracle_bytes() {
        let query_path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/query.faa"
        ));
        let subject_path =
            Path::new("docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna");
        let full_subject_path = Path::new(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../docs/evidence/tlosan_stage_c/multi_hsp_20260924/subjects.fna"
        ));
        // The records of NCBI's reader, as TBLASTN reads the query and the subjects.
        let read = |path: &Path, prefix: &str, protein: bool| {
            use crate::blastinput::fasta_reader::{read_all, FastaInputSource, ReaderConfig};
            let config = match prefix {
                "Query_" => ReaderConfig::query("TBLASTN", protein, false),
                _ => ReaderConfig::subject("TBLASTN", protein, false),
            };
            let file = std::fs::File::open(path).unwrap();
            read_all(&mut FastaInputSource::from_file(file, config), &mut |_| {
                Ok(())
            })
            .unwrap()
        };
        let queries = read(query_path, "Query_", true);
        let subjects = read(full_subject_path, "Subject_", false);
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let subject_seqs: Vec<&[u8]> = subjects.iter().map(|record| record.seq()).collect();
        let seg = SegParams::default();
        for mode in [0, 2] {
            let (mut results, parameters, ungapped_karlin, query_validity) = run_local_for_report(
                &query_seqs,
                &subject_seqs,
                LocalStageDProfile {
                    seg: Some(&seg),
                    soft_masking: false,
                    mask_lowercase: false,
                    genetic_code: 1,
                    expect_value: 10.0,
                    max_target_seqs: 500,
                    query_offsets: &[],
                    subject_set: Default::default(),
                    engine_positives: false,
                    culling_limit: 0,
                },
                mode == 2,
                true,
                LocalStageDScoring::default(),
            )
            .unwrap();
            const FORMATS: [&str; 3] = ["0", "6", "7"];
            let mut outputs_bytes = vec![Vec::new(); FORMATS.len()];
            let mut diagnostics = Vec::new();
            let mut outputs = ReportOutputs {
                formats: FORMATS
                    .iter()
                    .zip(outputs_bytes.iter_mut())
                    .map(|(outfmt, bytes)| FormatOutput {
                        outfmt,
                        sink: OutputSink::Writer(bytes),
                    })
                    .collect(),
                diagnostics: &mut diagnostics,
                hits: None,
                observer: None,
            };
            render(
                &mut outputs,
                &queries,
                &subjects,
                crate::algorithm::tblastn::x_subject_store::SubjectLetters::Records(&subject_seqs),
                &ShownSubjects::Fasta { path: subject_path },
                SubjectSetNumbers::default(),
                &mut results,
                &parameters,
                &ungapped_karlin,
                &query_validity,
                &vec![false; queries.len()],
                LocalStageDScoring::default(),
                "BLOSUM62",
                Default::default(),
                1,
                Some(&seg),
                false,
                &[],
                &ReportRanges {
                    queries: &Default::default(),
                    subjects: &Default::default(),
                    epilog: true,
                },
                super::super::blast_xml::XmlValues {
                    expect: 10.0,
                    soft_masking: false,
                    gapped: true,
                    cbs: 2,
                    db_gencode: 1,
                    identity_warned: 0,
                },
                Default::default(),
                None,
                None,
                &Default::default(),
            )
            .unwrap();
            drop(outputs);
            for (fmt, actual) in FORMATS.iter().zip(&outputs_bytes) {
                let expected = std::fs::read(format!(
                    "{}/../docs/evidence/tlosan_stage_e/initial_20260925/multi_hsp.mode{mode}.outfmt{fmt}",
                    env!("CARGO_MANIFEST_DIR"),
                ))
                .unwrap();
                if *actual != expected {
                    // No NCBI counterpart: test plumbing; the platform's temporary directory
                    // replaces the Unix-only "/tmp" that made a mismatch fail on Windows with
                    // "path not found" instead of reporting the byte difference.
                    let path = std::env::temp_dir()
                        .join(format!("tblastn_stagee_mode{mode}_fmt{fmt}.actual"))
                        .display()
                        .to_string();
                    std::fs::write(&path, actual).unwrap();
                    let offset = actual
                        .iter()
                        .zip(&expected)
                        .position(|(a, b)| a != b)
                        .unwrap_or(actual.len().min(expected.len()));
                    panic!("mode={mode} outfmt={fmt} first byte difference at {offset}, actual saved to {path}");
                }
            }
        }
    }
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:1490-1589;
// c++/src/objtools/align_format/showalign.cpp:1600-1625:
// AcknowledgeBlastQuery; x_DisplayDeflines; DisplaySeqalign;
// the formatted result contains full aligned protein strings and the
// translated subject nucleotide endpoints.
/// The final HSP list with rendered alignments, in the order of the tabular rows.
#[allow(clippy::too_many_arguments)]
fn pairwise_hits(
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    hitlists: &[KappaResultHitList],
    scoring: LocalStageDScoring,
    genetic_code: u8,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    ranges: &ReportRanges<'_>,
) -> Result<Vec<PairwiseHit>> {
    let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
    let mut hits = Vec::new();
    let mut x_window = Vec::new();
    for (q_idx, hitlist) in hitlists.iter().enumerate() {
        let query = query_records[q_idx].seq();
        // NCBI c++/src/algo/blast/format/blast_format.cpp:1541-1557:
        // results.GetMaskedQueryRegions(masklocs);
        // CDisplaySeqalign(..., &masklocs, ...);
        // display.SetSeqLocChar(CDisplaySeqalign::eLowerCase);
        let report_masks = encode_tblastn_lookup_query(query, seg, mask_lowercase).seg_masks;
        for list in hitlist.lists() {
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            let subject = subject_records
                .get(oid)
                .context("TBLASTN subject OID out of range")?;
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            // The subject's letters, or (EXPERIMENT, LOSAT_X_TBNSTORE) its packed form, of which
            // the translations and the rows below read windows.
            let subject_seq = subject_letters.nt(oid);
            let mut translation = TargetTranslation::of(subject_seq, &code);
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
                let (qseq, sseq) = aligned_sequences(
                    query,
                    subject_seq,
                    &mut x_window,
                    &code,
                    &mut translation,
                    hsp.frame,
                    hsp.q_start,
                    hsp.s_start,
                    hsp.s_end,
                    &payload.edit_script,
                    &report_masks,
                )?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2120-2154
                // ```c++
                // void CDisplaySeqalign::x_FillIdentityInfo(const string& sequence_standard,
                //                                           const string& sequence ,
                //                                           int& match, int& positive,
                //                                           string& middle_line)
                // {
                //     match = 0;
                //     positive = 0;
                //     int min_length=min<int>((int)sequence_standard.size(), (int)sequence.size());
                //     if(m_AlignOption & eShowMiddleLine){
                //         middle_line = sequence;
                //     }
                //     for(int i=0; i<min_length; i++){
                //         if(sequence_standard[i]==sequence[i]){
                //             if(m_AlignOption & eShowMiddleLine){
                //                 if(m_MidLineStyle == eBar ) {
                //                     middle_line[i] = '|';
                //                 } else if (m_MidLineStyle == eChar){
                //                     middle_line[i] = sequence[i];
                //                 }
                //             }
                //             match ++;
                //         } else {
                //             if ((m_AlignType&eProt)
                //                 && m_Matrix[(int)sequence_standard[i]][(int)sequence[i]] > 0){
                //                 positive ++;
                //                 if(m_AlignOption & eShowMiddleLine){
                //                     if (m_MidLineStyle == eChar){
                //                         middle_line[i] = '+';
                //                     }
                //                 }
                //             } else {
                //                 if (m_AlignOption & eShowMiddleLine){
                //                     middle_line[i] = ' ';
                //                 }
                //             }
                // ```
                let report_positives =
                    crate::utils::matrix::protein_display_positives(&qseq, &sseq, scoring.matrix);
                // The rows are read from the searched letters; the report prints record
                // coordinates (`record_coordinates`).
                let (q_start, q_end, s_start, s_end) = record_coordinates(
                    ranges,
                    q_idx,
                    oid,
                    (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                )?;
                // NCBI reference: c++/src/objtools/align_format/showalign.cpp:411-422
                // ```c++
                // static int s_GetFrame (int start, ENa_strand strand, const CSeq_id& id,
                //                        CScope& sp)
                // {
                //     int frame = 0;
                //     if (strand == eNa_strand_plus) {
                //         frame = (start % 3) + 1;
                //     } else if (strand == eNa_strand_minus) {
                //         frame = -(((int)sp.GetBioseqHandle(id).GetBioseqLength() - start - 1)
                //                   % 3 + 1);
                //
                //     }
                //     return frame;
                // }
                // ```
                // The report's frame comes from the record coordinate of the subject row's
                // start and the record's length (the search's frame without
                // `-subject_loc`).
                let subject_length = ranges.subjects.length(oid, subject_seq.len());
                let subject_frame = crate::blastinput::seq_range::record_frame(
                    frame > 0,
                    usize::try_from(s_start - 1)?,
                    subject_length,
                );
                let hit = Hit {
                    identity: if payload.align_length > 0 {
                        100.0 * payload.report_num_ident as f64 / payload.align_length as f64
                    } else {
                        0.0
                    },
                    length: payload.align_length,
                    mismatch: payload.report_mismatches as i64,
                    gapopen: payload.gap_opens,
                    q_start: usize::try_from(q_start)?,
                    q_end: usize::try_from(q_end)?,
                    s_start: usize::try_from(s_start)?,
                    s_end: usize::try_from(s_end)?,
                    e_value: linked.evalue,
                    bit_score: payload.bit_score,
                    num_ident: payload.report_num_ident,
                    query_frame: 0,
                    query_length: query.len(),
                    q_idx: u32::try_from(q_idx)?,
                    s_idx: u32::try_from(oid)?,
                    raw_score: hsp.score,
                    sort_query_offset: usize::try_from(hsp.q_start)?,
                    sort_query_end: usize::try_from(hsp.q_end)?,
                    sort_subject_offset: usize::try_from(hsp.s_start)?,
                    sort_subject_end: usize::try_from(hsp.s_end)?,
                    has_sort_offsets: true,
                    gap_info: Some(payload.edit_script.clone()),
                    num_positives: report_positives,
                };
                hits.push(PairwiseHit {
                    hit,
                    query_seq: Some(qseq),
                    subject_seq: Some(sseq),
                    query_frame: None,
                    subject_frame: Some(i8::try_from(subject_frame)?),
                    positives: Some(report_positives),
                    gaps: Some(payload.gap_letters),
                    subject_length: Some(subject_length),
                    // The title after its first word (the `bio` description of the inputs
                    // that `bio` read alike), as BLASTN's and TBLASTX's; the reports take the
                    // subject's title bytes by subject index (`write_tblastn_pairwise_report`).
                    subject_title: subject.title.iter().position(|&byte| byte == b' ').map(
                        |space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned(),
                    ),
                    // NCBI core/blast_kappa.c:331-342:
                    // eDontAdjustMatrix -> 0; eCompoScaleOldMatrix -> 1;
                    // every other adjusted matrix -> 2.
                    comp_adjust_method: Some(match payload.matrix_adjust_rule {
                        EMatrixAdjustRule::DontAdjustMatrix => 0,
                        EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                        _ => 2,
                    }),
                    // NCBI core/link_hsps.c:1765-1810 sets the HSP
                    // linked count; showalign.cpp:3595-3598 prints sum_n.
                    sum_n: Some(linked.num),
                    subject_stored_gaps: Vec::new(),
                });
            }
        }
    }
    Ok(hits)
}

// NCBI c++/src/algo/blast/format/blast_format.cpp:1490-1589;
// c++/src/objtools/align_format/showalign.cpp:1600-1625:
// AcknowledgeBlastQuery; x_DisplayDeflines; DisplaySeqalign;
#[allow(clippy::too_many_arguments)]
fn write_pairwise(
    writer: &mut impl Write,
    hits: &[PairwiseHit],
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    subject_set: SubjectSetNumbers,
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    report_counts: crate::blastinput::report_counts::ReportCounts,
    ranges: &ReportRanges<'_>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> Result<()> {
    ensure!(
        ungapped_karlin.len() == query_records.len() && query_validity.len() == query_records.len(),
        "TBLASTN query statistics count mismatch"
    );
    let spec = ProteinScoringSpec {
        matrix: scoring.matrix,
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
    };
    // The set's numbers: a database's, or those of the subject records (the whole set).
    let (num_sequences, total_nt) = match shown {
        ShownSubjects::Fasta { .. } => (
            subject_records.len(),
            subject_letters.lengths().sum::<usize>(),
        ),
        ShownSubjects::Database { pairwise, .. } => (
            usize::try_from(pairwise.num_sequences)?,
            usize::try_from(pairwise.total_letters)?,
        ),
    };
    // The Gumbel block's database length takes `-dbsize` (`SubjectSetNumbers`).
    let gumbel = subject_set.gumbel_block(&spec, total_nt as u64);
    if gumbel.is_none() && subject_set.searchsp.is_none() {
        anyhow::bail!("TBLASTN pairwise Spouge parameters are unavailable");
    }
    let report = BlastpPairwiseReport {
        version: "2.17.0+".into(),
        database_name: shown.name(),
        database_num_sequences: num_sequences,
        database_total_letters: total_nt,
        // NCBI blast_format.cpp:2267: options.GetMatrixName(), as typed.
        matrix_name: matrix_name.to_string(),
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        word_threshold: scoring.threshold,
        window_size: scoring.window,
        gapped_karlin: lookup_protein_params(&spec),
        gumbel,
        // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2913-2927
        // ```c
        //     m_NumDescriptions = m_DfltNumDescriptions;
        //     m_NumAlignments = m_DfltNumAlignments;
        // ...
        //     if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
        //         m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
        //         m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
        // ```
        // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221);
        // `blastinput::report_counts` applies the -num_descriptions, -num_alignments and
        // -max_target_seqs values.
        num_descriptions: report_counts.num_descriptions,
        num_alignments: report_counts.num_alignments,
    };
    let queries: Vec<_> = query_records
        .iter()
        .zip(ungapped_karlin)
        .zip(&parameters.lengths)
        .enumerate()
        .map(|(q_idx, ((query, &karlin), length))| BlastpPairwiseQuery {
            // TBLASTN's report reads the validity from its own arrays.
            valid: true,
            batch_skipped: false,
            // TBLASTN's report reads the query's title bytes by query index (`query_titles`).
            query_name: String::new(),
            // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:742-744
            // ```c++
            //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
            //             out << "\nLength=";
            //             out << cbs.GetInst().GetLength() <<"\n";
            // ```
            // The query record's length, also with `-query_loc`.
            query_length: ranges.queries.length(q_idx, query.seq().len()),
            ungapped_karlin: karlin,
            effective_search_space: length.eff_searchsp,
        })
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
    let subject_titles: Vec<Arc<[u8]>> = subject_records
        .iter()
        .map(|record| Arc::from(record.title()))
        .collect();
    let config = PairwiseConfig {
        program: "tblastn".into(),
        protein_matrix: scoring.matrix,
        ..PairwiseConfig::default()
    };
    write_tblastn_pairwise_report(
        hits,
        writer,
        &config,
        &queries,
        &query_titles,
        query_validity,
        query_batch_skipped,
        &subject_titles,
        shown.parsed_subjects(),
        &report,
        shown.pairwise_database(),
        probe,
        warnings,
        ranges.epilog,
    )?;
    Ok(())
}

// NCBI c++/src/objtools/align_format/tabular.cpp:971-1021;
// c++/src/objtools/align_format/showalign.cpp:2122-2149:
// alnVec->GetWholeAlnSeqString(0,m_QuerySeq);
// alnVec->GetWholeAlnSeqString(1,m_SubjectSeq);
// the edit script consumes query and translated subject residues in order.
// NCBI reference: c++/src/objtools/alnmgr/alnvec.cpp:911-918
// ```c++
//     int state = 0;
//     size_t aa_i = 0;
//     for (size_t na_i = 0; na_i < na_size; ) {
//         for (size_t i = 0; i < 3; i++) {
//             state = tbl.NextCodonState(state, na[na_i++]);
//         }
//         aa[aa_i++] = tbl.GetCodonResidue(state);
//     }
// ```
// The subject row is translated again from the nucleotides for the display, with the
// table of `CTrans_table` (an ambiguous codon is B, Z or J where its residues are D/N, E/Q
// or I/L; TBLASTX's `display_residue`), not taken from the search's translation.
#[allow(clippy::too_many_arguments)]
fn aligned_sequences(
    query: &[u8],
    // EXPERIMENT (LOSAT_X_TBNSTORE): the letters, or a packed subject.
    subject: SubjectNt<'_>,
    window_buffer: &mut Vec<u8>,
    code: &GeneticCode,
    target: &mut TargetTranslation<'_>,
    frame: i8,
    q_start: i32,
    s_start: i32,
    s_end: i32,
    script: &[GapEditOp],
    report_masks: &[(usize, usize)],
) -> Result<(String, String)> {
    ensure!(!script.is_empty(), "TBLASTN report HSP has no edit script");
    let (translated, _, base) = target.get(frame, s_start, s_end)?;
    let mut q = usize::try_from(q_start)?;
    let mut s = super::search_gapped::window_index(s_start, base)
        .context("TBLASTN report translation window starts after HSP")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1508-1516
    // ```c
    //     translation_frame = range->context;
    //     if (translation_frame > 0) {
    //         translation_start = 3 * range->begin;
    //     } else {
    //         translation_start =
    //             self->length - 3 * range->end + translation_frame + 1;
    //     }
    //     num_nucleotides =
    //         3 * (range->end - range->begin) + ABS(translation_frame) - 1;
    // ```
    // The nucleotides the row is translated from: the record's letters, or (EXPERIMENT,
    // LOSAT_X_TBNSTORE) the letters of the HSP's codons made from the store, the window NCBI
    // translates for the residues the script reads (from `s_start` on). Its codon of residue `o`
    // is the one `display_residue` reads in the whole subject: from `3 * (o - s_start) +
    // |frame| - 1` on the plus strand, and as far from the window's end as from the subject's on
    // the minus strand.
    let (nucleotides, first) = match subject {
        SubjectNt::Letters(letters) => (letters, 0),
        SubjectNt::Packed(_) => {
            let first = usize::try_from(s_start)?;
            let residues: usize = script
                .iter()
                .filter(|op| !matches!(op, GapEditOp::Ins(_)))
                .map(|op| op.num() as usize)
                .sum();
            let frame_shift = usize::from(frame.unsigned_abs());
            let length = 3 * residues + frame_shift - 1;
            let start = if frame > 0 {
                Some(3 * first)
            } else {
                (subject.len() + 1).checked_sub(3 * (first + residues) + frame_shift)
            }
            .filter(|start| start + length <= subject.len())
            .context("TBLASTN report HSP outside the subject")?;
            (subject.fetch(start, start + length, window_buffer), first)
        }
    };
    let display = crate::algorithm::tblastx::report::DisplayCode::new(code);
    let displayed_subject = |s: usize| -> u8 {
        display.residue(nucleotides, frame, (base + s as isize) as usize - first)
    };
    const NCBISTDAA_TO_AA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    let mut query_string = String::new();
    let mut subject_string = String::new();
    for &op in script {
        for _ in 0..op.num() {
            match op {
                GapEditOp::Sub(_) => {
                    let &qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    let &sa = translated
                        .get(s)
                        .context("TBLASTN report subject offset overflow")?;
                    // NCBI blast_format.cpp:1541-1557 renders filter locations
                    // as lowercase in the pairwise query sequence.
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:184-187
                    // ```c++
                    // inline unsigned char s_ASCII_MustBeLowerToUpper(unsigned char c)
                    // {
                    //     return c + ('A' - 'a');
                    // }
                    // ```
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:945-949
                    // ```c++
                    //         case eCharType_MaskedNonGap:
                    //             CloseGap(pos == 0);
                    //             m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
                    //             OpenMask();
                    //             ++m_CurrentPos;
                    // ```
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
                    // ```c++
                    //                     } else if (m_SeqLocChar==eLowerCase){
                    //                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
                    // ```
                    let qa = qa.to_ascii_uppercase();
                    let displayed = if report_masks
                        .iter()
                        .any(|&(left, right)| left <= q && q < right)
                    {
                        qa.to_ascii_lowercase()
                    } else {
                        qa
                    };
                    query_string.push(char::from(displayed));
                    ensure!(
                        NCBISTDAA_TO_AA.get(sa as usize).is_some(),
                        "invalid TBLASTN subject amino acid"
                    );
                    subject_string.push(char::from(displayed_subject(s)));
                    q += 1;
                    s += 1;
                }
                GapEditOp::Del(_) => {
                    let &sa = translated
                        .get(s)
                        .context("TBLASTN report subject offset overflow")?;
                    query_string.push('-');
                    ensure!(
                        NCBISTDAA_TO_AA.get(sa as usize).is_some(),
                        "invalid TBLASTN subject amino acid"
                    );
                    subject_string.push(char::from(displayed_subject(s)));
                    s += 1;
                }
                GapEditOp::Ins(_) => {
                    let &qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:184-187
                    // ```c++
                    // inline unsigned char s_ASCII_MustBeLowerToUpper(unsigned char c)
                    // {
                    //     return c + ('A' - 'a');
                    // }
                    // ```
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:945-949
                    // ```c++
                    //         case eCharType_MaskedNonGap:
                    //             CloseGap(pos == 0);
                    //             m_SeqData[m_CurrentPos] = s_ASCII_MustBeLowerToUpper(c);
                    //             OpenMask();
                    //             ++m_CurrentPos;
                    // ```
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
                    // ```c++
                    //                     } else if (m_SeqLocChar==eLowerCase){
                    //                         actualSeq[i-start]=tolower((unsigned char) actualSeq[i-start]);
                    // ```
                    let qa = qa.to_ascii_uppercase();
                    let displayed = if report_masks
                        .iter()
                        .any(|&(left, right)| left <= q && q < right)
                    {
                        qa.to_ascii_lowercase()
                    } else {
                        qa
                    };
                    query_string.push(char::from(displayed));
                    subject_string.push('-');
                    q += 1;
                }
            }
        }
    }
    Ok((query_string, subject_string))
}

/// `aligned_sequences` into two strings the caller keeps, for the outfmt 5 rows, which are
/// made for one HSP at a time. It calls `aligned_sequences` as it is (the rows of 0/6/7 are
/// made by that function alone, unchanged since the base of OF4b, DECISIONS 23) and moves its
/// strings into the caller's.
#[allow(clippy::too_many_arguments)]
pub(super) fn aligned_sequences_into(
    query: &[u8],
    subject: SubjectNt<'_>,
    window_buffer: &mut Vec<u8>,
    code: &GeneticCode,
    target: &mut TargetTranslation<'_>,
    frame: i8,
    q_start: i32,
    s_start: i32,
    s_end: i32,
    script: &[GapEditOp],
    report_masks: &[(usize, usize)],
    query_string: &mut String,
    subject_string: &mut String,
) -> Result<()> {
    (*query_string, *subject_string) = aligned_sequences(
        query,
        subject,
        window_buffer,
        code,
        target,
        frame,
        q_start,
        s_start,
        s_end,
        script,
        report_masks,
    )?;
    Ok(())
}
