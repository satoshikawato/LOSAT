//! NCBI TBLASTN local-subject result formatting.

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, ensure, Context, Result};
use bio::io::fasta;

use super::search_gapped::TargetTranslation;
use super::search_seed::encode_tblastn_lookup_query;
use super::stage_d_pipeline::LocalStageDScoring;
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use crate::api::local_blast::{FormatProbe, HspIndex, ReportOutputs};
use crate::common::{GapEditOp, Hit};
use crate::config::ProteinScoringSpec;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
};
use crate::report::pairwise::{
    write_tblastn_pairwise_report, BlastpPairwiseQuery, BlastpPairwiseReport, PairwiseConfig,
    PairwiseHit,
};
use crate::stats::spouge::lookup_protein_gumbel_params;
use crate::stats::tables::lookup_protein_params;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;

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
#[allow(clippy::too_many_arguments)]
/// The format number of a checked `-outfmt` value, as NCBI's `ParseFormattingString` reads
/// it (a custom specification of outfmt 0 is ignored); -1 for a value it does not read.
pub(super) fn format_number(outfmt: &str) -> i32 {
    crate::blastinput::app::parse_formatting_string(outfmt).map_or(-1, |choice| choice.number)
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
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    subject_path: &Path,
    hitlists: &mut [KappaResultHitList],
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    max_target_seqs_given: Option<usize>,
    genetic_code: u8,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    query_warning_lines: &[Vec<u8>],
    ranges: &ReportRanges<'_>,
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
    // format (`QueryWarnings`).
    let mut format_warnings = Some(crate::report::query_warnings::QueryWarnings {
        before: query_warning_lines,
        sink: &mut *outputs.diagnostics,
    });
    let observer = &mut outputs.observer;
    for (format_index, format) in outputs.formats.iter_mut().enumerate() {
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        let mut warnings = format_warnings.take();
        let mut writer = format.sink.open()?;
        match format_number(format.outfmt) {
            6 | 7 => write_tabular(
                &mut writer,
                format_number(format.outfmt) == 7,
                query_records,
                subject_records,
                subject_path,
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
                subject_path,
                parameters,
                ungapped_karlin,
                query_validity,
                query_batch_skipped,
                scoring,
                matrix_name,
                max_target_seqs_given,
                ranges,
                probe.as_mut(),
                warnings.as_mut(),
            )?,
            _ => bail!("unsupported TBLASTN outfmt {}", format.outfmt),
        }
        writer.flush()?;
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
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    subject_path: &Path,
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
            // The title's bytes (`write_outfmt7_query_line`): `bio`'s ID, a space and the
            // rest of the defline.
            let mut title = query.id().as_bytes().to_vec();
            if let Some(desc) = query.desc() {
                title.push(b' ');
                title.extend_from_slice(desc.as_bytes());
            }
            crate::report::outfmt6::write_outfmt7_query_line(writer, &title)?;
            writeln!(
                writer,
                "# Database: User specified sequence set (Input: {})",
                subject_path.display()
            )?;
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
                let subject_len = i32::try_from(subject.seq().len())?;
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
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    query.id(),
                    subject.id(),
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
fn record_coordinates(
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::algorithm::tblastn::stage_d_pipeline::{
        run_local_for_report, LocalStageDProfile, LocalStageDScoring,
    };
    use crate::api::local_blast::{FormatOutput, OutputSink};
    use crate::utils::seg::SegParams;

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
        let read = |path: &Path| {
            fasta::Reader::from_file(path)
                .unwrap()
                .records()
                .collect::<std::result::Result<Vec<_>, _>>()
                .unwrap()
        };
        let queries = read(query_path);
        let subjects = read(full_subject_path);
        let query_seqs: Vec<_> = queries.iter().map(|record| record.seq().to_vec()).collect();
        let subject_seqs: Vec<_> = subjects
            .iter()
            .map(|record| record.seq().to_vec())
            .collect();
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
                subject_path,
                &mut results,
                &parameters,
                &ungapped_karlin,
                &query_validity,
                &vec![false; queries.len()],
                LocalStageDScoring::default(),
                "BLOSUM62",
                None,
                1,
                Some(&seg),
                false,
                &[],
                &ReportRanges {
                    queries: &Default::default(),
                    subjects: &Default::default(),
                    epilog: true,
                },
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
                    let path = format!("/tmp/tblastn_stagee_mode{mode}_fmt{fmt}.actual");
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
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    hitlists: &[KappaResultHitList],
    scoring: LocalStageDScoring,
    genetic_code: u8,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    ranges: &ReportRanges<'_>,
) -> Result<Vec<PairwiseHit>> {
    let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
    let mut hits = Vec::new();
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
            let mut translation = TargetTranslation::new(subject.seq(), &code);
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                let hsp = &linked.hsp;
                let frame = i32::from(hsp.frame);
                // NCBI c++/src/algo/blast/core/blast_hits.c:1087-1105:
                // translated coordinates use three nucleotide letters per
                // residue and reverse strand endpoints run in descending order.
                let subject_len = i32::try_from(subject.seq().len())?;
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
                    subject.seq(),
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
                let subject_length = ranges.subjects.length(oid, subject.seq().len());
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
                    mismatch: payload.report_mismatches,
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
                    subject_title: subject.desc().map(str::to_owned),
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
    query_records: &[fasta::Record],
    subject_records: &[fasta::Record],
    subject_path: &Path,
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    max_target_seqs_given: Option<usize>,
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
    let total_nt: usize = subject_records
        .iter()
        .map(|record| record.seq().len())
        .sum();
    let gumbel = lookup_protein_gumbel_params(&spec, i64::try_from(total_nt / 3)?)
        .context("TBLASTN pairwise Spouge parameters are unavailable")?;
    let report = BlastpPairwiseReport {
        version: "2.17.0+".into(),
        database_name: format!(
            "User specified sequence set (Input: {})",
            subject_path.display()
        ),
        database_num_sequences: subject_records.len(),
        database_total_letters: total_nt,
        // NCBI blast_format.cpp:2267: options.GetMatrixName(), as typed.
        matrix_name: matrix_name.to_string(),
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        word_threshold: f64::from(scoring.threshold),
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
        // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221).
        num_descriptions: max_target_seqs_given.unwrap_or(500),
        num_alignments: max_target_seqs_given.unwrap_or(250),
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
            query_name: match query.desc() {
                Some(desc) => format!("{} {desc}", query.id()),
                None => query.id().to_string(),
            },
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
    let subject_ids: Vec<Arc<str>> = subject_records
        .iter()
        .map(|record| Arc::from(record.id()))
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
        query_validity,
        query_batch_skipped,
        &subject_ids,
        &report,
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
    subject: &[u8],
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
    let displayed_subject = |s: usize| -> u8 {
        crate::algorithm::tblastx::report::display_residue(
            subject,
            frame,
            (base + s as isize) as usize,
            code,
        )
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
