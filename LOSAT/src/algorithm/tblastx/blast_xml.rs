//! TBLASTX's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit list
//! (`report::blast_xml`), one query at a time, once per run.
//!
//! TBLASTX is an ungapped search, so NCBI makes one `diags` Seq-align per subject of a query's
//! final hit list, one Std-seg per HSP (`BLASTUngappedHspListToSeqAlign`, see `seqalign`), and
//! the XML writer takes each Std-seg as an `Hsp` of the subject's `Hit`
//! (`PrepareBlastUngappedSeqalign`): the final list (`report_order`) holds the HSPs in that
//! order, query by query, one run of HSPs per subject. What `s_SeqAlignSetToXMLHsps` reads of an
//! HSP is its Seq-align (the coordinates in the records, lower first, the score, the bit score,
//! the e-value or the sum e-value of a linked set) and the two rows of the alignment, from which
//! it counts identities and positives again: nothing of the engine's counts and nothing of
//! outfmt 0's midline. The rows are the ones outfmt 0 shows (`SubjectLetters::rows_into`, both
//! translated again from the nucleotides with B, Z and J for the ambiguous codons) without the
//! lower case of the masks; the masks are applied to the query row as NCBI's `s_MaskQuerySeq`
//! does (`report::blast_xml::hsp`), over the SEG intervals of each frame joined as NCBI joins
//! them (`report::xml_query_masks`).
//!
//! The XML2 and JSON documents of outfmt 15 and 16 (`report::blast_xml2`) are written from the
//! same per-query records by `Xml2Writer`: `write_tblastx_xml` is generic over the writer
//! (`StructuredXml`). They add the masks of every searched query (also one without hits; a
//! query of a batch NCBI did not search has none), the `searched` flag of the query's batch,
//! the taxid of a database sequence, and the genetic codes (TBLASTX is ungapped and has no
//! composition-based statistics: no `gap-open`, `gap-extend` or `cbs`).
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

use anyhow::Result;

use super::blast_engine::{TblastxHsp, TblastxQueryStats, TblastxReportRun};
use super::ncbi_cutoffs::try_compute_eff_lengths_tblastx;
use super::subjects::SubjectSetNumbers;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    filter, QueryMask, StructuredXml, SubjectNameCache, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
    XmlSubjects, REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::{Xml2Run, Xml2Target};
use crate::report::query_warnings::{QueryWarnings, INVALID_QUERY_MESSAGE};
use crate::report::serial::seqalign::seqalign_evalue;
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTX";

/// `-gapopen 11 -gapextend 1`: TBLASTX has no options for them, and NCBI's XML prints its option
/// defaults (blastxml_format.cpp:602-607; DECISIONS 19, Q7), or the best costs of `-matrix`
/// (`gap_costs`).
const GAP_OPEN: i32 = 11;
const GAP_EXTEND: i32 = 1;

/// The gap costs of the XML's parameters: with `-matrix`, the matrix's best costs.
///
/// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:257-273
/// ```c
///     int gap_open=0, gap_extend=0;
///     if (args.Exist(kArgMatrixName) && args[kArgMatrixName])
///          BLAST_GetProteinGapExistenceExtendParams
///              (args[kArgMatrixName].AsString().c_str(), &gap_open, &gap_extend);
///
///     if (args.Exist(kArgGapOpen) && args[kArgGapOpen]) {
///         opt.SetGapOpeningCost(args[kArgGapOpen].AsInteger());
///     }
///     else if (args.Exist(kArgMatrixName) && args[kArgMatrixName]) {
///         opt.SetGapOpeningCost(gap_open);
///     }
/// ```
fn gap_costs(matrix: Option<&str>) -> (i32, i32) {
    match matrix {
        Some(name) => crate::stats::protein_options::protein_gap_existence_extend_params(name)
            .unwrap_or((0, 0)),
        None => (GAP_OPEN, GAP_EXTEND),
    }
}

/// The messages of a query's result (the warnings the other formats write to standard error,
/// `report::query_warning`): the Karlin-Altschul message of a batch NCBI did not search (a
/// translated query posts it only when no context of the batch is valid) and the message of a
/// record without letters. TBLASTX posts nothing for an invalid query of a searched batch
/// (blast_stat.c:2783-2790).
fn query_messages(query: &FastaRecord, searched: bool) -> Vec<String> {
    let mut messages = Vec::new();
    if !searched {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    if query.seq().is_empty() {
        messages.push(crate::algorithm::blastn::blast_engine::NO_DATA_MESSAGE.to_string());
    }
    messages
}

/// The `Statistics` of a query: the first valid context's length adjustment, effective search
/// space and ideal ungapped Karlin block (TBLASTX is ungapped); the unsearched values for a query
/// without a valid context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-100
/// ```c++
///     // find the first valid context corresponding to this query
///     for (i = 0; i < context_per_query; i++) {
///         ...
///         if (ctx->is_valid) {
///             m_SearchSpace = ctx->eff_searchsp;
/// 	    m_LengthAdjustment = ctx->length_adjustment;
///             break;
/// ```
/// The search keeps that context's block, search space and length, not its length adjustment: it
/// is computed again here with the search's own function (`try_compute_eff_lengths_tblastx`,
/// the call of `search_query_batch` with the same subject set), and the search space it gives
/// must be the kept one.
fn query_stats(stats: &TblastxQueryStats, numbers: &SubjectSetNumbers) -> XmlStats {
    let Some(karlin) = &stats.karlin else {
        return XmlStats::NONE;
    };
    let lengths = try_compute_eff_lengths_tblastx(
        stats.aa_len as i64,
        numbers.effective_length(),
        numbers.count,
        karlin,
    );
    debug_assert_eq!(
        lengths.map(|lengths| lengths.eff_searchsp),
        Some(stats.eff_searchsp),
        "the effective search space of the first valid context"
    );
    XmlStats {
        hsp_len: lengths.map_or(0, |lengths| lengths.length_adjustment as i32),
        eff_space: stats.eff_searchsp,
        kappa: karlin.k,
        lambda: karlin.lambda,
        entropy: karlin.h,
    }
}

/// What the XML2 and JSON documents need of the run beyond outfmt 5's (`Xml2Writer::new`):
/// the target, and the genetic codes (`Parameters`; TBLASTX is ungapped and has no
/// composition-based statistics).
pub(crate) fn xml2_run(run: &TblastxReportRun<'_>) -> Xml2Run {
    Xml2Run {
        target: match &run.database {
            Some(_) => Xml2Target::Database,
            None => Xml2Target::Subjects(run.subjects.len()),
        },
        gapped: false,
        cbs: 0,
        query_gencode: Some(i32::from(run.args.query_gencode)),
        db_gencode: Some(i32::from(run.args.db_gencode)),
    }
}

/// Writes the document of every query of `run` (`hits` is the final list after
/// `shift_to_records`: record coordinates, in query order) through `document` (BLAST XML, or
/// the XML2 / JSON documents) to `writer`, each query's warnings before it, with the document's
/// header before the first query and its end after the last when `run.epilog`. Everything NCBI
/// cannot write is refused before the first byte: the `-db` name, the query titles and
/// messages, the titles of the subjects with hits. `all_masks`: the masks of every searched
/// query are listed (XML2 and JSON), not those of the queries with hits alone (outfmt 5).
/// `run.xml_taxa` holds the taxids of the database sequences that XML2 and JSON show.
pub(crate) fn write_tblastx_xml<W: StructuredXml>(
    document: &mut W,
    all_masks: bool,
    writer: &mut dyn Write,
    hits: &[TblastxHsp],
    run: &TblastxReportRun<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let args = run.args;
    let taxa = run.xml_taxa;
    let (db_name, db_num, db_len) = match &run.database {
        Some(database) => (
            database.name,
            i32::try_from(database.pairwise.num_sequences)?,
            i64::try_from(database.pairwise.total_letters)?,
        ),
        None => ("", 0, 0),
    };
    names::check_database_name(PROGRAM, db_name)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2206-2208 (XML2
    // and JSON) and :2242-2243 (BLAST XML)
    // ```c++
    //     if ((m_FormatType == CFormattingArgs::eXml2) || (m_FormatType == CFormattingArgs::eJson) ||
    //         (m_FormatType == CFormattingArgs::eXml2_S) || (m_FormatType == CFormattingArgs::eJson_S)) {
    //         if(!m_AccumulatedResults.empty()) {
    // ...
    //     if (m_FormatType == CFormattingArgs::eXml) {
    //         m_Outfile << m_BlastXMLIncremental->m_SerialXmlEnd << endl;
    // ```
    // The batch loop refused a name that the formatting table does not know when a batch has
    // letters (`run_batch_loop`); a run without queries makes no report data, and the document
    // is its header and end.
    let matrix = args.matrix_name();
    let positives = FormattingMatrix::new(matrix);
    anyhow::ensure!(
        positives.is_some() || hits.is_empty(),
        "no formatting matrix for {matrix}"
    );
    let (gap_open, gap_extend) = gap_costs(args.matrix.as_deref());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:414-417
    // ```c
    //     if (SBlastFilterOptionsMaskAtHash(filtering_options)) {
    //         if (strlen(retval) != 0) {
    //             /* Add mask at hash as a modifier for other filtering options */
    //             if (!s_SafeStrCat(&retval, &retval_size, "m;")) {
    // ```
    // `-soft_masking true` (`mask_at_hash`) ends the filter string with `m;`.
    let filter = filter::filter_string(
        None,
        args.seg_spec()?.filter_options(),
        None,
        args.soft_masking(),
    );
    let xml_run = XmlRun {
        program: "tblastx",
        version: "TBLASTX 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name: db_name.as_bytes(),
        matrix: Some(matrix),
        expect: args.evalue,
        gap_open,
        gap_extend,
        reward_penalty: None,
        filter: &filter,
        db_num,
        db_len,
        positives: positives.as_ref(),
        mask_char: b'X',
    };
    let queries = run.queries;
    let query_placements = &run.ranges.input.placements;
    let subject_placements = run.ranges.subject_placements;
    // The hits of each query: the list is in query order.
    let mut query_hits: Vec<std::ops::Range<usize>> = Vec::with_capacity(queries.len());
    let mut next = 0usize;
    for q_idx in 0..queries.len() {
        let start = next;
        while next < hits.len() && hits[next].hit.q_idx as usize == q_idx {
            next += 1;
        }
        query_hits.push(start..next);
    }
    debug_assert_eq!(next, hits.len(), "the final list is in query order");
    // The definition line and message of every query, and the names of the subjects with
    // hits, made (and refused) before the first byte.
    let texts = queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| {
            names::query_text(
                PROGRAM,
                &query.local_id,
                query.title(),
                false,
                &query_messages(query, run.searched[q_idx]),
                !query_hits[q_idx].is_empty(),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let subject_length =
        |s_idx: usize| subject_placements.length(s_idx, run.subject_letters.length(s_idx));
    let subjects = match &run.database {
        Some(database) => XmlSubjects::Database {
            ordinals: &database.reported.ordinals,
            titles: crate::report::blast_xml::XmlTitles::Records(run.subjects),
            // Only XML2 shows the taxid (`taxa` is empty for outfmt 5).
            taxids: &taxa.taxids,
            taxonomy: &taxa.taxonomy,
            length: &subject_length,
            parsed: &[],
            parsed_taxonomy: &[],
        },
        None => XmlSubjects::Records {
            records: run.subjects,
            length: &subject_length,
            protein: false,
        },
    };
    let mut cache = SubjectNameCache::new(&subjects);
    let mut last_subject = None;
    for hsp in hits {
        if last_subject != Some(hsp.hit.s_idx) {
            last_subject = Some(hsp.hit.s_idx);
            cache.get(hsp.hit.s_idx as usize, PROGRAM)?;
        }
    }
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let mut taxdb_printer = taxa.notes.printer();
    let mut buf: Vec<u8> = Vec::new();
    let mut masks: Vec<QueryMask> = Vec::new();
    let (mut window, mut query_row, mut subject_row): (Vec<u8>, Vec<u8>, Vec<u8>) =
        (Vec::new(), Vec::new(), Vec::new());
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let text = &texts[q_idx];
        let searched_length = query.seq().len();
        let record_length = query_placements.length(q_idx, searched_length) as i32;
        let query_hits = &hits[query_hits[q_idx].clone()];
        masks.clear();
        // A query of a batch NCBI did not search has no masks (its result is made without the
        // search); XML2 and JSON list the masks of every other query, also one without hits.
        if run.searched[q_idx] && (all_masks || !query_hits.is_empty()) {
            crate::algorithm::tblastx::report::xml_query_masks(
                &run.query_stats[q_idx],
                searched_length,
                query_placements.offset(q_idx),
                &mut masks,
            );
        }
        let xml_query = XmlQuery {
            id: query.local_id.as_bytes(),
            title: &text.title,
            len: i32::try_from(searched_length)?,
            stats: query_stats(&run.query_stats[q_idx], run.numbers),
            message: text.message.as_deref(),
            // Read by XML2 alone.
            searched: run.searched[q_idx],
            masks: &masks,
        };
        let written = document.query(&xml_run, &xml_query, &mut buf, |sink| {
            let mut first = 0usize;
            while first < query_hits.len() {
                let s_idx = query_hits[first].hit.s_idx;
                let mut end = first;
                while end < query_hits.len() && query_hits[end].hit.s_idx == s_idx {
                    end += 1;
                }
                let subject = cache.get(s_idx as usize, PROGRAM)?;
                let subject_record_length = subjects.length(s_idx as usize);
                let hit = XmlHit::of(subject, subject_record_length);
                // NCBI's messages about the taxonomy database: opening and the volume's lookup
                // at the first `Hit` that needs them, the lookup of the `Hit` every time
                // (nothing for the other formats).
                if let Some(warnings) = warnings.as_deref_mut() {
                    taxdb_printer.at_hit(
                        s_idx as usize,
                        subject.taxid,
                        &subject.taxonomy,
                        &mut *warnings.sink,
                    )?;
                }
                sink.hit(&hit, |hsps| {
                    for hsp in &query_hits[first..end] {
                        run.subject_letters.rows_into(
                            hsp,
                            query.seq(),
                            &query_code,
                            &db_code,
                            &mut window,
                            &mut query_row,
                            &mut subject_row,
                        );
                        let hit = &hsp.hit;
                        // The Seq-align's endpoints, lower first: the tabular order of a
                        // minus frame is not the XML's.
                        let (query_from, query_to) = (
                            hit.q_start.min(hit.q_end) as i32,
                            hit.q_start.max(hit.q_end) as i32,
                        );
                        let (hit_from, hit_to) = (
                            hit.s_start.min(hit.s_end) as i32,
                            hit.s_start.max(hit.s_end) as i32,
                        );
                        let plus = hit.query_frame > 0;
                        hsps.hsp(&XmlHsp {
                            score: hit.raw_score,
                            bit_score: hit.bit_score,
                            evalue: seqalign_evalue(hit.e_value),
                            query_from,
                            query_to,
                            hit_from,
                            hit_to,
                            query_frame: translation_frame(
                                plus,
                                query_from,
                                query_to,
                                record_length,
                            ),
                            hit_frame: translation_frame(
                                hsp.subject_frame > 0,
                                hit_from,
                                hit_to,
                                subject_record_length,
                            ),
                            qrow: &query_row,
                            srow: &subject_row,
                            ops: None,
                            // The first residue of the row: the lowest nucleotide of the
                            // interval on the plus strand, the highest on the minus strand,
                            // where the row runs backwards.
                            query_first: if plus { query_from - 1 } else { query_to - 1 },
                            query_step: if plus { 3 } else { -3 },
                        })?;
                    }
                    Ok(())
                })?;
                first = end;
            }
            Ok(())
        });
        // A failed query leaves its partial bytes in `buf`, which the error return drops.
        written?;
        writer.write_all(&buf)?;
        buf.clear();
    }
    if run.epilog {
        document.finish(&mut buf)?;
        writer.write_all(&buf)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stats::lookup_protein_params_ungapped;

    // NCBI c++/src/algo/blast/api/blast_results.cpp:72-115: a query without a valid context has
    // the -1 blocks and zero lengths.
    #[test]
    fn a_query_without_a_valid_context_has_the_unsearched_statistics() {
        let numbers = SubjectSetNumbers {
            total: 3_000,
            count: 3,
            dbsize: None,
        };
        assert_eq!(
            query_stats(&TblastxQueryStats::default(), &numbers),
            XmlStats::NONE
        );
    }

    // The length adjustment is the search's own function of the first valid context's length and
    // the subject set (`-dbsize` replaces the length, not the number of sequences); the search
    // space is the one the search kept.
    #[test]
    fn the_length_adjustment_is_computed_again_from_the_context_and_the_subject_set() {
        let karlin = lookup_protein_params_ungapped(crate::config::ScoringMatrix::Blosum62);
        for (dbsize, count) in [(None, 3), (Some(5_000_000), 3), (None, 1)] {
            let numbers = SubjectSetNumbers {
                total: 3_000,
                count,
                dbsize,
            };
            let lengths = try_compute_eff_lengths_tblastx(
                130,
                numbers.effective_length(),
                numbers.count,
                &karlin,
            )
            .expect("fits");
            let stats = TblastxQueryStats {
                karlin: Some(karlin.clone()),
                eff_searchsp: lengths.eff_searchsp,
                aa_len: 130,
                context0_search_space: 0,
                seg_masks: Vec::new(),
            };
            let xml = query_stats(&stats, &numbers);
            assert_eq!(i64::from(xml.hsp_len), lengths.length_adjustment);
            assert!(xml.hsp_len > 0);
            assert_eq!(xml.eff_space, lengths.eff_searchsp);
            assert_eq!(
                (xml.kappa, xml.lambda, xml.entropy),
                (karlin.k, karlin.lambda, karlin.h)
            );
        }
    }

    // TBLASTX posts the Karlin-Altschul message only for a batch that was not searched; an empty
    // record has its own message (blast_setup.hpp:190-198).
    #[test]
    fn the_messages_of_a_query_follow_its_batch_and_its_letters() {
        let query = FastaRecord::new("Query_1", b"q", b"ACGTACGTA");
        assert!(query_messages(&query, true).is_empty());
        assert_eq!(
            query_messages(&query, false),
            vec![INVALID_QUERY_MESSAGE.to_string()]
        );
        let empty = FastaRecord::new("Query_2", b"e", b"");
        assert_eq!(query_messages(&empty, true).len(), 1);
        assert_eq!(query_messages(&empty, false).len(), 2);
    }
}
