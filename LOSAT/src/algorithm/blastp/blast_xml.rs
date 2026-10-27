//! BLASTP's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit lists
//! (`report::blast_xml`), one query at a time.
//!
//! NCBI makes the Iteration of each query from its Seq-aligns (`BlastHitList2SeqAlign_OMF`,
//! see `seqalign`); the XML walks the lists as the Seq-annot formats do. What NCBI's
//! `s_SeqAlignSetToXMLHsps` reads of an HSP is its Seq-align: the coordinates in the records,
//! the score, the e-value and the bit score, and the two rows of the alignment from which it
//! counts identities and positives again (so, unlike the Seq-annot, nothing of the engine's
//! counts).
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

use anyhow::{bail, Context, Result};

use super::hsp::BlastpHitList;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    QueryMask, StructuredXml, SubjectNameCache, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
    XmlSubjects,
};
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::QueryNames;
use crate::report::serial::seqalign::seqalign_evalue;

/// What the BLASTP document needs besides the final hit lists.
pub(crate) struct BlastpXml<'a> {
    pub run: XmlRun<'a>,
    pub queries: &'a QueryNames<'a>,
    /// By query: the statistics (`XmlStats::NONE` for a query without a valid context).
    pub stats: Vec<XmlStats>,
    /// By query: the warning messages of its result (the O characters read as X, a query
    /// without letters, an invalid query).
    pub messages: Vec<Vec<String>>,
    pub subjects: XmlSubjects<'a>,
    /// NCBI's messages about the taxonomy database that come before the `Hit`s' own
    /// (XML2 and JSON with `-db`; nothing otherwise).
    pub taxdb_notes: crate::algorithm::common::database::TaxdbNotes,
    /// By query, the masked regions of the query in record coordinates (`QueryMask`).
    pub masks: Vec<Vec<QueryMask>>,
    /// By query: whether NCBI searched the query's batch (it has a valid context in it). Only
    /// XML2 reads it.
    pub searched: Vec<bool>,
    /// The program, for the refusals' messages (`BLASTP`).
    pub program: &'static str,
    /// Whether the document ends with this call.
    pub epilog: bool,
}

/// Writes the document of every query of `hit_lists` (indexed by query) through `document`
/// (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before it,
/// with the document's header before the first query and its end after the last when
/// `xml.epilog`. `rows(q)` gives the rows of query `q`'s HSPs, in the lists' order.
pub(crate) fn write_blastp_xml<W: StructuredXml>(
    xml: &BlastpXml<'_>,
    document: &mut W,
    hit_lists: &[Option<BlastpHitList>],
    rows: &mut dyn FnMut(usize) -> Result<Vec<PairwiseHit>>,
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let program = xml.program;
    // Everything NCBI cannot write is refused before the first byte of the document: the
    // database name, each query's definition line and message, the names of the subjects with
    // hits. This includes the last query of a `-subject` call, whose report NCBI holds back
    // until the next query's arrives and never writes when the next batch fails to read
    // (outfmt 15 and 16): LOSAT refuses a title it would not have written then. It cannot
    // wait for the report to be written, which would refuse after a byte (the reports before
    // it are written by now), and the next batch is not read yet.
    names::check_text(xml.run.db_name)
        .map_err(|refusal| refusal.error(program, "a database name"))?;
    let has_hits = |q_idx: usize| {
        hit_lists[q_idx].as_ref().is_some_and(|list| {
            list.hsplist_array
                .iter()
                .take(list.hsplist_count)
                .any(|hsp_list| !hsp_list.hsps.is_empty())
        })
    };
    let texts = names::query_texts(
        program,
        &xml.queries.records[..hit_lists.len()],
        true,
        &xml.messages,
        has_hits,
    )?;
    // The names of the subjects with hits, made once each; a title NCBI cannot write is refused
    // here, before the first byte.
    let mut subjects = SubjectNameCache::new(&xml.subjects);
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
            if let Some(hsp) = hsp_list.hsps.first() {
                subjects.get(hsp.s_idx as usize, program)?;
            }
        }
    }
    let mut taxdb_printer = xml.taxdb_notes.printer();
    let mut buf: Vec<u8> = Vec::new();
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let record = &xml.queries.records[q_idx];
        let text = &texts[q_idx];
        let query = XmlQuery {
            id: record.local_id.as_bytes(),
            title: &text.title,
            len: xml.queries.searched_lengths[q_idx] as i32,
            stats: xml.stats[q_idx],
            message: text.message.as_deref(),
            searched: xml.searched[q_idx],
            masks: &xml.masks[q_idx],
        };
        let mut hit_rows: Vec<PairwiseHit> = Vec::new();
        if hit_list.is_some() {
            hit_rows = rows(q_idx)?;
        }
        let mut next_row = 0usize;
        let start = buf.len();
        let written = document.query(&xml.run, &query, &mut buf, |sink| {
            let Some(hit_list) = hit_list else {
                return Ok(());
            };
            for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
                let Some(first) = hsp_list.hsps.first() else {
                    continue;
                };
                let s_idx = first.s_idx as usize;
                let subject = subjects.get(s_idx, program)?;
                let hit = XmlHit::of(subject, xml.subjects.length(s_idx));
                // NCBI's messages about the taxonomy database: opening and the volume's lookup
                // at the first `Hit` that needs them, the lookup of the `Hit` every time.
                if let Some(warnings) = warnings.as_deref_mut() {
                    taxdb_printer.at_names(s_idx, subject, &mut *warnings.sink)?;
                }
                sink.hit(&hit, |hsps| {
                    for hsp in &hsp_list.hsps {
                        let row = hit_rows
                            .get(next_row)
                            .context("BLASTP report rows do not cover the hit list")?;
                        next_row += 1;
                        let (Some(qrow), Some(srow)) = (&row.query_seq, &row.subject_seq) else {
                            bail!("unimplemented BLASTP report alignment for final HSP");
                        };
                        hsps.hsp(&XmlHsp {
                            score: hsp.raw_score,
                            bit_score: hsp.bit_score,
                            evalue: seqalign_evalue(hsp.e_value),
                            query_from: hsp.q_start as i32,
                            query_to: hsp.q_end as i32,
                            hit_from: hsp.s_start as i32,
                            hit_to: hsp.s_end as i32,
                            query_frame: 0,
                            hit_frame: 0,
                            qrow: qrow.as_bytes(),
                            srow: srow.as_bytes(),
                            ops: hsp.gap_info.as_deref(),
                            query_first: hsp.q_start as i32 - 1,
                            query_step: 1,
                        })?;
                    }
                    Ok(())
                })?;
            }
            Ok(())
        });
        if let Err(error) = written {
            buf.truncate(start);
            return Err(error);
        }
        writer.write_all(&buf)?;
        buf.clear();
    }
    if xml.epilog {
        document.finish(&mut buf)?;
        writer.write_all(&buf)?;
    }
    Ok(())
}
