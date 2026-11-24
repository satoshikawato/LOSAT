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
//! (`PrintOneResultSet` is called once per query's result set, the accumulated lists hold that
//! one query.)

use std::io::Write;

use anyhow::{bail, Context, Result};

use super::blast_engine::query_rows::QueryRows;
use super::hsp::BlastpHitList;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    hit_ranges, written_queries, QueryMask, StructuredXml, SubjectNameCache, XmlHit, XmlHsp,
    XmlQuery, XmlRun, XmlStats, XmlSubjects,
};
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::QueryWarnings;
use crate::report::record::QueryNames;
use crate::report::serial::seqalign::{read_num_ident, seqalign_evalue};

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
    /// The app's program type (`BlastpApp::program_type`, `RunRecord::program_type`): psiblast's
    /// Seq-align of an HSP of 0 identities has no `num_ident` (`serial::seqalign::read_num_ident`).
    pub program_type: Option<crate::report::record::ProgramType>,
    /// Whether the document ends with this call.
    pub epilog: bool,
    /// Whether the run stops after this call: no epilog and no later report call (the block
    /// search's next block) follow it (`written_queries`).
    pub stops: bool,
    /// Every structured format of the call is XML2 or JSON (outfmt 5, which writes every report,
    /// is not among them; the library takes several formats).
    pub only_xml2: bool,
    /// psiblast's PSSM query of the molecule `other` (`PsiSearch::query_other_coding`): the
    /// first `Hsp`'s rows throw (`report::tabular::query_sequence_exception`).
    pub query_other_coding: Option<i32>,
}

/// Writes the document of every query of `hit_lists` (indexed by query) through `document`
/// (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before it,
/// with the document's header before the first query and its end after the last when
/// `xml.epilog`. `rows` gives the rows of each query's HSPs, in the lists' order.
pub(crate) fn write_blastp_xml<W: StructuredXml>(
    xml: &BlastpXml<'_>,
    document: &mut W,
    hit_lists: &[Option<BlastpHitList>],
    rows: &mut QueryRows<'_>,
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let program = xml.program;
    // A title NCBI crashes on (`names::Refusal`) is refused before the first byte of the
    // document: for the reports that the document writes, each query's definition line and the
    // names of the subjects with hits. A run that stops never writes the
    // `-subject` report of its last query, which NCBI holds back, unless another format of the
    // call (outfmt 5) writes it (`written_queries`).
    let reports = written_queries(
        hit_lists.len(),
        xml.stops,
        xml.only_xml2 && document.holds_back(),
    );
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
        &xml.queries.records[..reports],
        true,
        &xml.messages,
        has_hits,
    )?;
    // The names of the subjects with hits in the written reports, made once each; a title NCBI
    // crashes on is refused here, before the first byte.
    let mut subjects = SubjectNameCache::new(&xml.subjects);
    for hit_list in hit_lists.iter().take(reports).flatten() {
        for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
            if let Some(hsp) = hsp_list.hsps.first() {
                subjects.get(hsp.s_idx as usize, program)?;
            }
        }
    }
    let mut taxdb_printer = xml.taxdb_notes.printer();
    let mut buf: Vec<u8> = Vec::new();
    // The document's start (`StructuredXml::start`), before the first query's warnings, after the
    // refusal of the database counts that the written reports would print.
    document.start((reports > 0).then_some(&xml.run), &mut buf)?;
    writer.write_all(&buf)?;
    buf.clear();
    for (q_idx, hit_list) in hit_lists.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        if q_idx == reports {
            // The held report that is never written (`written_queries`).
            document.unwritten(&mut buf)?;
            crate::report::query_warnings::write_document_posts(
                &buf,
                document.take_posts(),
                writer,
                warnings.as_deref_mut(),
            )?;
            buf.clear();
            break;
        }
        let record = &xml.queries.records[q_idx];
        let text = &texts[q_idx];
        let query = XmlQuery {
            id: text.query_id(record.local_id.as_bytes()),
            title: &text.title,
            len: xml.queries.searched_lengths[q_idx] as i32,
            stats: xml.stats[q_idx],
            message: text.message.as_deref(),
            searched: xml.searched[q_idx],
            masks: &xml.masks[q_idx],
        };
        let hit_rows: &[PairwiseHit] = match hit_list {
            Some(_) => rows.query(q_idx)?,
            None => &[],
        };
        let mut next_row = 0usize;
        let start = buf.len();
        // What NCBI posts where it loads the subject of each `Hit` (`HitPosts`): the texts of the
        // subjects that the query's report reads one by one (`QueryWarnings::subject_read`) and
        // the taxonomy database's messages.
        let mut posts = crate::report::query_warnings::HitPosts::new(
            warnings
                .as_deref()
                .is_some_and(|warnings| warnings.reads_with_texts(q_idx)),
        );
        let posting = warnings.is_some();
        let written = document.query(&xml.run, &query, &mut buf, |sink| {
            let Some(hit_list) = hit_list else {
                return Ok(());
            };
            let lists =
                &hit_list.hsplist_array[..hit_list.hsplist_count.min(hit_list.hsplist_array.len())];
            // A copy of one record that follows another is shown in its `Hit`.
            let ranges = hit_ranges(
                lists,
                |hsp_list| hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize),
                |previous, s_idx| xml.subjects.joins(previous, s_idx),
            );
            for range in ranges {
                let s_idx = lists[range.start].hsps[0].s_idx as usize;
                let subject = subjects.get(s_idx, program)?.names()?;
                let hit = XmlHit::of(subject, xml.subjects.length(s_idx));
                // NCBI's messages about the taxonomy database: opening and the volume's lookup
                // at the first `Hit` that needs them, the lookup of the `Hit` every time.
                posts.hit(sink.position(), &hit, s_idx as u32, |texts| {
                    if posting {
                        taxdb_printer.at_names(s_idx, subject, texts)
                    } else {
                        Ok(())
                    }
                })?;
                sink.hit(&hit, |hsps| {
                    for hsp in lists[range].iter().flat_map(|hsp_list| &hsp_list.hsps) {
                        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:280-300
                        // The rows read the query's letters, which throw for a query of the
                        // molecule `other`; `s_BlastXMLAddIteration` catches it.
                        if let Some(coding) = xml.query_other_coding {
                            return Err(crate::report::blast_xml::CaughtException(
                                crate::report::tabular::query_sequence_exception(coding),
                            )
                            .into());
                        }
                        let row = hit_rows
                            .get(next_row)
                            .context("BLASTP report rows do not cover the hit list")?;
                        next_row += 1;
                        let identities =
                            read_num_ident(row.seqalign_num_ident as i64, xml.program_type);
                        if hsps.unaligned([XmlHsp::unaligned(
                            hsp.raw_score,
                            hsp.bit_score,
                            seqalign_evalue(hsp.e_value),
                            identities,
                        )])? {
                            continue;
                        }
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
                            alignment_identities: identities,
                        })?;
                    }
                    Ok(())
                })?;
            }
            Ok(())
        });
        // A query whose report raised an exception that NCBI's report function catches
        // (`CaughtException`: the `Hit` of a subject whose lookup throws, `HitNames`): the XML2 and
        // JSON documents leave its report out (`Xml2Document::report`), and the post of the
        // exception follows the texts of the loads before it (`HitPosts::caught`). Another failed
        // query leaves its partial bytes in `buf`, which the error return drops.
        match crate::report::blast_xml::writer::caught_post(written, &*document, xml.run.program) {
            Ok(Some(post)) => posts.caught(post),
            Ok(None) => {}
            Err(error) => {
                buf.truncate(start);
                return Err(error);
            }
        }
        posts.document(document.take_posts());
        crate::report::query_warnings::write_hit_posts(
            &buf,
            &posts,
            q_idx,
            writer,
            warnings.as_deref_mut(),
        )?;
        buf.clear();
    }
    if xml.epilog {
        document.finish(&mut buf)?;
        crate::report::query_warnings::write_document_posts(
            &buf,
            document.take_posts(),
            writer,
            warnings,
        )?;
    }
    Ok(())
}
