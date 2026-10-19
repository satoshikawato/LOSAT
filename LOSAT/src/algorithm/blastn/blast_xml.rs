//! BLASTN's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit lists
//! (`report::blast_xml`), one query at a time, for every task (megablast, blastn,
//! dc-megablast, blastn-short).
//!
//! What NCBI's `s_SeqAlignSetToXMLHsps` reads of an HSP is its Seq-align: the coordinates in the
//! records (the query always ascending; the subject descending for a minus-strand HSP, whose
//! Dense-seg NCBI reverses so that the query is shown on its plus strand), the score, e-value
//! and bit score, and the two rows from which it counts identities again. The masks never apply
//! to a nucleotide query: blastn's masks have no frame (`eFrameNotSet`, api/blast_aux.cpp:
//! 892) and the HSP's query frame is 1, so `s_MaskQuerySeq` skips them all.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:249-257
//! ```c++
//!         if (!kTranslated && query_is_na && subject_is_na) {
//!             q_frame = s_frame = 1;
//!             // For reverse strand alignment, set subject frame to -1 and
//!             // swap start and end coordinates.
//!             if (eNa_strand_minus == kAlign.GetSeqStrand(0)){
//!                 s_frame = -1;
//!                 int tmp = s_start;
//!                 s_start = s_end;
//!                 s_end = tmp;
//!             }
//! ```

use std::io::Write;

use anyhow::{Context, Result};

use super::hsp::BlastnHitList;
use super::pairwise::xml_rows;
use super::x_subject_store::SubjectLetters;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    BlastXmlWriter, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats, XmlSubjects,
};
use crate::report::query_warnings::QueryWarnings;

/// What the BLASTN document needs besides the final hit lists.
pub(crate) struct BlastnXml<'a> {
    pub run: XmlRun<'a>,
    /// The queries of this call and their searched lengths.
    pub queries: &'a [FastaRecord],
    pub searched_lengths: &'a [usize],
    /// By query: the statistics (`XmlStats::NONE` for a query without a valid context).
    pub stats: Vec<XmlStats>,
    /// By query: the warning messages of its result (a query without letters, an invalid query).
    pub messages: Vec<Vec<String>>,
    pub subjects: XmlSubjects<'a>,
    /// The letters of the subjects, which the rows read.
    pub letters: SubjectLetters<'a>,
    /// The number of queries formatted before this call, and whether the document ends with it.
    pub queries_before: usize,
    pub epilog: bool,
}

/// Writes the document of every query of `hit_lists` (indexed by query) to `writer`, each
/// query's warnings before it, with the document's header before the first query and its end
/// after the last when `xml.epilog`.
pub(crate) fn write_blastn_xml(
    xml: &BlastnXml<'_>,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut dyn Write,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    const PROGRAM: &str = "BLASTN";
    names::check_text(xml.run.db_name)
        .map_err(|refusal| refusal.error(PROGRAM, "a database name"))?;
    let has_hits = |q_idx: usize| {
        hit_lists
            .get(q_idx)
            .and_then(Option::as_ref)
            .is_some_and(|list| list.hsplist_array.iter().any(|list| !list.hsps.is_empty()))
    };
    let texts = names::query_texts(PROGRAM, xml.queries, false, &xml.messages, has_hits)?;
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            if let Some(hsp) = hsp_list.hsps.first() {
                xml.subjects.names(hsp.s_idx as usize, PROGRAM)?;
            }
        }
    }
    let mut document = BlastXmlWriter::new(xml.queries_before);
    let mut buf: Vec<u8> = Vec::new();
    for (q_idx, record) in xml.queries.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let text = &texts[q_idx];
        let query = XmlQuery {
            id: record.local_id.as_bytes(),
            def: &text.def,
            len: xml.searched_lengths[q_idx] as i32,
            stats: xml.stats[q_idx],
            message: text.message.as_deref(),
        };
        let hit_list = hit_lists.get(q_idx).and_then(Option::as_ref);
        let start = buf.len();
        let written = document.query(&xml.run, &query, &mut buf, |sink| {
            let Some(hit_list) = hit_list else {
                return Ok(());
            };
            for hsp_list in &hit_list.hsplist_array {
                let Some(first) = hsp_list.hsps.first() else {
                    continue;
                };
                let s_idx = first.s_idx as usize;
                let subject = xml.subjects.names(s_idx, PROGRAM)?;
                let hit = XmlHit {
                    id: &subject.id,
                    def: &subject.def,
                    accession: &subject.accession,
                    len: xml.subjects.length(s_idx),
                };
                sink.hit(&hit, |hsps| {
                    for hsp in &hsp_list.hsps {
                        let (qrow, srow) = xml_rows(hsp, xml.queries, xml.letters)
                            .context("BLASTN report rows of a final HSP")?;
                        let minus = hsp.query_frame < 0;
                        hsps.hsp(&XmlHsp {
                            score: hsp.raw_score,
                            bit_score: hsp.bit_score,
                            evalue: seqalign_evalue(hsp.e_value),
                            query_from: hsp.q_start as i32,
                            query_to: hsp.q_end as i32,
                            hit_from: hsp.s_start as i32,
                            hit_to: hsp.s_end as i32,
                            query_frame: 1,
                            hit_frame: if minus { -1 } else { 1 },
                            qrow: &qrow,
                            srow: &srow,
                            ops: hsp.gap_info.as_deref(),
                            query_first: hsp.q_start as i32 - 1,
                            masks: &[],
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
        document.finish(&mut buf);
        writer.write_all(&buf)?;
    }
    Ok(())
}

/// The e-value of the Seq-align: 0 below `SMALLEST_EVALUE`
/// (api/blast_seqalign.cpp:1186).
fn seqalign_evalue(evalue: f64) -> f64 {
    if evalue < 1.0e-180 {
        0.0
    } else {
        evalue
    }
}
