//! BLASTX's BLAST XML (outfmt 5): the `CBlastOutput` document of one query batch's final hit
//! lists (`report::blast_xml`), one query at a time.
//!
//! NCBI makes the Iteration of each query from its Seq-aligns, in the order the Seq-annot
//! formats show them (`seqalign`): subject by subject in hit-list order, each subject's HSPs
//! sorted by e-value (`render_queries` sorts them). What `s_SeqAlignSetToXMLHsps` reads of an
//! HSP is its Seq-align (coordinates in the records, score, e-value and bit score) and the two
//! rows of the alignment, from which it counts identities and positives again; nothing of the
//! engine's counts and nothing of the outfmt 0 midline. The rows are the ones outfmt 0 shows
//! (`report::alignment`), the masks the context's nucleotide intervals
//! (`QueryContext::dna_masks`, after the core's protein-to-DNA conversion).
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
//! (`PrintOneResultSet` is called once per query's result set; the accumulated lists hold that
//! one query.)

use std::io::Write;

use anyhow::{Context, Result};

use super::args::ResolvedOptions;
use super::input::FastaRecord;
use super::report::{alignment, has_seqalign, reported, ShownSubjects};
use super::runtime::BatchResults;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names::{self, QueryText};
use crate::report::blast_xml::subjects::SubjectNames;
use crate::report::blast_xml::{
    filter, BlastXmlWriter, QueryMask, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
    REFERENCE_GAPPED_BLAST,
};
use crate::report::query_warnings::INVALID_QUERY_MESSAGE;
use crate::report::serial::seqalign::seqalign_evalue;

const PROGRAM: &str = "BLASTX";

/// LOSAT's refusal of a title or name that NCBI's XML stream cannot write (not an error of
/// NCBI's engine: the run ends with exit status 1 and the message, `native::report_error`
/// passes it on).
#[derive(Debug)]
pub(super) struct XmlRefusal(String);

impl std::fmt::Display for XmlRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for XmlRefusal {}

fn refuse(refusal: names::Refusal, what: &str) -> anyhow::Error {
    XmlRefusal(refusal.error(PROGRAM, what).to_string()).into()
}

/// `s_GetTranslationFrame` (blastxml_format.cpp:85-98) on the 1-based endpoints of the HSP's
/// query interval (`start` the lower, `end` the higher) in a record of `length` letters.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:90-98
/// ```c++
///     int frame;
///
///     if (plus_strand) {
///         frame = (start - 1) % 3 + 1;
///     } else {
///         frame = -((seq_length - end) % 3 + 1);
///     }
///
///     return frame;
/// ```
fn translation_frame(plus_strand: bool, start: i32, end: i32, length: i32) -> i32 {
    if plus_strand {
        (start - 1) % 3 + 1
    } else {
        -((length - end) % 3 + 1)
    }
}

/// The masked regions of a query as NCBI hands them to the XML writer
/// (`Blast_GetSeqLocInfoVector`, blast_aux.cpp:930-959): per context, in context order, the
/// context's nucleotide intervals other than an empty one and the whole query, with the
/// context's frame. LOSAT's BLASTX has no `-query_loc`, so the interval is the whole record.
fn query_masks(batch: &BatchResults, qidx: usize, length: usize) -> Vec<QueryMask> {
    let mut masks = Vec::new();
    for context in &batch.prepared.contexts[qidx * 6..qidx * 6 + 6] {
        for &(from, to) in &context.dna_masks {
            if from > to || (from == 0 && to == length as i32 - 1) {
                continue;
            }
            masks.push(QueryMask {
                from,
                to,
                frame: i32::from(context.frame),
            });
        }
    }
    masks
}

/// The `Statistics` of a query: the first valid context's length adjustment and search
/// space, with the gapped Karlin block (the table value; the ideal ungapped block of the
/// context for `-ungapped`); the unsearched values for a query without a valid context.
fn query_stats(
    batch: &BatchResults,
    qidx: usize,
    gapped: bool,
    gapped_karlin: &crate::stats::tables::KarlinParams,
) -> XmlStats {
    let Some(parameters) = batch.parameters[qidx * 6..qidx * 6 + 6]
        .iter()
        .find(|parameters| parameters.valid)
    else {
        return XmlStats::NONE;
    };
    let karlin = if gapped {
        gapped_karlin
    } else {
        &parameters.ungapped
    };
    XmlStats {
        hsp_len: parameters.length_adjustment,
        eff_space: parameters.search_space,
        kappa: karlin.k,
        lambda: karlin.lambda,
        entropy: karlin.h,
    }
}

/// The messages of a query's result (the warnings the other formats write to standard error,
/// `report::query_warning`): the Karlin-Altschul message of a batch NCBI did not search, and
/// the message of a record without letters.
fn query_messages(query: &FastaRecord, search_skipped: bool) -> Vec<String> {
    let mut messages = Vec::new();
    if search_skipped {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    if query.sequence.is_empty() {
        messages.push("Sequence contains no data".to_string());
    }
    messages
}

/// The names of a `Hit`, for `-subject` (`Subject_N` and the generated definition line) or a
/// database's sequence (`gnl|BL_ORD_ID|n`, the stored title as stored).
fn subject_names(
    shown: &ShownSubjects<'_>,
    subject: &FastaRecord,
    oid: usize,
) -> Result<SubjectNames> {
    match shown {
        ShownSubjects::Fasta { .. } => {
            let (id, def) = names::xml_subject_names(
                subject.internal_id.as_bytes(),
                subject.title.as_bytes(),
                true,
            )
            .map_err(|refusal| refuse(refusal, "a subject title"))?;
            Ok(SubjectNames {
                id,
                def,
                accession: subject.internal_id.as_bytes().to_vec(),
            })
        }
        ShownSubjects::Database { reported, .. } => {
            let ordinal = reported.ordinal(oid);
            let def = names::xml_db_def(&reported.stored_titles()[oid])
                .map_err(|refusal| refuse(refusal, "a database title"))?;
            Ok(SubjectNames {
                id: format!("gnl|BL_ORD_ID|{ordinal}").into_bytes(),
                def,
                accession: ordinal.to_string().into_bytes(),
            })
        }
    }
}

/// Refuses a `-db` name that NCBI's XML stream cannot write (a byte outside US-ASCII); the
/// run checks it before the search, as the Seq-annot formats check theirs.
pub(super) fn check_database_name(name: &str) -> Result<()> {
    names::check_text(name.as_bytes()).map_err(|refusal| refuse(refusal, "a database name"))
}

/// Writes the documents of every query of `batches` (one query batch) to `writer`: the
/// document's header before the first query of the run (`queries_before` is the number of
/// queries formatted before this call), each query's Iteration, and nothing else; the end of
/// the document is `write_epilog`'s. Everything NCBI cannot write is refused before this
/// call's first byte: the query titles and messages, the titles of the subjects with hits.
pub(super) fn write_blastx_xml(
    writer: &mut impl Write,
    options: &ResolvedOptions,
    queries: &[FastaRecord],
    subjects: &[FastaRecord],
    shown: &ShownSubjects<'_>,
    batches: &[BatchResults],
    queries_before: usize,
) -> Result<()> {
    use crate::config::{ProteinScoringSpec, ScoringMatrix};
    let gapped_karlin = crate::stats::tables::lookup_protein_params(&ProteinScoringSpec {
        matrix: ScoringMatrix::Blosum62,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
    });
    let positives = FormattingMatrix::new(&options.matrix)
        .with_context(|| format!("no formatting matrix for {}", options.matrix))?;
    let (db_name, db_num, db_len) = match shown {
        ShownSubjects::Fasta { .. } => (String::new(), 0, 0),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.to_string(),
            i32::try_from(pairwise.num_sequences)?,
            i64::try_from(pairwise.total_letters)?,
        ),
    };
    check_database_name(&db_name)?;
    let filter = filter::filter_string(
        None,
        options
            .seg
            .enabled
            .then_some((options.seg.window, options.seg.locut, options.seg.hicut)),
        None,
        options.soft_masking,
    );
    let run = XmlRun {
        program: "blastx",
        version: "BLASTX 2.17.0+",
        reference: REFERENCE_GAPPED_BLAST,
        db_name: db_name.as_bytes(),
        matrix: Some(&options.matrix),
        expect: options.evalue,
        gap_open: options.gap_open,
        gap_extend: options.gap_extend,
        reward_penalty: None,
        filter: &filter,
        db_num,
        db_len,
        positives: Some(&positives),
        mask_char: b'X',
    };
    // The definition line and message of every query, and the names of the subjects with hits,
    // made (and refused) before the first byte.
    let mut texts: Vec<QueryText> = Vec::new();
    let mut subject_cache: Vec<Option<SubjectNames>> = Vec::new();
    subject_cache.resize_with(subjects.len(), || None);
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let has_hits = lists.iter().any(|list| reported(list, options.gapped));
            texts.push(
                names::query_text(
                    PROGRAM,
                    &query.internal_id,
                    query.title.as_bytes(),
                    false,
                    &query_messages(query, batch.search_skipped),
                    has_hits,
                )
                .map_err(|error| anyhow::Error::new(XmlRefusal(error.to_string())))?,
            );
            for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                if subject_cache[list.oid].is_none() {
                    subject_cache[list.oid] =
                        Some(subject_names(shown, &subjects[list.oid], list.oid)?);
                }
            }
        }
    }
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    let mut document = BlastXmlWriter::new(queries_before);
    let mut buf: Vec<u8> = Vec::new();
    let mut texts = texts.iter();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let text = texts.next().context("one text per query")?;
            let length = query.sequence.len();
            let masks = query_masks(batch, qidx, length);
            let xml_query = XmlQuery {
                id: query.internal_id.as_bytes(),
                def: &text.def,
                len: i32::try_from(length)?,
                stats: query_stats(batch, qidx, options.gapped, &gapped_karlin),
                message: text.message.as_deref(),
            };
            let start = buf.len();
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                    let subject = &subjects[list.oid];
                    let names = subject_cache[list.oid]
                        .as_ref()
                        .context("names of a subject with hits")?;
                    let hit = XmlHit {
                        id: &names.id,
                        def: &names.def,
                        accession: &names.accession,
                        len: i32::try_from(subject.sequence.len())?,
                    };
                    sink.hit(&hit, |hsps| {
                        for hsp in list
                            .hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                        {
                            let rows = alignment(batch, hsp, subject, query, &code)?;
                            // The Seq-align's endpoints, lower first: the minus frame's
                            // outfmt 6 order is not the XML's.
                            let (query_from, query_to) =
                                (rows.qstart.min(rows.qend), rows.qstart.max(rows.qend));
                            let plus = rows.frame > 0;
                            hsps.hsp(&XmlHsp {
                                score: hsp.hsp.score,
                                bit_score: hsp.bit_score,
                                evalue: seqalign_evalue(hsp.evalue),
                                query_from,
                                query_to,
                                hit_from: rows.sstart,
                                hit_to: rows.send,
                                query_frame: translation_frame(
                                    plus,
                                    query_from,
                                    query_to,
                                    length as i32,
                                ),
                                hit_frame: 0,
                                qrow: rows.qseq.as_bytes(),
                                srow: rows.sseq.as_bytes(),
                                ops: Some(&hsp.edit_script),
                                // The first residue of the row: the lowest nucleotide of the
                                // interval on the plus strand, the highest on the minus
                                // strand, where the row runs backwards.
                                query_first: if plus { query_from - 1 } else { query_to - 1 },
                                query_step: if plus { 3 } else { -3 },
                                masks: &masks,
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
    }
    Ok(())
}

/// The end of the document (`PrintEpilog`): its trailer when `queries` queries were formatted,
/// else the epilog's end of line alone.
pub(super) fn write_xml_epilog(writer: &mut impl Write, queries: usize) -> Result<()> {
    let mut buf = Vec::new();
    BlastXmlWriter::new(queries).finish(&mut buf);
    writer.write_all(&buf)?;
    Ok(())
}
