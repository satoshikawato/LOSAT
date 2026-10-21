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
//!
//! Every refusal (`XmlRefusal`: a byte NCBI's XML stream cannot write in a title) comes before
//! the first byte of the document, though the query file is read batch by batch:
//!
//! - `prescan_queries` reads the query file once more before the search (only with outfmt 5)
//!   and refuses the first query whose title the writer would refuse;
//! - when some `-subject` title would be refused if its subject had a hit
//!   (`subjects_may_be_refused`), or some query's message would be (a title that is refused only
//!   in the message of a query NCBI did not search), the run holds the document in memory
//!   (`HeldOutput`) and writes it at the epilog; a refusal found while the batches are reported
//!   drops it, with nothing written.

use std::cell::RefCell;
use std::collections::HashMap;
use std::io::Write;
use std::rc::Rc;

use anyhow::{ensure, Context, Result};

use super::args::ResolvedOptions;
use super::input::FastaRecord;
use super::report::{
    display_query_residue, has_seqalign, query_endpoints, reported, ShownSubjects,
};
use super::results::Hsp;
use super::runtime::BatchResults;
use crate::common::GapEditOp;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names::{self, QueryText};
use crate::report::blast_xml::subjects::SubjectNames;
use crate::report::blast_xml::{
    filter, BlastXmlWriter, QueryMask, StructuredXml, XmlHit, XmlHsp, XmlQuery, XmlRun, XmlStats,
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
        ShownSubjects::Fasta { .. } => SubjectNames::of_record(
            subject.internal_id.as_bytes(),
            subject.title.as_bytes(),
            true,
        )
        .map_err(|refusal| refuse(refusal, "a subject title")),
        ShownSubjects::Database { reported, .. } => SubjectNames::of_database_sequence(
            reported.ordinal(oid),
            &reported.stored_titles()[oid],
            // XML2 alone shows the taxid (outfmt 15 and 16 of BLASTX: OF5b).
            0,
            Default::default(),
        )
        .map_err(|refusal| refuse(refusal, "a database title")),
    }
}

/// The rows of an HSP's alignment into `qrow` and `srow` (cleared first): the letters
/// `report::alignment` shows in outfmt 0 (the query's codons translated, the subject's
/// residues, `-` for a gap), without the rest of what it builds (the masked row, the midline
/// counts), which the XML counts again from the rows itself. Returns the HSP's endpoints
/// `(query_from, query_to)` as the Seq-align has them, lower first.
fn xml_rows(
    hsp: &Hsp,
    subject: &FastaRecord,
    query: &FastaRecord,
    code: &crate::utils::genetic_code::GeneticCode,
    qrow: &mut Vec<u8>,
    srow: &mut Vec<u8>,
) -> Result<(i32, i32)> {
    let raw = &hsp.hsp;
    qrow.clear();
    srow.clear();
    let (mut qi, mut si) = (usize::try_from(raw.q_start)?, usize::try_from(raw.s_start)?);
    let whole = [GapEditOp::Sub(u32::try_from(raw.q_end - raw.q_start)?)];
    let script: &[GapEditOp] = if hsp.edit_script.is_empty() {
        &whole
    } else {
        &hsp.edit_script
    };
    for &op in script {
        let (count, consume_q, consume_s) = match op {
            GapEditOp::Sub(n) => (n, true, true),
            GapEditOp::Ins(n) => (n, true, false),
            GapEditOp::Del(n) => (n, false, true),
        };
        for _ in 0..count {
            qrow.push(if consume_q {
                display_query_residue(&query.sequence, raw.frame, qi, code)
            } else {
                b'-'
            });
            srow.push(if consume_s {
                subject.sequence[si]
            } else {
                b'-'
            });
            qi += usize::from(consume_q);
            si += usize::from(consume_s);
        }
    }
    ensure!(
        qi == raw.q_end as usize && si == raw.s_end as usize,
        "BLASTX edit script does not reach final endpoints"
    );
    let (qstart, qend) = query_endpoints(query.sequence.len(), raw.frame, raw.q_start, raw.q_end);
    Ok((qstart.min(qend), qstart.max(qend)))
}

/// The output of a run with outfmt 5: it passes the bytes on, or, once `hold` is called, keeps
/// them in memory until `release` (the epilog), so that a refusal found while the batches are
/// reported leaves nothing written. A run that ends with any other error writes what it held
/// (the bytes of the batches before it, as the run without a hold does); `discard` is for the
/// refusal alone.
#[derive(Clone)]
pub(super) struct HeldOutput(Rc<RefCell<HeldState>>);

struct HeldState {
    inner: Box<dyn Write>,
    holding: bool,
    held: Vec<u8>,
    discarded: bool,
}

impl HeldOutput {
    pub(super) fn new(inner: Box<dyn Write>) -> Self {
        HeldOutput(Rc::new(RefCell::new(HeldState {
            inner,
            holding: false,
            held: Vec::new(),
            discarded: false,
        })))
    }

    /// Keeps everything written from now on in memory.
    pub(super) fn hold(&self) {
        self.0.borrow_mut().holding = true;
    }

    /// Drops what is held, and everything written after (the run ends with a refusal).
    pub(super) fn discard(&self) {
        let mut state = self.0.borrow_mut();
        state.discarded = true;
        state.held.clear();
    }

    /// Writes what is held and passes the bytes on again.
    pub(super) fn release(&self) -> std::io::Result<()> {
        let mut state = self.0.borrow_mut();
        state.holding = false;
        let held = std::mem::take(&mut state.held);
        state.inner.write_all(&held)?;
        state.inner.flush()
    }
}

impl Write for HeldOutput {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        let mut state = self.0.borrow_mut();
        if state.discarded {
            Ok(bytes.len())
        } else if state.holding {
            state.held.extend_from_slice(bytes);
            Ok(bytes.len())
        } else {
            state.inner.write(bytes)
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        let mut state = self.0.borrow_mut();
        if state.holding || state.discarded {
            Ok(())
        } else {
            state.inner.flush()
        }
    }
}

impl Drop for HeldState {
    fn drop(&mut self) {
        if !self.discarded && !self.held.is_empty() {
            let _ = self.inner.write_all(&self.held);
            let _ = self.inner.flush();
        }
    }
}

/// The `QueryText` of a query, refused as the writer refuses it. `search_skipped`: NCBI did not
/// search the query's batch (its message then carries the query's id, made of its title).
fn query_text(query: &FastaRecord, search_skipped: bool, has_hits: bool) -> Result<QueryText> {
    names::query_text(
        PROGRAM,
        &query.internal_id,
        query.title.as_bytes(),
        false,
        &query_messages(query, search_skipped),
        has_hits,
    )
    .map_err(|error| anyhow::Error::new(XmlRefusal(error.to_string())))
}

/// Reads the query file once before the search and refuses the first query the writer would
/// refuse (a title with a byte NCBI's XML stream cannot write), so that no refusal can follow
/// a written byte. Returns whether some query would be refused if NCBI had not searched its
/// batch (a title that only its message cannot carry); the run then holds the document.
///
/// The reader is the search's own (same records and titles). A record it refuses for every
/// format (a title that is not UTF-8) ends the scan: the run reads the file again and fails
/// there as it always did.
pub(super) fn prescan_queries(file: std::fs::File, lowercase: bool) -> Result<bool> {
    let mut stream = super::input::FastaStream::from_file(file);
    let mut may_be_refused = false;
    let scanned = super::input::parse_fasta_stream_with_warnings(
        &mut stream,
        false,
        lowercase,
        super::query_setup::BATCH_SIZE,
        &mut |batch| {
            for query in batch {
                query_text(query, false, true)?;
                may_be_refused |= query_text(query, true, true).is_err();
            }
            Ok(())
        },
        &mut |_| Ok(()),
    );
    match scanned {
        Err(error) if error.is::<XmlRefusal>() => Err(error),
        _ => Ok(may_be_refused),
    }
}

/// Whether the title of some `-subject` would be refused if the subject had a hit.
pub(super) fn subjects_may_be_refused(subjects: &[FastaRecord]) -> bool {
    subjects.iter().any(|subject| {
        names::xml_subject_names(
            subject.internal_id.as_bytes(),
            subject.title.as_bytes(),
            true,
        )
        .is_err()
    })
}

/// Refuses a `-db` name that NCBI's XML stream cannot write (a byte outside US-ASCII); the
/// run checks it before the search, as the Seq-annot formats check theirs.
pub(super) fn check_database_name(name: &str) -> Result<()> {
    names::check_database_name(PROGRAM, name).map_err(|error| XmlRefusal(error.to_string()).into())
}

/// Writes the documents of every query of `batches` (one query batch) to `writer`: the
/// document's header before the first query of the run (`queries_before` is the number of
/// queries formatted before this call), each query's Iteration, and nothing else; the end of
/// the document is `write_epilog`'s. Everything NCBI cannot write is refused before this
/// call's first byte: the query titles and messages, the titles of the subjects with hits
/// (and, with a hold, the whole run's bytes stay in memory, `HeldOutput`).
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
    let mut subject_cache: HashMap<usize, SubjectNames> = HashMap::new();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let has_hits = lists.iter().any(|list| reported(list, options.gapped));
            texts.push(query_text(query, batch.search_skipped, has_hits)?);
            for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                if let std::collections::hash_map::Entry::Vacant(slot) =
                    subject_cache.entry(list.oid)
                {
                    slot.insert(subject_names(shown, &subjects[list.oid], list.oid)?);
                }
            }
        }
    }
    let code = crate::utils::genetic_code::GeneticCode::from_id(options.query_gencode);
    let mut document = BlastXmlWriter::new(queries_before);
    let mut buf: Vec<u8> = Vec::new();
    // The rows of one HSP at a time, reused.
    let (mut qrow, mut srow): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    let mut texts = texts.iter();
    for batch in batches {
        for (qidx, lists) in batch.queries.iter().enumerate() {
            let query = &queries[batch.query_ordinal + qidx];
            let text = texts.next().context("one text per query")?;
            let length = query.sequence.len();
            let masks = query_masks(batch, qidx, length);
            let xml_query = XmlQuery {
                id: query.internal_id.as_bytes(),
                title: &text.title,
                len: i32::try_from(length)?,
                stats: query_stats(batch, qidx, options.gapped, &gapped_karlin),
                message: text.message.as_deref(),
                // Read by XML2 alone (outfmt 15 and 16 of BLASTX: OF5b).
                searched: !batch.search_skipped,
                masks: &masks,
            };
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                for list in lists.iter().filter(|list| reported(list, options.gapped)) {
                    let subject = &subjects[list.oid];
                    let names = subject_cache
                        .get(&list.oid)
                        .context("names of a subject with hits")?;
                    let hit = XmlHit {
                        id: &names.id,
                        title: &names.title,
                        accession: &names.accession,
                        len: i32::try_from(subject.sequence.len())?,
                        taxid: names.taxid,
                        taxonomy: &names.taxonomy,
                    };
                    sink.hit(&hit, |hsps| {
                        for hsp in list
                            .hsps
                            .iter()
                            .filter(|hsp| has_seqalign(hsp, options.gapped))
                        {
                            // The Seq-align's endpoints, lower first: the minus frame's
                            // outfmt 6 order is not the XML's.
                            let (query_from, query_to) =
                                xml_rows(hsp, subject, query, &code, &mut qrow, &mut srow)?;
                            let plus = hsp.hsp.frame > 0;
                            hsps.hsp(&XmlHsp {
                                score: hsp.hsp.score,
                                bit_score: hsp.bit_score,
                                evalue: seqalign_evalue(hsp.evalue),
                                query_from,
                                query_to,
                                hit_from: hsp.hsp.s_start + 1,
                                hit_to: hsp.hsp.s_end,
                                query_frame: translation_frame(
                                    plus,
                                    query_from,
                                    query_to,
                                    length as i32,
                                ),
                                hit_frame: 0,
                                qrow: &qrow,
                                srow: &srow,
                                ops: Some(&hsp.edit_script),
                                // The first residue of the row: the lowest nucleotide of the
                                // interval on the plus strand, the highest on the minus
                                // strand, where the row runs backwards.
                                query_first: if plus { query_from - 1 } else { query_to - 1 },
                                query_step: if plus { 3 } else { -3 },
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
    }
    Ok(())
}

/// The end of the document (`PrintEpilog`): its trailer when `queries` queries were formatted,
/// else the epilog's end of line alone.
pub(super) fn write_xml_epilog(writer: &mut impl Write, queries: usize) -> Result<()> {
    let mut buf = Vec::new();
    BlastXmlWriter::new(queries).finish(&mut buf)?;
    writer.write_all(&buf)?;
    Ok(())
}
