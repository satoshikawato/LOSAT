//! BLASTN's BLAST XML (outfmt 5) and XML2 / JSON (outfmt 15 and 16): the `CBlastOutput` document
//! of its final hit lists (`report::blast_xml`, `report::blast_xml2`), one query at a time, for
//! every task (megablast, blastn, dc-megablast, blastn-short).
//!
//! What NCBI's `s_SeqAlignSetToXMLHsps` reads of an HSP is its Seq-align: the coordinates in the
//! records (the query always ascending; the subject descending for a minus-strand HSP, whose
//! Dense-seg NCBI reverses so that the query is shown on its plus strand), the score, e-value
//! and bit score, and the two rows from which it counts identities again. The masks never apply
//! to a nucleotide query in outfmt 5: blastn's masks have no frame (`eFrameNotSet`, api/
//! blast_aux.cpp:892) and the HSP's query frame is 1, so `s_MaskQuerySeq` skips them all. XML2
//! compares them with the frame 0 of an untranslated query (`blastxml2_format.cpp:211, :253`), so
//! its `qseq` is masked and the `Search` lists them (`BlastnXml::masks`, kept empty for outfmt 5).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:249-257

use anyhow::{Context, Result};

use super::hsp::BlastnHitList;
use super::pairwise::xml_rows_into;
use super::x_subject_store::SubjectLetters;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    hit_ranges, written_queries, QueryMask, StructuredXml, SubjectNameCache, XmlHit, XmlHsp,
    XmlQuery, XmlRun, XmlStats, XmlSubjects,
};
use crate::report::query_warnings::QueryWarnings;
use crate::report::serial::seqalign::{read_num_ident, seqalign_evalue};

/// The program that BLASTN's documents and refusals name.
const PROGRAM: &str = "BLASTN";

/// What the BLASTN document needs besides the final hit lists.
pub(crate) struct BlastnXml<'a> {
    pub run: XmlRun<'a>,
    /// The queries of this call and their searched lengths.
    pub queries: &'a [FastaRecord],
    /// The records whose letters the rows read (`pairwise::QueryRecords`).
    pub query_records: super::pairwise::QueryRecords<'a>,
    pub searched_lengths: &'a [usize],
    /// By query: the statistics (`XmlStats::NONE` for a query without a valid context).
    pub stats: Vec<XmlStats>,
    /// By query: the warning messages of its result (a query without letters, an invalid query).
    pub messages: Vec<Vec<String>>,
    /// By query, the masked regions in record coordinates (`QueryMask`, frame 0): the XML2 /
    /// JSON documents only; empty for outfmt 5, which does not show them.
    pub masks: Vec<Vec<QueryMask>>,
    /// By query: whether NCBI searched the query's batch (it has a valid context in it); empty
    /// when every query was. Only XML2 shows it.
    pub searched: Vec<bool>,
    pub subjects: XmlSubjects<'a>,
    /// NCBI's messages about the taxonomy database that come before the `Hit`s' own
    /// (XML2 and JSON with `-db`; nothing otherwise).
    pub taxdb_notes: crate::algorithm::common::database::TaxdbNotes,
    /// The letters of the subjects, which the rows read.
    pub letters: SubjectLetters<'a>,
    /// Whether the search is ungapped (`-ungapped`, Dense-diag Seq-aligns).
    pub ungapped: bool,
    /// Whether the document ends with this call.
    pub epilog: bool,
    /// Whether the run stops after this call without the epilog (the reports of the batches
    /// before a failed batch); a call without the epilog that a later call of the run follows
    /// (BLASTN's report stream, `x_rpt_stream`) does not stop it.
    pub stops: bool,
    /// Every structured format of the call is XML2 or JSON (outfmt 5, which writes every report,
    /// is not among them; the library takes several formats).
    pub only_xml2: bool,
    /// Under `-parse_deflines`, the scopes of the call's queries, which say how each result is
    /// shown (`query_batch::xml_show_of`); `None` otherwise.
    pub query_scopes: Option<crate::blastinput::query_batch::QueryScopesFrom<'a>>,
}

/// Writes the document of every query of `hit_lists` (indexed by query) through `document`
/// (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before it,
/// with the document's header before the first query and its end after the last when
/// `xml.epilog`.
pub(crate) fn write_blastn_xml<W: StructuredXml>(
    xml: &BlastnXml<'_>,
    document: &mut W,
    hit_lists: &[Option<BlastnHitList>],
    writer: &mut dyn crate::report::ostream::ReportWrite,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    // A title NCBI crashes on (`names::Refusal`) is refused where NCBI makes the report that
    // shows it, after only the bytes that NCBI has flushed by then (remaining-gaps DECISIONS D-SO-6,
    // Owner-delegated 2026-11-27, under the Owner's Q-X2c of 2026-11-16). A document that writes
    // each report at once (outfmt 5; XML2 and JSON with `-db`) makes the call's reports in one
    // step, so their refusals come before the call's first byte, after the reports of the calls
    // before. A `-subject` XML2 or JSON document (`holds_back`) makes a report only when the
    // result of another query arrives or in the epilog, and writes it at once: the refusal of a
    // report comes there, after the reports before it (`StructuredXml::hold_refused`), in one call
    // or across the calls of BLASTN's report stream alike; the report of the last query of a run
    // that stops is never made, nor refused (`written_queries`).
    let queries = xml.queries.len();
    let shows = |q_idx: usize| crate::blastinput::query_batch::xml_show_of(xml.query_scopes, q_idx);
    let holds_back = xml.only_xml2 && document.holds_back();
    let reports = written_queries(
        queries,
        xml.stops,
        holds_back,
        crate::blastinput::query_batch::held_start_of(xml.query_scopes, queries),
    );
    let has_hits = |q_idx: usize| {
        hit_lists
            .get(q_idx)
            .and_then(Option::as_ref)
            .is_some_and(|list| list.hsplist_array.iter().any(|list| !list.hsps.is_empty()))
    };
    let unshown = |q_idx: usize| {
        crate::blastinput::query_batch::names_unshown_in(
            xml.query_scopes,
            reports,
            holds_back,
            q_idx,
        )
    };
    // The texts of query `q_idx` (its refusal: `names::Refusal`).
    let text_of = |q_idx: usize| {
        names::query_texts(
            PROGRAM,
            &xml.queries[q_idx..=q_idx],
            false,
            &xml.messages[q_idx..=q_idx],
            |_| has_hits(q_idx),
            |_| unshown(q_idx),
        )
        .map(|mut text| text.pop().expect("the text of one query"))
    };
    let mut subjects = SubjectNameCache::new(&xml.subjects);
    let texts: Vec<Result<names::QueryText>> = if holds_back {
        // Each report's refusal waits for the place where NCBI makes the report (the loop below).
        (0..reports).map(text_of).collect()
    } else {
        // The names of the subjects with hits in the written reports are made once each; a title
        // NCBI crashes on is refused here, before the first byte.
        let texts = names::query_texts(
            PROGRAM,
            &xml.queries[..reports],
            false,
            &xml.messages[..reports],
            has_hits,
            unshown,
        )?;
        name_subjects(&mut subjects, &hit_lists[..reports.min(hit_lists.len())])?;
        texts.into_iter().map(Ok).collect()
    };
    let mut taxdb_printer = xml.taxdb_notes.printer();
    let (mut query_row, mut subject_row): (Vec<u8>, Vec<u8>) = (Vec::new(), Vec::new());
    // The call's output buffer, that of the call before where the writer keeps it.
    let mut buf: Vec<u8> = document.take_out();
    // The document's start (`StructuredXml::start`), before the first query's warnings, after the
    // refusal of the database counts that the written reports would print. A `-subject` XML2 or
    // JSON document's start waits for the output's first write or flush (`StartFirst`), as NCBI's
    // start waits in its stream's buffer, so that a refused first report leaves nothing.
    document.start((reports > 0).then_some(&xml.run), &mut buf)?;
    let mut writer = StartFirst {
        inner: writer,
        start: buf.split_off(0),
    };
    if !holds_back {
        writer.emit()?;
    }
    let writer = &mut writer;
    for q_idx in 0..xml.queries.len() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        let show = shows(q_idx);
        if q_idx == reports {
            // The held report that is never written (`written_queries`); the one held before it is
            // written, or fails (`CaughtException`: its post follows), unless this result joins
            // it (`XmlShow::joins_held`).
            let written = document.unwritten(&mut buf, show.joins_held);
            crate::report::query_warnings::write_document_end(
                written,
                &buf,
                &mut *document,
                xml.run.app,
                writer,
                warnings.as_deref_mut(),
            )?;
            buf.clear();
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2392
            // Queries whose results join the unwritten held report were still read, exported and
            // searched before the stopping batch, so their warnings and requests follow.
            if let Some(warnings) = warnings.as_deref_mut() {
                for later in (q_idx + 1)..xml.queries.len() {
                    warnings.before_query(later, &mut *writer)?;
                }
            }
            break;
        }
        // The query whose id, title and length the result shows (`XmlShow::label`).
        let label = show.label.unwrap_or(q_idx);
        if holds_back {
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2403
            // A report that LOSAT refuses (the title it shows, or that of a subject of its hits):
            // the result is not made, and what the document writes at its arrival (the report
            // held before it) is written; the refusal comes where the report would be written.
            let refusal = match (&texts[label], &texts[q_idx]) {
                (Err(_), _) => text_of(label).err(),
                (_, Err(_)) => text_of(q_idx).err(),
                _ => {
                    name_subjects(&mut subjects, hit_lists.get(q_idx..=q_idx).unwrap_or(&[])).err()
                }
            };
            if let Some(refusal) = refusal {
                let written = document.hold_refused(&mut buf, show.joins_held, refusal);
                crate::report::query_warnings::write_document_end(
                    written,
                    &buf,
                    &mut *document,
                    xml.run.app,
                    writer,
                    warnings.as_deref_mut(),
                )?;
                buf.clear();
                continue;
            }
        }
        let (Ok(text), Ok(shown)) = (&texts[q_idx], &texts[label]) else {
            anyhow::bail!("LOSAT internal error: a refused query text of BLAST XML was written");
        };
        let query = XmlQuery {
            id: shown.query_id(xml.queries[label].local_id.as_bytes()),
            title: &shown.title,
            len: xml.searched_lengths[label] as i32,
            stats: xml.stats[q_idx],
            message: text.message.as_deref(),
            searched: xml.searched.get(q_idx).copied().unwrap_or(true),
            masks: xml.masks.get(q_idx).map_or(&[], Vec::as_slice),
            joins_held: show.joins_held,
        };
        let hit_list = hit_lists.get(q_idx).and_then(Option::as_ref);
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
            let lists = &hit_list.hsplist_array[..];
            // A copy of one record that follows another is shown in its `Hit`.
            let ranges = hit_ranges(
                lists,
                |hsp_list| hsp_list.hsps.first().map(|hsp| hsp.s_idx as usize),
                |previous, s_idx| xml.subjects.joins(previous, s_idx),
            );
            for range in ranges {
                let s_idx = lists[range.start].hsps[0].s_idx as usize;
                let subject = subjects.get(s_idx, PROGRAM)?.names()?;
                let hit = XmlHit::of(subject, xml.subjects.length(s_idx));
                // NCBI's messages about the taxonomy database: opening and the volume's lookup
                // at the first `Hit` that needs them, the lookup of the `Hit` every time.
                posts.hit(sink.position(), &hit, lists[range.start].oid, |texts| {
                    if posting {
                        taxdb_printer.at_names(subject, texts)
                    } else {
                        Ok(())
                    }
                })?;
                sink.hit(&hit, |hsps| {
                    let q_offset = xml.query_records.record(q_idx).map_or(0, |r| r.offset);
                    for hsp in lists[range].iter().flat_map(|hsp_list| &hsp_list.hsps) {
                        let identities = read_num_ident(hsp.num_ident as i64, None);
                        if hsps.unaligned([XmlHsp::unaligned(
                            hsp.raw_score,
                            hsp.bit_score,
                            seqalign_evalue(hsp.e_value),
                            identities,
                        )])? {
                            continue;
                        }
                        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:249-257
                        // NCBI reference (598d8ae6): c++/src/objects/seqalign/Dense_seg.cpp:197-199
                        // An alignment whose query row starts below 0 in record coordinates becomes
                        // an all-gap row 0 in the Dense-seg, so `GetSeqStart` throws `Row is empty`.
                        if hsp.starts_left_of_record(q_offset, xml.ungapped) {
                            return Err(crate::report::blast_xml::CaughtException(
                                super::left_of_query::dense_seg_empty_row_start(),
                            )
                            .into());
                        }
                        xml_rows_into(
                            hsp,
                            xml.query_records,
                            xml.letters,
                            &mut query_row,
                            &mut subject_row,
                        )
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
                            qrow: &query_row,
                            srow: &subject_row,
                            // NCBI reference (598d8ae6):
                            // c++/src/algo/blast/format/blastxml_format.cpp:225-235
                            // (blastxml2_format.cpp:212-222 the same.) An HSP without an edit
                            // script (`-ungapped`) is a Dense-diag: no gaps, one segment over
                            // the row (`count_hsp` of no script), whatever `-` its query row
                            // shows past the query's record (`pairwise::context_letters`).
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
        match crate::report::blast_xml::writer::caught_post(written, &*document, xml.run.app) {
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2203-2230
        // The held report of a `-subject` run is written here, or fails (`CaughtException`: its
        // post comes before the document's end).
        let written = document.finish(&mut buf);
        crate::report::query_warnings::write_document_end(
            written,
            &buf,
            &mut *document,
            xml.run.app,
            writer,
            warnings,
        )?;
    }
    if xml.epilog || xml.stops {
        // The run ends here, normally or with NCBI's error, and NCBI's exit flushes its stream.
        writer.emit()?;
    } else if !writer.start.is_empty() {
        // The call wrote nothing (it holds its report): a later call writes the start.
        document.start_not_written();
    }
    document.keep_out(buf);
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:360-367
/// The output of a report call whose document's start (`start`) is not written yet: NCBI writes the
/// start of the XML2 and JSON documents into its stream's buffer before the search
/// (`PrintProlog`), and the start reaches the output with the stream's first flush (a report's, or
/// the one before a post on the diagnostics, `query_warnings::write_after`), so a crash before
/// that (LOSAT's refusal of the first report) leaves nothing. The start is written before the
/// call's first write, flush or region of the output (`emit`), and the call's other writes pass
/// through unchanged.
struct StartFirst<'w> {
    inner: &'w mut dyn crate::report::ostream::ReportWrite,
    start: Vec<u8>,
}

impl StartFirst<'_> {
    /// Writes the start, if it is not written yet.
    fn emit(&mut self) -> std::io::Result<()> {
        if !self.start.is_empty() {
            let start = std::mem::take(&mut self.start);
            self.inner.write_all(&start)?;
        }
        Ok(())
    }
}

impl std::io::Write for StartFirst<'_> {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.is_empty() {
            return Ok(0);
        }
        self.emit()?;
        self.inner.write(bytes)
    }

    fn write_all(&mut self, bytes: &[u8]) -> std::io::Result<()> {
        if bytes.is_empty() {
            return Ok(());
        }
        self.emit()?;
        self.inner.write_all(bytes)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.emit()?;
        self.inner.flush()
    }
}

impl crate::report::ostream::ReportWrite for StartFirst<'_> {
    fn site(&mut self, site: crate::report::ostream::Site) -> std::io::Result<()> {
        self.emit()?;
        self.inner.site(site)
    }

    fn ncbi_flush(&mut self) -> std::io::Result<()> {
        self.emit()?;
        self.inner.ncbi_flush()
    }

    fn resolve_object(&mut self, pieces: &[usize]) -> std::io::Result<()> {
        self.emit()?;
        self.inner.resolve_object(pieces)
    }

    fn badbit(&mut self) -> std::io::Result<()> {
        self.emit()?;
        self.inner.badbit()
    }

    fn write_fault(&mut self, position: usize) -> std::io::Result<()> {
        self.emit()?;
        self.inner.write_fault(position)
    }

    /// The room of the stream with the start in its buffer, where it waits.
    fn ncbi_room(&self) -> usize {
        self.inner.ncbi_room().saturating_sub(self.start.len())
    }

    fn try_with_inserter(
        &mut self,
        f: &mut dyn FnMut(&mut crate::report::ostream::StreamInserter<'_>) -> std::io::Result<()>,
    ) -> std::io::Result<bool> {
        self.emit()?;
        self.inner.try_with_inserter(f)
    }

    fn insert_cut_run(
        &mut self,
        report: &[u8],
        from: usize,
        cuts: &[(usize, crate::report::ostream::Cut)],
        to: usize,
    ) -> std::io::Result<()> {
        self.emit()?;
        self.inner.insert_cut_run(report, from, cuts, to)
    }
}

/// Makes the names of the subjects with hits in `hit_lists` (`SubjectNameCache::get`): a title that
/// NCBI crashes on is the error.
fn name_subjects(
    subjects: &mut SubjectNameCache<'_, '_>,
    hit_lists: &[Option<BlastnHitList>],
) -> Result<()> {
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            if let Some(hsp) = hsp_list.hsps.first() {
                subjects.get(hsp.s_idx as usize, PROGRAM)?;
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "blast_xml_tests.rs"]
mod tests;
