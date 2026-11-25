//! TBLASTX's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit list
//! (`report::blast_xml`), one query at a time, once per run.
//!
//! TBLASTX is an ungapped search, so NCBI makes one `diags` Seq-align per subject of a query's
//! final hit list, one Std-seg per HSP (`BLASTUngappedHspListToSeqAlign`, see `seqalign`), and
//! the XML writer takes each Std-seg as an `Hsp` of the subject's `Hit`
//! (`PrepareBlastUngappedSeqalign`): the final list (`hit_collector.rs`) holds the HSPs in that
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
//! (`PrintOneResultSet` is called once per query's result set, the accumulated lists hold that
//! one query.)

use anyhow::Result;

use super::blast_engine::{TblastxHsp, TblastxQueryStats, TblastxReportRun};
#[cfg(test)]
use super::ncbi_cutoffs::try_compute_eff_lengths_tblastx;
use super::subjects::SubjectSetNumbers;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names;
use crate::report::blast_xml::{
    filter, written_queries, QueryMask, StructuredXml, SubjectNameCache, XmlHit, XmlHsp, XmlQuery,
    XmlRun, XmlStats, XmlSubjects, REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::{Xml2Run, Xml2Target};
use crate::report::query_warnings::{QueryWarnings, EMPTY_DATABASE_MESSAGE, INVALID_QUERY_MESSAGE};
use crate::report::serial::seqalign::{read_num_ident, seqalign_evalue};
use crate::utils::genetic_code::GeneticCode;

const PROGRAM: &str = "TBLASTX";

/// `-gapopen 11 -gapextend 1`: TBLASTX has no options for them, and NCBI's XML prints its option
/// defaults (blastxml_format.cpp:602-607; DECISIONS 19, Q7), or the best costs of `-matrix`
/// (`gap_costs`).
const GAP_OPEN: i32 = 11;
const GAP_EXTEND: i32 = 1;

/// The gap costs of the XML's parameters: with `-matrix`, the matrix's best costs.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:257-273
fn gap_costs(matrix: Option<&std::ffi::OsStr>) -> (i32, i32) {
    match matrix {
        Some(name) => crate::stats::protein_options::protein_gap_existence_extend_params(
            crate::utils::os_bytes::os_str_bytes(name),
        )
        .unwrap_or((0, 0)),
        None => (GAP_OPEN, GAP_EXTEND),
    }
}

/// The messages of a query's result (the warnings the other formats write to standard error,
/// `report::query_warning`): the Karlin-Altschul message of a batch NCBI did not search (a
/// translated query posts it only when no context of the batch is valid) and the message of a
/// record without letters. TBLASTX posts nothing for an invalid query of a searched batch
/// (blast_stat.c:2783-2790).
fn query_messages(query: &FastaRecord, searched: bool, filtered_empty: bool) -> Vec<String> {
    let mut messages = Vec::new();
    if !searched {
        messages.push(INVALID_QUERY_MESSAGE.to_string());
    }
    if query.seq().is_empty() {
        messages.push(crate::algorithm::blastn::blast_engine::NO_DATA_MESSAGE.to_string());
    }
    if filtered_empty && searched {
        messages.push(EMPTY_DATABASE_MESSAGE.to_string());
    }
    messages
}

/// The `Statistics` of a query: the first valid context's length adjustment, effective search
/// space and ideal ungapped Karlin block (TBLASTX is ungapped); the unsearched values for a query
/// without a valid context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-100
/// The search keeps that context's block, search space and length, not its length adjustment: it
/// is computed again here with the search's own function (`try_compute_eff_lengths_tblastx`,
/// the call of `context_cutoffs` with the same subject set), and the search space it gives
/// must be the kept one.
pub(crate) fn query_stats(stats: &TblastxQueryStats, numbers: &SubjectSetNumbers) -> XmlStats {
    let Some(karlin) = &stats.karlin else {
        return XmlStats::NONE;
    };
    let lengths = numbers.context_lengths(stats.aa_len as i64, karlin).ok();
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
        subject_names: match &run.database {
            Some(_) => None,
            None => crate::report::blast_xml2::subject_names(
                run.subjects.iter().map(FastaRecord::believed),
            ),
        },
        gapped: false,
        cbs: 0,
        query_gencode: Some(i32::from(run.args.query_gencode)),
        db_gencode: Some(i32::from(run.args.db_gencode)),
        program_type: None,
    }
}

/// Writes the document of every query of `run` (`hits` is the final list after
/// `shift_to_records`: record coordinates, in query order, with the Seq-align's `num_ident` of
/// each HSP where the display counts replaced its identities, `report::seqalign_identities`, which
/// `run_impl::report_counts` makes) through
/// `document` (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before it, with the document's
/// header before the first query and its end after the last when `run.epilog`. A query or
/// subject title that NCBI crashes on (`names::Refusal`) is refused before the first byte; every
/// other text is written as NCBI's streams write it. `all_masks`: the masks of every searched
/// query are listed (XML2 and JSON), not those of the queries with hits alone (outfmt 5).
/// `run.xml_taxa` holds the taxids of the database sequences that XML2 and JSON show.
pub(crate) fn write_tblastx_xml<W: StructuredXml>(
    document: &mut W,
    all_masks: bool,
    only_xml2: bool,
    writer: &mut dyn crate::report::ostream::ReportWrite,
    (hits, identities): (&[TblastxHsp], Option<&[u32]>),
    run: &TblastxReportRun<'_>,
    mut warnings: Option<&mut QueryWarnings<'_>>,
) -> Result<()> {
    let args = run.args;
    let taxa = run.xml_taxa;
    let (db_name, db_num, db_len, db_len_overflow) = match &run.database {
        Some(database) => (
            database.name,
            database.pairwise.num_sequences,
            database.pairwise.total_letters,
            database.pairwise.total_letters_overflow,
        ),
        None => (b"".as_slice(), 0, 0, false),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2206-2243 (XML2
    // and JSON) and :2242-2243 (BLAST XML)
    // The batch loop refused a name that the formatting table does not know when a batch has
    // letters (`run_batch_loop`); a run without queries makes no report data, and the document
    // is its header and end.
    let matrix = args.matrix_name();
    let positives = FormattingMatrix::new(matrix.as_bytes());
    anyhow::ensure!(
        positives.is_some() || hits.is_empty(),
        "no formatting matrix for {matrix}"
    );
    let (gap_open, gap_extend) = gap_costs(args.matrix.as_deref());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:414-417
    // `-soft_masking true` (`mask_at_hash`) ends the filter string with `m;`.
    let filter = filter::filter_string(
        None,
        args.seg_spec()?.filter_options(),
        None,
        None,
        args.soft_masking(),
    );
    let xml_run = XmlRun {
        program: "tblastx",
        version: concat!("TBLASTX ", crate::wording::report_version!()),
        reference: REFERENCE_GAPPED_BLAST,
        db_name,
        matrix: Some(&matrix),
        expect: args.evalue,
        gap_open,
        gap_extend,
        reward_penalty: None,
        filter: &filter,
        db_num,
        db_len,
        db_len_overflow,
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
    // hits, made (and refused) before the first byte, for the reports that the document writes:
    // a run that stops (the one report call without the epilog: the XML formats are not streamed,
    // `report_stream`) never writes the `-subject` report of its last query, which NCBI holds
    // back, unless another format of the call (outfmt 5) writes it (`written_queries`; `only_xml2`:
    // every structured format of the call is XML2 or JSON).
    let reports = written_queries(
        queries.len(),
        !run.epilog,
        only_xml2 && document.holds_back(),
    );
    let texts = queries[..reports]
        .iter()
        .enumerate()
        .map(|(q_idx, query)| {
            names::query_text(
                PROGRAM,
                query.local_id.as_bytes(),
                query.warning_id(),
                query.believed(),
                query.title(),
                false,
                &query_messages(query, run.searched[q_idx], run.filtered_empty),
                !query_hits[q_idx].is_empty(),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    let subject_length =
        |s_idx: usize| subject_placements.length(s_idx, run.subject_letters.length(s_idx));
    let subjects = match &run.database {
        Some(database) => XmlSubjects::Database {
            ordinals: &database.reported.ordinals,
            oids: database.oids,
            titles: crate::report::blast_xml::XmlTitles::Records(run.subjects),
            // Only XML2 shows the taxid (`taxa` is empty for outfmt 5).
            taxids: &taxa.taxids,
            taxonomy: &taxa.taxonomy,
            length: &subject_length,
            parsed: &database.reported.parsed,
            parsed_taxonomy: &taxa.parsed_taxonomy,
            loads: &database.reported.loads,
        },
        None => XmlSubjects::Records {
            records: run.subjects,
            length: &subject_length,
            protein: false,
        },
    };
    let mut cache = SubjectNameCache::new(&subjects);
    let mut last_subject = None;
    let written_hits = query_hits
        .get(reports)
        .map_or(hits.len(), |range| range.start);
    for hsp in &hits[..written_hits] {
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
    // The document's start (`StructuredXml::start`), before the first query's warnings, after the
    // refusal of the database counts that the written reports would print.
    document.start((reports > 0).then_some(&xml_run), &mut buf)?;
    writer.write_all(&buf)?;
    buf.clear();
    for (q_idx, query) in queries.iter().enumerate() {
        // The query's warnings come before its document (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
            // NCBI's prefetch of the query's batch, where it begins (`TaxdbPrinter::batch_texts`).
            if q_idx < reports {
                warnings.prefetch_lines(&mut *writer, taxdb_printer.batch_texts(q_idx))?;
            }
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
        let text = &texts[q_idx];
        let searched_length = query.seq().len();
        let record_length = query_placements.length(q_idx, searched_length) as i32;
        let query_start = query_hits[q_idx].start;
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
            id: text.query_id(query.local_id.as_bytes()),
            title: &text.title,
            len: i32::try_from(searched_length)?,
            stats: query_stats(&run.query_stats[q_idx], run.numbers),
            message: text.message.as_deref(),
            // Read by XML2 alone.
            searched: run.searched[q_idx],
            masks: &masks,
        };
        // What NCBI posts where it loads the subject of each `Hit` (`HitPosts`): the texts of the
        // subjects that the query's report reads one by one (`QueryWarnings::subject_read`) and
        // the taxonomy database's messages.
        let mut posts = crate::report::query_warnings::HitPosts::new(
            warnings
                .as_deref()
                .is_some_and(|warnings| warnings.reads_with_texts(q_idx)),
        );
        let posting = warnings.is_some();
        let written = document.query(&xml_run, &xml_query, &mut buf, |sink| {
            let mut first = 0usize;
            while first < query_hits.len() {
                let s_idx = query_hits[first].hit.s_idx;
                // The subject's HSPs and those of the copies of one record that follow it
                // (`XmlSubjects::joins`, as `hit_ranges`).
                let mut end = first + 1;
                while end < query_hits.len() {
                    let (previous, next) =
                        (query_hits[end - 1].hit.s_idx, query_hits[end].hit.s_idx);
                    if previous != next && !subjects.joins(previous as usize, next as usize) {
                        break;
                    }
                    end += 1;
                }
                let subject = cache.get(s_idx as usize, PROGRAM)?.names()?;
                let subject_record_length = subjects.length(s_idx as usize);
                let hit = XmlHit::of(subject, subject_record_length);
                // NCBI's messages about the taxonomy database: opening and the volume's lookup
                // at the first `Hit` that needs them, the lookup of the `Hit` every time
                // (nothing for the other formats).
                posts.hit(sink.position(), &hit, s_idx, |texts| {
                    if posting {
                        taxdb_printer.at_names((q_idx, s_idx as usize), subject, texts)
                    } else {
                        Ok(())
                    }
                })?;
                sink.hit(&hit, |hsps| {
                    for (index, hsp) in (first..end).zip(&query_hits[first..end]) {
                        let identities = read_num_ident(
                            identities.map_or(hsp.hit.num_ident as i64, |identities| {
                                i64::from(identities[query_start + index])
                            }),
                            None,
                        );
                        if hsps.unaligned([XmlHsp::unaligned(
                            hsp.hit.raw_score,
                            hsp.hit.bit_score,
                            seqalign_evalue(hsp.hit.e_value),
                            identities,
                        )])? {
                            continue;
                        }
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
                            alignment_identities: identities,
                        })?;
                    }
                    Ok(())
                })?;
                first = end;
            }
            Ok(())
        });
        // A query whose report raised an exception that NCBI's report function catches
        // (`CaughtException`: the `Hit` of a subject whose lookup throws, `HitNames`): the XML2 and
        // JSON documents leave its report out (`Xml2Document::report`), and the post of the
        // exception follows the texts of the loads before it (`HitPosts::caught`). Another failed
        // query leaves its partial bytes in `buf`, which the error return drops.
        if let Some(post) =
            crate::report::blast_xml::writer::caught_post(written, &*document, xml_run.program)?
        {
            posts.caught(post);
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
    if run.epilog {
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

#[cfg(test)]
#[path = "blast_xml_tests.rs"]
mod tests;
