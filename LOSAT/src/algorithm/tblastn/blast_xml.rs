//! TBLASTN's BLAST XML (outfmt 5): the `CBlastOutput` document of its final hit lists
//! (`report::blast_xml`), one query at a time, once per run.
//!
//! NCBI makes the Iteration of each query from its Seq-aligns (`BlastHitList2SeqAlign_OMF`,
//! see `seqalign`); the XML walks the lists as the Seq-annot formats do: subject by subject in
//! hit-list order, each subject's HSPs sorted by e-value (`render` sorts them), the HSPs of
//! score 0 dropped. What NCBI's `s_SeqAlignSetToXMLHsps` reads of an HSP is its Seq-align
//! (the coordinates in the records, the score, the bit score, the e-value or the sum e-value
//! of a linked set) and the two rows of the alignment, from which it counts identities and
//! positives again: nothing of the engine's counts and nothing of outfmt 0's midline. The rows
//! are the ones outfmt 0 shows (`aligned_sequences`, the subject translated again with B, Z
//! and J for the ambiguous codons) without the lower case of the query's masks; the masks are
//! applied to the query row as NCBI's `s_MaskQuerySeq` does (`report::blast_xml::hsp`).
//!
//! The XML2 and JSON documents of outfmt 15 and 16 (`report::blast_xml2`) are written from the
//! same per-query records by `Xml2Writer`: `write` is generic over the writer
//! (`StructuredXml`). They add the masks of every query (also one without hits), the
//! `searched` flag of the query's batch, the taxid of a database sequence, and `cbs`,
//! `db-gencode` and the gap costs of `XmlValues`.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
//! (`PrintOneResultSet` is called once per query's result set, the accumulated lists hold that
//! one query.)

use anyhow::{ensure, Context, Result};

use super::database::XmlTaxa;
use super::search_gapped::TargetTranslation;
use super::search_seed::tblastn_report_mask_intervals;
use super::stage_d_pipeline::LocalStageDScoring;
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::stage_e_report::{
    aligned_sequences_into, record_coordinates, ReportRanges, ShownSubjects,
};
use super::x_subject_store::SubjectLetters;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::config::ProteinScoringSpec;
use crate::report::blast_xml::hsp::translation_frame;
use crate::report::blast_xml::matrix::FormattingMatrix;
use crate::report::blast_xml::names::{self, QueryText};
use crate::report::blast_xml::subjects::{HitNames, SubjectNames};
use crate::report::blast_xml::{
    filter, hit_ranges, written_queries, QueryMask, StructuredXml, XmlHit, XmlHsp, XmlQuery,
    XmlRun, XmlStats, REFERENCE_GAPPED_BLAST,
};
use crate::report::blast_xml2::{Xml2Run, Xml2Target};
use crate::report::query_warnings::{query_messages, QueryWarnings};
use crate::report::serial::seqalign::{read_num_ident, seqalign_evalue_read};
use crate::stats::tables::lookup_protein_params;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;

const PROGRAM: &str = "TBLASTN";

/// What the document says of the run besides the hit lists (`Parameters`), made once per run.
#[derive(Clone, Copy, Debug)]
pub(super) struct XmlValues {
    /// `Parameters_expect`: the e-value threshold.
    pub expect: f64,
    /// `-soft_masking`: the filter string ends with `m;`.
    pub soft_masking: bool,
    /// Whether the search is gapped (`Parameters_gap-open` and `-gap-extend` of XML2).
    pub gapped: bool,
    /// `-comp_based_stats` as NCBI's number (0 or 2; `Parameters_cbs` of XML2 when non-zero).
    pub cbs: i32,
    /// `-db_gencode` (`Parameters_db-gencode` of XML2).
    pub db_gencode: i32,
    /// IDENTITY searched in mode 0 for another mode: the messages of the searched queries before
    /// this index, the first query batch's, have the warning (`IDENTITY_CBS_RESET_MESSAGE`,
    /// `tblastn_app::search`).
    pub identity_warned: usize,
    /// A list kept no sequence of the database: the messages of the queries of the searched
    /// batches have `EMPTY_DATABASE_MESSAGE`.
    pub filtered_empty: bool,
    /// `-in_pssm` (PSI-TBLASTN): `BlastOutput_program` and `_version` name `psitblastn`
    /// (`Blast_ProgramNameFromType` of the PSSM search's program type).
    pub pssm: bool,
}

/// The document's texts, made (and refused) before the first byte of any format is written:
/// the definition line and message of every query, the names of the subjects with hits, the
/// statistics of every query.
pub(super) struct PreparedXml {
    values: XmlValues,
    /// `None` for a matrix name the formatting table does not know (with a query that is an
    /// error; `prepare`): only a run without a query reaches `write` with it.
    positives: Option<FormattingMatrix>,
    filter: Vec<u8>,
    matrix_name: String,
    gap_open: i32,
    gap_extend: i32,
    db_name: String,
    db_num: i64,
    db_len: i64,
    db_len_overflow: bool,
    /// By query, the texts of the queries whose reports the documents write: every query, or
    /// all but the last of a run that stops when the documents hold that report back
    /// (`written_queries`).
    texts: Vec<QueryText>,
    stats: Vec<XmlStats>,
    /// By subject index: the names of the subjects with hits.
    subject_names: Vec<Option<HitNames>>,
    /// By subject index: `Hit_len`, the length of the subject's record.
    subject_lengths: Vec<i32>,
    /// By rank: the reported database sequences with parsed Seq-ids (empty for `-subject` and
    /// without them), whose copies of one record are shown in one `Hit`
    /// (`db_subject::joins_as_copy`).
    parsed: Vec<Option<std::sync::Arc<crate::report::db_subject::ParsedSubject>>>,
    /// By query: whether NCBI searched the query's batch (it has a valid context in it).
    /// Only XML2 shows it.
    searched: Vec<bool>,
    /// XML2 and JSON: the masks of every query are made, not those of the queries with hits
    /// alone (a query without hits lists its masks too).
    all_masks: bool,
    /// What XML2 needs of the run beyond outfmt 5: the target, the gap costs, `cbs` and
    /// `db-gencode`.
    xml2: Xml2Run,
    /// NCBI's messages about the taxonomy database (XML2 and JSON with `-db`).
    taxa_notes: crate::algorithm::common::database::TaxdbNotes,
    /// The PSSM query's Bioseq letters where they are not its record's
    /// (`PsiSearch::bioseq_letters`), which the rows show and the counts compare as they are
    /// (`row_letters::QueryLetters`); a run of such a query has no other.
    query_bioseq: Option<std::sync::Arc<[u8]>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89,101-103
/// The `Statistics` of every query, with the gapped Karlin block of the scoring (a PSSM query's
/// `kbp_gap_psi`, `PsiSearch::gapped_karlin`, as its pairwise report; an ungapped search: the
/// query's own ungapped block, `ungapped_karlin`); the unsearched values for a query without a
/// valid context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:241-248
/// The XML's Karlin-Altschul values are those of the query's gapped block, or of its ungapped
/// block when it has no gapped one. An ungapped search makes no gapped block, so its statistics
/// are the ungapped block's (XML2's are read the same way, data4xml2format.cpp:258-273).
pub(super) fn query_stats(
    parameters: &LocalSubjectParameters,
    query_validity: &[bool],
    scoring: LocalStageDScoring,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
) -> Vec<XmlStats> {
    let gapped = scoring.gapped.then(|| {
        let matrix = lookup_protein_params(&ProteinScoringSpec {
            matrix: scoring.matrix,
            gap_open: scoring.gap_open,
            gap_extend: scoring.gap_extend,
        });
        match pssm {
            Some(psi) => psi.gapped_karlin(&matrix),
            None => matrix,
        }
    });
    query_validity
        .iter()
        .enumerate()
        .map(|(q_idx, &valid)| match parameters.lengths.get(q_idx) {
            Some(length) if valid => {
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:236-251
                // A gapped block that the ancillary data does not keep (a PSSM file's,
                // `stats::tables::ancillary_kept`): the XML takes the query's ungapped block.
                let karlin = gapped
                    .filter(|karlin| crate::stats::tables::ancillary_kept(karlin.lambda))
                    .unwrap_or(ungapped_karlin[q_idx]);
                XmlStats {
                    hsp_len: length.length_adjustment as i32,
                    eff_space: length.eff_searchsp,
                    kappa: karlin.k,
                    lambda: karlin.lambda,
                    entropy: karlin.h,
                }
            }
            _ => XmlStats::NONE,
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
pub(super) fn prepare(
    values: XmlValues,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    hitlists: &[KappaResultHitList],
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    seg_filter: Option<(i32, f64, f64)>,
    ranges: &ReportRanges<'_>,
    xml2: bool,
    only_xml2: bool,
    taxa: XmlTaxa,
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
) -> Result<PreparedXml> {
    ensure!(
        query_validity.len() == query_records.len()
            && parameters.lengths.len() >= query_records.len(),
        "TBLASTN query statistics count mismatch"
    );
    // A matrix name that the formatting table does not know stops a run with queries before its
    // reports (`stage_e_report::xml_matrix_failure`, `render`): here only a run without a query.
    let positives = FormattingMatrix::new(matrix_name);
    let stats = query_stats(parameters, query_validity, scoring, ungapped_karlin, pssm);
    let has_hits = |hitlist: &KappaResultHitList| {
        hitlist
            .lists()
            .iter()
            .any(|list| !list.hsps.hsps.is_empty())
    };
    let target = match shown {
        ShownSubjects::Fasta { .. } => Xml2Target::Subjects(subject_records.len()),
        ShownSubjects::Database { .. } => Xml2Target::Database,
    };
    // The reports that the documents write (the one report call, so a run that stops is one
    // without the epilog): a run that stops never writes the `-subject` report of its last query
    // when every document of the call holds it back (`written_queries`; outfmt 5 writes it).
    let reports = written_queries(
        query_records.len(),
        !ranges.epilog,
        only_xml2 && target.holds_back(),
    );
    let texts = query_records[..reports]
        .iter()
        .zip(hitlists)
        .zip(query_validity)
        .enumerate()
        .map(|(q_idx, ((query, hitlist), &valid))| {
            names::query_text(
                PROGRAM,
                query.local_id.as_bytes(),
                query.warning_id(),
                query.believed(),
                query.title(),
                true,
                &query_messages(
                    query.seq(),
                    valid,
                    q_idx < values.identity_warned,
                    values.filtered_empty && !query_batch_skipped[q_idx],
                ),
                has_hits(hitlist),
            )
        })
        .collect::<Result<Vec<_>>>()?;
    // The names of the subjects with hits, made once each.
    let mut subject_names: Vec<Option<HitNames>> = Vec::new();
    let count = match shown {
        ShownSubjects::Fasta { .. } => subject_records.len(),
        ShownSubjects::Database { reported, .. } => reported.records.len(),
    };
    subject_names.resize_with(count, || None);
    let mut subject_lengths = vec![0i32; count];
    // What the reports do at each subject's `Hit` (`database::ReportedSubjects::loads`, the one
    // answer; a subject set's records are always loaded).
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let loads = match shown {
        ShownSubjects::Fasta { .. } => &every_loaded,
        ShownSubjects::Database { reported, .. } => &reported.loads,
    };
    for list in hitlists[..reports]
        .iter()
        .flat_map(KappaResultHitList::lists)
    {
        if list.hsps.hsps.is_empty() {
            continue;
        }
        let s_idx = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
        if subject_names[s_idx].is_some() {
            continue;
        }
        let record = subject_records
            .get(s_idx)
            .context("TBLASTN subject OID out of range")?;
        // The `Hit` of a subject whose lookup throws has the names of NCBI's `catch` with the
        // exception, and no names are made for a subject that no report loads (`HitNames::of`).
        subject_names[s_idx] = Some(HitNames::of(
            loads,
            s_idx,
            shown.parsed_subject(s_idx),
            || match shown.oid(s_idx) {
                None => match record.believed() {
                    Some(ids) => SubjectNames::of_believed(ids, record.title(), false),
                    None => {
                        SubjectNames::of_record(record.local_id.as_bytes(), record.title(), false)
                    }
                }
                .map_err(|refusal| refusal.error(PROGRAM, "a subject title")),
                // The taxid is shown by XML2 alone (`taxa` is empty for the other formats).
                Some(db_oid) => Ok(match shown.parsed_subject(s_idx) {
                    Some(parsed) => SubjectNames::of_parsed(
                        parsed,
                        taxa.parsed_taxonomy.get(s_idx).map_or(&[], Vec::as_slice),
                    ),
                    None => SubjectNames::of_database_sequence(
                        db_oid,
                        record.title(),
                        taxa.taxids.get(s_idx).copied().unwrap_or(0),
                        taxa.taxonomy.get(s_idx).cloned().unwrap_or_default(),
                    ),
                }),
            },
        )?);
        // The whole record (`-subject_loc` searches part of it).
        subject_lengths[s_idx] =
            i32::try_from(ranges.subjects.length(s_idx, subject_letters.length(s_idx)))?;
    }
    let (db_name, db_num, db_len, db_len_overflow) = match shown {
        ShownSubjects::Fasta { .. } => (String::new(), 0, 0, false),
        ShownSubjects::Database { name, pairwise, .. } => (
            name.to_string(),
            pairwise.num_sequences,
            pairwise.total_letters,
            pairwise.total_letters_overflow,
        ),
    };
    let xml2_run = Xml2Run {
        target,
        subject_names: match shown {
            ShownSubjects::Fasta { .. } => crate::report::blast_xml2::subject_names(
                subject_records.iter().map(FastaRecord::believed),
            ),
            ShownSubjects::Database { .. } => None,
        },
        gapped: values.gapped,
        cbs: values.cbs,
        query_gencode: None,
        db_gencode: Some(values.db_gencode),
    };
    Ok(PreparedXml {
        values,
        positives,
        // NCBI prints the filtering options (`BlastFilteringOptionsToString` of `SSegOptions`), as
        // blastp, blastx and tblastx do, not the parameters `SegParametersNewAa` makes of them.
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1146-1154
        filter: filter::filter_string(None, seg_filter, None, None, values.soft_masking),
        matrix_name: matrix_name.to_string(),
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        db_name,
        db_num,
        db_len,
        db_len_overflow,
        texts,
        stats,
        subject_names,
        subject_lengths,
        parsed: match shown {
            ShownSubjects::Fasta { .. } => Vec::new(),
            ShownSubjects::Database { reported, .. } => reported.parsed.clone(),
        },
        searched: query_batch_skipped
            .iter()
            .map(|&skipped| !skipped)
            .collect(),
        all_masks: xml2,
        xml2: xml2_run,
        taxa_notes: taxa.notes,
        query_bioseq: pssm.and_then(|psi| psi.bioseq_letters.clone()),
    })
}

/// The masked regions of a query in NCBI's form (`Blast_GetSeqLocInfoVector`): the SEG (and
/// lower-case) intervals of the searched letters (half-open, not yet joined), joined where
/// they overlap (`BlastSeqLocCombine` of `filter_out` with 0, blast_filter.c:1256; adjacent intervals
/// stay apart), other than an empty one and the whole searched interval, as inclusive
/// positions in the record (`offset` is the interval's start), frame 0.
fn query_masks(
    intervals: &[(usize, usize)],
    length: usize,
    offset: usize,
    out: &mut Vec<QueryMask>,
) {
    out.clear();
    let ranges = intervals
        .iter()
        .map(|&(start, end)| (start as i32, end as i32 - 1))
        .collect();
    for (from, to) in crate::algorithm::common::masks::combine_masks(ranges) {
        if from > to || (from == 0 && to == length as i32 - 1) {
            continue;
        }
        out.push(QueryMask {
            from: from + offset as i32,
            to: to + offset as i32,
            frame: 0,
        });
    }
}

/// The masks of a query's outfmt 0 rows and query-anchored sections, half-open in the searched
/// letters: NCBI's list of the query (`query_masks` of its SEG and lower-case intervals,
/// `tblastn_report_mask_intervals`, without an offset), in which a SEG run and a lower-case run
/// that touch stay two masks (two `<font>` spans under `-html`, `x_OutputSeq`) and a mask of the
/// whole searched interval is not listed (the formatter gets `TMaskedQueryRegions`).
pub(super) fn row_masks(
    query: &[u8],
    seg: Option<&SegParams>,
    mask_lowercase: bool,
) -> Vec<(usize, usize)> {
    let mut masks = Vec::new();
    query_masks(
        &tblastn_report_mask_intervals(query, seg, mask_lowercase),
        query.len(),
        0,
        &mut masks,
    );
    masks
        .iter()
        .map(|mask| (mask.from as usize, mask.to as usize + 1))
        .collect()
}

impl PreparedXml {
    /// What the XML2 and JSON documents need of the run beyond outfmt 5's (`Xml2Writer::new`).
    pub(super) fn xml2_run(&self) -> Xml2Run {
        self.xml2.clone()
    }

    /// Writes the document of every query of `hitlists` (indexed by query) through `document`
    /// (BLAST XML, or the XML2 / JSON documents) to `writer`, each query's warnings before
    /// it, with the document's header before the first query and its end after the last when
    /// `ranges.epilog`.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn write<W: StructuredXml>(
        &self,
        document: &mut W,
        writer: &mut dyn crate::report::ostream::ReportWrite,
        query_records: &[FastaRecord],
        subject_letters: SubjectLetters<'_>,
        hitlists: &[KappaResultHitList],
        genetic_code: u8,
        seg: Option<&SegParams>,
        mask_lowercase: bool,
        ranges: &ReportRanges<'_>,
        mut warnings: Option<&mut QueryWarnings<'_>>,
    ) -> Result<()> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blastxml_format.cpp:682-684
        // The program's name, and its upper case with the version (outfmt 7's program line).
        let program_type = super::in_pssm::program_type(self.values.pssm);
        let program = program_type.map_or("tblastn", crate::report::record::ProgramType::name);
        let version = super::tabular::program_version(self.values.pssm);
        let run = self.positives.as_ref().map(|positives| XmlRun {
            program,
            version,
            reference: REFERENCE_GAPPED_BLAST,
            db_name: &self.db_name,
            matrix: Some(&self.matrix_name),
            expect: self.values.expect,
            gap_open: self.gap_open,
            gap_extend: self.gap_extend,
            reward_penalty: None,
            filter: &self.filter,
            db_num: self.db_num,
            db_len: self.db_len,
            db_len_overflow: self.db_len_overflow,
            positives: Some(positives),
            mask_char: b'X',
        });
        // The document's start (`StructuredXml::start`), before the first query's warnings, after
        // the refusal of the database counts that the written reports would print, also for a run
        // that stops before its first report (its first batch fails) and for a run without a
        // query.
        let mut buf: Vec<u8> = Vec::new();
        document.start(run.as_ref().filter(|_| !self.texts.is_empty()), &mut buf)?;
        writer.write_all(&buf)?;
        buf.clear();
        let Some(run) = run else {
            // A run without a query (with one, `render` stops before) and a matrix the formatting
            // table does not know (`-matrix blosum62`: TBLASTN runs it, the table does not):
            // NCBI never builds the table and writes the epilog alone (for outfmt 5 its end of
            // line).
            if ranges.epilog {
                document.finish(&mut buf)?;
                crate::report::query_warnings::write_document_posts(
                    &buf,
                    document.take_posts(),
                    writer,
                    warnings.as_deref_mut(),
                )?;
            }
            return Ok(());
        };
        let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
        let mut taxdb_printer = self.taxa_notes.printer();
        let mut masks: Vec<QueryMask> = Vec::new();
        let mut x_window: Vec<u8> = Vec::new();
        let mut qrow = String::new();
        let mut srow = String::new();
        for (q_idx, (query, hitlist)) in query_records.iter().zip(hitlists).enumerate() {
            // The query's warnings come before its document (`QueryWarnings`).
            if let Some(warnings) = warnings.as_deref_mut() {
                warnings.before_query(q_idx, &mut *writer)?;
                // NCBI's prefetch of the query's batch, where it begins (`TaxdbPrinter::batch_texts`).
                if q_idx < self.texts.len() {
                    warnings.prefetch_lines(&mut *writer, taxdb_printer.batch_texts(q_idx))?;
                }
            }
            if q_idx == self.texts.len() {
                // The held report that is never written (`prepare`, `written_queries`).
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
            let text = &self.texts[q_idx];
            let query_letters = query.seq();
            masks.clear();
            // A query of a batch NCBI did not search has no masks (its result is made
            // without the search, probe `BATCH_SIZE=10` with `-lcase_masking`); XML2 and JSON
            // list the masks of every other query, also one without hits.
            // A PSSM query's results have no masks (`report::record::carries_query_masks`).
            if self.searched[q_idx]
                && crate::report::record::carries_query_masks(program_type)
                && (self.all_masks
                    || hitlist
                        .lists()
                        .iter()
                        .any(|list| !list.hsps.hsps.is_empty()))
            {
                // The masks of the report: the SEG and lower-case intervals, joined as NCBI
                // joins them (`query_masks`).
                let intervals = tblastn_report_mask_intervals(query_letters, seg, mask_lowercase);
                query_masks(
                    &intervals,
                    query_letters.len(),
                    ranges.queries.offset(q_idx),
                    &mut masks,
                );
            }
            let xml_query = XmlQuery {
                id: text.query_id(query.local_id.as_bytes()),
                title: &text.title,
                len: i32::try_from(query_letters.len())?,
                stats: self.stats[q_idx],
                message: text.message.as_deref(),
                // Read by XML2 alone.
                searched: self.searched[q_idx],
                masks: &masks,
            };
            // What NCBI posts where it loads the subject of each `Hit` (`HitPosts`): the texts of
            // the subjects that the query's report reads one by one
            // (`QueryWarnings::subject_read`) and the taxonomy database's messages.
            let mut posts = crate::report::query_warnings::HitPosts::new(
                warnings
                    .as_deref()
                    .is_some_and(|warnings| warnings.reads_with_texts(q_idx)),
            );
            let posting = warnings.is_some();
            let written = document.query(&run, &xml_query, &mut buf, |sink| {
                let lists = hitlist.lists();
                ensure!(
                    lists
                        .iter()
                        .all(|list| list.hsps.hsps.is_empty() || list.oid >= 0),
                    "negative TBLASTN subject OID"
                );
                // A copy of one record that follows another is shown in its `Hit`.
                let hits = hit_ranges(
                    lists,
                    |list| (!list.hsps.hsps.is_empty()).then_some(list.oid as usize),
                    |previous, rank| {
                        crate::report::db_subject::rank_joins_as_copy(&self.parsed, previous, rank)
                    },
                );
                for range in hits {
                    let first = usize::try_from(lists[range.start].oid)
                        .context("negative TBLASTN subject OID")?;
                    let names = self.subject_names[first]
                        .as_ref()
                        .context("names of a subject with hits")?
                        .names()?;
                    let hit = XmlHit::of(names, self.subject_lengths[first]);
                    // NCBI's messages about the taxonomy database: opening and the volume's
                    // lookup at the first `Hit` that needs them, the lookup of the `Hit` every
                    // time (nothing for the other formats).
                    posts.hit(sink.position(), &hit, first as u32, |texts| {
                        if posting {
                            taxdb_printer.at_names((q_idx, first), names, texts)
                        } else {
                            Ok(())
                        }
                    })?;
                    sink.hit(&hit, |hsps| {
                        for list in &lists[range] {
                            if list.hsps.hsps.is_empty() {
                                continue;
                            }
                            let oid = usize::try_from(list.oid)
                                .context("negative TBLASTN subject OID")?;
                            ensure!(
                                list.payloads.len() == list.hsps.hsps.len(),
                                "TBLASTN report payload count mismatch"
                            );
                            let unaligned = list.hsps.hsps.iter().zip(&list.payloads).map(
                                |(linked, payload)| {
                                    XmlHsp::unaligned(
                                        linked.hsp.score,
                                        payload.bit_score,
                                        seqalign_evalue_read(linked.evalue),
                                        read_num_ident(i64::from(payload.num_ident), program_type),
                                    )
                                },
                            );
                            if hsps.unaligned(unaligned)? {
                                continue;
                            }
                            let record_length = self.subject_lengths[oid];
                            let subject_seq = subject_letters.nt(oid);
                            let mut translation = TargetTranslation::of(subject_seq, &code);
                            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                                let hsp = &linked.hsp;
                                let frame = i32::from(hsp.frame);
                                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1087-1105:
                                // translated coordinates use three nucleotide letters per
                                // residue and reverse strand endpoints run in descending order.
                                let subject_len = i32::try_from(subject_seq.len())?;
                                let (s_start, s_end) =
                                    super::stage_e_report::subject_endpoints(hsp, subject_len);
                                // The rows, without the lower case of the masks (a PSSM query's
                                // Bioseq letters as they are: alnvec.cpp:116-121, cited at
                                // `in_pssm::query_residues`).
                                aligned_sequences_into(
                                    crate::report::row_letters::QueryLetters::new(
                                        query_letters,
                                        self.query_bioseq.as_deref(),
                                    ),
                                    subject_seq,
                                    &mut x_window,
                                    &code,
                                    &mut translation,
                                    hsp.frame,
                                    hsp.q_start,
                                    hsp.s_start,
                                    hsp.s_end,
                                    &payload.edit_script,
                                    &mut qrow,
                                    &mut srow,
                                )?;
                                // The coordinates in the searched letters, moved into the records.
                                let (q_start, q_end, s_start, s_end) = record_coordinates(
                                    ranges,
                                    q_idx,
                                    oid,
                                    (hsp.q_start + 1, hsp.q_end, s_start, s_end),
                                )?;
                                // The Seq-align's endpoints, lower first: the minus frame's
                                // outfmt 6 order is not the XML's.
                                let (hit_from, hit_to) = (s_start.min(s_end), s_start.max(s_end));
                                hsps.hsp(&XmlHsp {
                                    score: hsp.score,
                                    bit_score: payload.bit_score,
                                    evalue: seqalign_evalue_read(linked.evalue),
                                    query_from: q_start,
                                    query_to: q_end,
                                    hit_from,
                                    hit_to,
                                    query_frame: 0,
                                    hit_frame: translation_frame(
                                        frame > 0,
                                        hit_from,
                                        hit_to,
                                        record_length,
                                    ),
                                    qrow: qrow.as_bytes(),
                                    srow: srow.as_bytes(),
                                    ops: Some(&payload.edit_script),
                                    query_first: q_start - 1,
                                    query_step: 1,
                                    alignment_identities: read_num_ident(
                                        i64::from(payload.num_ident),
                                        program_type,
                                    ),
                                })?;
                            }
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
            if let Some(post) =
                crate::report::blast_xml::writer::caught_post(written, &*document, run.program)?
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
        if ranges.epilog {
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
}

#[cfg(test)]
#[path = "blast_xml_tests.rs"]
mod tests;
