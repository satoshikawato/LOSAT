//! NCBI TBLASTN local-subject result formatting.

use std::io::Write;
use std::path::Path;
use std::sync::Arc;

use anyhow::{bail, ensure, Context, Result};

use super::search_gapped::TargetTranslation;
use super::search_seed::encode_tblastn_lookup_query;
use super::stage_d_pipeline::{LocalStageDScoring, SubjectSetNumbers};
use super::stage_d_results::KappaResultHitList;
use super::stage_d_stats::LocalSubjectParameters;
use super::x_subject_store::SubjectLetters;
use crate::api::local_blast::{FormatProbe, HspIndex, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::common::{GapEditOp, Hit};
use crate::config::ProteinScoringSpec;
use crate::core::composition_adjustment::redo_alignment::EMatrixAdjustRule;
use crate::report::outfmt6::{
    format_bitscore_ncbi, format_evalue_ncbi_tabular, format_percent_identity_ncbi,
};
use crate::report::pairwise::{
    write_tblastn_pairwise_report, BlastpPairwiseQuery, BlastpPairwiseReport, PairwiseConfig,
    PairwiseDatabase, PairwiseHit,
};
use crate::report::row_letters::{stored_gap_columns, QueryLetters, RowSide};
use crate::sequence::packed_subject::SubjectNt;
use crate::stats::tables::lookup_protein_params;
use crate::utils::genetic_code::GeneticCode;
use crate::utils::seg::SegParams;

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:385-391
/// The subject set as the reports name it.
pub(super) enum ShownSubjects<'a> {
    /// `-subject`: `User specified sequence set (Input: PATH)`; the records are the whole set,
    /// and `parsed` the ones whose ids `-parse_deflines` read (empty without the option).
    Fasta {
        path: &'a Path,
        parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    },
    /// `-db`: the `-db` value (outfmt 7's `# Database:` line), the database's `Database:`
    /// lines (outfmt 0) and the reported sequences by rank (their tabular ids).
    Database {
        name: &'a str,
        pairwise: &'a PairwiseDatabase,
        reported: &'a super::database::ReportedSubjects,
        /// The OID of each ordinal id, which names a sequence without parsed Seq-ids.
        oids: crate::seqdb::OidMap<'a>,
        /// `-db_hard_mask` (series NP, M4b PR 2): the search read `B` over the masked bases
        /// (`XDbStream::read_stage`), and the reports count the stored letters
        /// (`recounts_shown_counts`).
        hard_masked: bool,
    },
}

impl ShownSubjects<'_> {
    /// The name after `# Database: ` (outfmt 7) and in the pairwise report of a subject set.
    pub(super) fn name(&self) -> String {
        match self {
            Self::Fasta { path, .. } => {
                format!("User specified sequence set (Input: {})", path.display())
            }
            Self::Database { name, .. } => name.to_string(),
        }
    }

    /// Whether the search read other letters of the database sequences than the reports show
    /// (`-db_hard_mask`: `ShownSubjects::Database::hard_masked`).
    fn hard_masked(&self) -> bool {
        matches!(
            self,
            Self::Database {
                hard_masked: true,
                ..
            }
        )
    }

    /// The `Database:` lines of a database search's pairwise report.
    fn pairwise_database(&self) -> Option<&PairwiseDatabase> {
        match self {
            Self::Fasta { .. } => None,
            Self::Database { pairwise, .. } => Some(pairwise),
        }
    }

    /// The database title of the Seq-annots (`RunRecord::db_title`; NCBI's
    /// `m_DbInfo.front().definition`, the `Database:` line's title), `None` for `-subject`.
    pub(super) fn db_title(&self) -> Option<&[u8]> {
        self.pairwise_database()
            .map(|pairwise| pairwise.definition.as_slice())
    }

    /// The OID of the reported database sequence of rank `rank`, `None` for `-subject`.
    pub(super) fn oid(&self, rank: usize) -> Option<crate::seqdb::Oid> {
        match self {
            Self::Fasta { .. } => None,
            Self::Database { reported, oids, .. } => Some(oids.oid(reported.ordinals[rank])),
        }
    }

    /// The reported database sequences with parsed Seq-ids by rank, or the records whose ids
    /// `-parse_deflines` read by index (empty without them).
    pub(super) fn parsed_subjects(
        &self,
    ) -> &[Option<Arc<crate::report::db_subject::ParsedSubject>>] {
        match self {
            Self::Fasta { parsed, .. } => parsed,
            Self::Database { reported, .. } => &reported.parsed,
        }
    }

    /// The sources from which NCBI makes the empty stored titles of the reported database
    /// sequences, by rank (empty for the records of a subject set, and where none has one).
    pub(super) fn subject_sources(&self) -> &[crate::report::defline::SourceTaxname] {
        match self {
            Self::Fasta { .. } => &[],
            Self::Database { reported, .. } => &reported.sources,
        }
    }
}

impl<'a> ShownSubjects<'a> {
    /// The reported database sequence of rank `rank` when it has parsed Seq-ids, or the record
    /// of index `rank` whose ids `-parse_deflines` read.
    pub(super) fn parsed_subject(
        &self,
        rank: usize,
    ) -> Option<&'a crate::report::db_subject::ParsedSubject> {
        match *self {
            Self::Fasta { parsed, .. } => parsed.get(rank).and_then(Option::as_deref),
            Self::Database { reported, .. } => reported.parsed_subject(rank),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:759-835,1411-1458:
// The formatter receives already sorted query hitlists from Stage D.
// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467:
// Each requested format prints the same final result without searching again.
/// The number of the BLAST archive among the formats (outfmt 11).
pub(super) const ARCHIVE_FORMAT: i32 = 11;

/// The format number of a checked `-outfmt` value, as NCBI's `ParseFormattingString` reads
/// it (a custom specification of outfmt 0 is ignored); -1 for a value it does not read.
pub(super) fn format_number(outfmt: &str) -> i32 {
    crate::blastinput::app::parse_formatting_string(outfmt).map_or(-1, |choice| choice.number)
}

/// The XML2 / JSON document of a checked `-outfmt` value (16, 15, 14 or 13; `None` for the
/// others). Outfmt 13 and 14 write side files next to `-out`, or to the caller's receiver of them
/// (`OutputSink::Files`); a `run_local` call with neither stops at NCBI's XInclude check.
pub(super) fn xml2_kind(outfmt: &str) -> Option<crate::report::blast_xml2::Xml2Kind> {
    use crate::blastinput::app::ReportFormat;
    crate::report::blast_xml2::Xml2Kind::of(ReportFormat::from_number(format_number(outfmt))?)
}

/// Whether a checked `-outfmt` value is a format of the structured XML family that takes the
/// path of outfmt 5: BLAST XML (5) and the XML2 and JSON documents (16, 15, 14, 13).
pub(super) fn is_structured_xml(outfmt: &str) -> bool {
    format_number(outfmt) == 5 || xml2_kind(outfmt).is_some()
}

/// Whether the reports of `outputs` read the letters of the subjects
/// (`formats_read_subject_letters`), for the PSSM query `pssm` (`-in_pssm`) when there is one
/// and in a search with `-db_hard_mask` (`hard_masked`) (`recounts_shown_counts`). The one
/// decision of whether the reported database sequences keep their letters
/// (`database::reported_subjects`) and whether `render` counts again.
pub(super) fn reads_subject_letters(
    outputs: &ReportOutputs<'_>,
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
    hard_masked: bool,
) -> bool {
    formats_read_subject_letters(
        outputs.hits.is_some(),
        recounts_shown_counts(pssm, hard_masked),
        outputs.formats.iter().map(|format| format.outfmt),
    )
}

/// Whether the reports `outfmts` (with `hits`, the caller's hit records) read the letters of the
/// subjects: outfmt 0 and the hit records (`pairwise_hits`), the XML family
/// (`blast_xml::XmlReport::write`) and the tabular rows (`tabular::reads_rows`), and with
/// `recounted` (a query whose reports count its HSPs again from their rows,
/// `recounts_shown_counts`) every tabular report, whose identities, positives and mismatches are
/// then the rows'. The Seq-aligns and the archive (8, 9, 11, 12) and the other tabular reports
/// (`write_tabular`, the default rows of 6 and 7, has no spec) read the lengths alone: their
/// counts are the search's. A value that this does not know reads them.
fn formats_read_subject_letters<'a>(
    hits: bool,
    recounted: bool,
    mut outfmts: impl Iterator<Item = &'a str>,
) -> bool {
    hits || outfmts.any(|outfmt| match format_number(outfmt) {
        8 | 9 | 11 | 12 => false,
        6 | 7 | 10 | 20 => {
            recounted
                || super::tabular::tabular_spec(outfmt)
                    .is_some_and(|spec| super::tabular::reads_rows(&spec))
        }
        // The organism report (`report::tax_report`) shows scores and definition lines.
        18 => false,
        _ => true,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:121, :196 (cited at
// `blast_xml2::reports_held`)
/// Where TBLASTN's run of the formats `formats` stops for a matrix name `matrix_name` that NCBI's
/// BLAST XML report data does not know (`blast_xml::matrix::xml_formatting_failure`), with the
/// facts of the search: no database (`-subject`, `m_IsBl2Seq`) and `-in_pssm` (PSI-TBLASTN,
/// `m_IsIterative`). The one decision behind the search strategy's last request
/// (`tblastn_app::search`) and the error of `render`.
pub(super) fn xml_matrix_failure(
    formats: &[crate::api::local_blast::FormatOutput<'_>],
    matrix_name: &str,
    database: bool,
    pssm: bool,
) -> Option<crate::report::blast_xml::matrix::XmlMatrixFailure> {
    crate::report::blast_xml::matrix::xml_formatting_failure(
        formats.iter().map(|format| format.outfmt),
        matrix_name,
        !database,
        pssm,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4650-4653,
// c++/src/algo/blast/api/blast_seqalign.cpp:181-182,201-204 and
// c++/src/objects/seq/seq_align_mapper_base.cpp:1142-1151
/// The refusal of a reported alignment whose traceback started on its frame's NULLB (TBLASTN
/// `-in_pssm -window_size 0`, a PSSM whose '-' column scores above BLAST_SCORE_MIN): the only one
/// whose edit script holds more query letters than its range, as NCBI's ends stop one letter
/// before the start (the right extension did not run) while the left one holds the start's
/// column, the NULLB. NCBI's Seq-align, made from the start and the edit script, places that
/// column past the subject's last base: with `-db` the reports take their ends from it and read
/// the subject past its end (C++ exceptions where the frame's last codon is not the subject's
/// last), and with `-subject` the remapping onto the subject's location truncates it and drops
/// every score of the alignment. LOSAT does not reproduce those reports yet (line RG item X2d-2).
#[cold]
#[inline(never)]
fn x_refuse_a_column_of_the_nullb(hitlists: &[KappaResultHitList]) -> Result<()> {
    for list in hitlists.iter().flat_map(KappaResultHitList::lists) {
        for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
            let query_letters: u64 = payload
                .edit_script
                .iter()
                .map(|op| match *op {
                    GapEditOp::Sub(n) | GapEditOp::Ins(n) => u64::from(n),
                    GapEditOp::Del(_) => 0,
                })
                .sum();
            let range = i64::from(linked.hsp.q_end) - i64::from(linked.hsp.q_start);
            if i64::try_from(query_letters)? != range {
                bail!(
                    "TBLASTN -in_pssm with -window_size 0: a reported alignment of a one-hit traceback that starts on the NULLB after its frame (the PSSM's '-' column scores it) holds a column of that NULLB, which NCBI BLAST+ reports past the end of the subject; LOSAT does not support it yet"
                );
            }
        }
    }
    Ok(())
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
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    subject_set: SubjectSetNumbers,
    hitlists: &mut [KappaResultHitList],
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    report_counts: crate::blastinput::report_counts::ReportCounts,
    formatting: crate::blastinput::report_counts::FormattingArgs,
    genetic_code: u8,
    seg: Option<&SegParams>,
    // The `-seg` options as NCBI's filter string prints them (`SegSpec::filter_options`), not the
    // SEG parameters `seg` (a value not above 0 or NaN there is SEG's default, a hicut below the
    // locut the locut): BLAST XML's `Parameters_filter`.
    seg_filter: Option<(i32, f64, f64)>,
    mask_lowercase: bool,
    query_posts: &[crate::report::diag_log::Posts],
    // The texts of the subjects that each query's report reads one by one (its batch's prefetch
    // off), by query, which the first format's writer places where it reads them
    // (`QueryWarnings::subject_read`).
    query_reads: &[crate::report::query_warnings::QueryReads],
    ranges: &ReportRanges<'_>,
    xml_values: super::blast_xml::XmlValues,
    xml_taxa: super::database::XmlTaxa,
    out: Option<&std::path::Path>,
    archive: Option<&super::archive::ArchiveRender<'_>>,
    tabular: &super::tabular::TabularInputs,
    // `-in_pssm` (PSI-TBLASTN): the query's PSSM, whose reports (outfmt 0, 6, 7 and 10) name the
    // program PSITBLASTN, show no query masks and print the PSSM's gapped Karlin-Altschul values.
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
    // `-show_gis`: the labels of the query-anchored sections (outfmt 1 to 4).
    show_gis: bool,
    // The first searched query of each query batch (NCBI's scope of the batch's queries, which
    // the query-anchored formats' refusals read); empty: one batch.
    batch_starts: &[usize],
) -> Result<()> {
    ensure!(
        query_records.len() == hitlists.len(),
        "TBLASTN report query count mismatch"
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:94
    let program = if pssm.is_some() {
        super::in_pssm::PROGRAM
    } else {
        "TBLASTN"
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1572-1577:
    // Blast_HSPListSortByEvalue(hsp_list); /* during Seq-align conversion */
    for hitlist in hitlists.iter_mut() {
        hitlist.sort_hsps_for_report();
        hitlist.drop_zero_score_hsps();
    }
    // A PSSM query's Bioseq letters where they are not its record's, which the reports read; its
    // counts are the rows' where its reports read other letters than its search and a report
    // reads the subjects' letters (`reads_subject_letters`; the Seq-aligns keep the search's).
    let query_bioseq = pssm.and_then(|psi| psi.bioseq_letters.as_deref());
    let hard_masked = shown.hard_masked();
    if recounts_shown_counts(pssm, hard_masked) && reads_subject_letters(outputs, pssm, hard_masked)
    {
        recount_shown_counts(
            query_records,
            query_bioseq,
            subject_letters,
            hitlists,
            genetic_code,
            hard_masked,
        )?;
    }
    let hitlists: &[KappaResultHitList] = hitlists;
    // `-in_pssm` with `-window_size 0`: an alignment whose traceback started on its frame's NULLB
    // (`search_gapped::x_traceback_start_past_frame`) is refused here, where NCBI reports it, and
    // not where NCBI drops it after its traceback (the hit list, `-max_target_seqs`).
    if pssm.is_some() {
        x_refuse_a_column_of_the_nullb(hitlists)?;
    }
    // BLAST XML (outfmt 5, and the XML2 and JSON documents of 15 and 16): its texts are made,
    // and a title that NCBI crashes on (`names::Refusal`) is refused, before any format writes
    // (`-out` is already created, empty, as NCBI's).
    //
    // NCBI writes the XML2 and JSON documents' start before the search (`PrintProlog`) and
    // fails for a matrix of its formatting table's that it does not know when it makes a report's
    // data (`xml_matrix_failure`): that of the run's first query, or of its second where the
    // report is held (`XmlMatrixFailure::SecondQueryOrEnd`: NCBI's next query of another Seq-id,
    // the second as LOSAT refuses two queries of one id before any output), in the batch of that
    // query, or for a run of one query whose report is held
    // in `PrintEpilog` (blast_format.cpp:2204-2212, cited at
    // `blast_xml::matrix::XmlMatrixFailure`) unless the run stops before (`ranges.epilog`; the
    // caller writes that stop). The start alone is written then (probe `tblastn_sub_lcmatrix`),
    // and the posts of the queries up to the failing one: the reader's messages of their batches
    // and NCBI's posts of the bytes replaced in the batches' search strategy requests.
    // (The search strategy's file holds the requests up to the failing batch's: `search` stops
    // the requests after it with the same decision, `StrategyExport::stop_after`.)
    if let Some(failure) = xml_matrix_failure(
        &outputs.formats,
        matrix_name,
        matches!(shown, ShownSubjects::Database { .. }),
        pssm.is_some(),
    )
    .filter(|_| !query_records.is_empty())
    {
        for format in outputs.formats.iter_mut() {
            if let Some(kind) = xml2_kind(format.outfmt) {
                let mut writer = format.sink.open()?;
                writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                writer.flush()?;
            }
        }
        let failing = (failure.query() < query_records.len()).then(|| failure.query());
        for posts in query_posts
            .iter()
            .take(failing.map_or(query_records.len(), |query| query + 1))
        {
            posts.write(
                &mut *outputs.diagnostics,
                archive.map(|render| render.handle),
            )?;
        }
        if failing.is_some() || (ranges.epilog && failure.fails_at_end(query_records.len())) {
            return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
        }
        return Ok(());
    }
    let xml_is_xml2 = outputs
        .formats
        .iter()
        .any(|format| xml2_kind(format.outfmt).is_some());
    // Every structured format is XML2 or JSON (outfmt 5 is not among them; the library takes
    // several formats).
    let xml_only_xml2 = outputs
        .formats
        .iter()
        .filter(|format| is_structured_xml(format.outfmt))
        .all(|format| xml2_kind(format.outfmt).is_some());
    let xml = if outputs
        .formats
        .iter()
        .any(|format| is_structured_xml(format.outfmt))
    {
        Some(super::blast_xml::prepare(
            xml_values,
            query_records,
            subject_records,
            subject_letters,
            shown,
            hitlists,
            parameters,
            ungapped_karlin,
            query_validity,
            query_batch_skipped,
            scoring,
            matrix_name,
            seg_filter,
            ranges,
            xml_is_xml2,
            xml_only_xml2,
            xml_taxa,
            pssm.and_then(|psi| psi.bioseq_letters.clone()),
        )?)
    } else {
        None
    };
    // The tabular subject IDs, made before any format is written: a subject that NCBI
    // writes part of a row for and then stops at (RP-20) is rejected first
    // (`tabular_subject_ids`). The tabular formats 6, 7, 10 and 20 show them, not the organism
    // report (`report::tabular::writes_rows`).
    let shows_ids = |format: &crate::api::local_blast::FormatOutput<'_>| {
        matches!(format_number(format.outfmt), 6 | 7 | 10 | 18 | 20)
            && crate::report::tabular::writes_rows(
                super::tabular::tabular_spec(format.outfmt).as_ref(),
            )
    };
    let subject_ids = if outputs.formats.iter().any(shows_ids) {
        match shown {
            ShownSubjects::Fasta { .. } => tabular_subject_ids(subject_records),
            ShownSubjects::Database { reported, oids, .. } => reported.tabular_subject_ids(
                *oids,
                outputs.formats.iter().any(|format| {
                    shows_ids(format)
                        && crate::report::tabular::reads_subject(
                            super::tabular::tabular_spec(format.outfmt).as_ref(),
                            crate::blastinput::app::Program::Tblastn,
                        )
                }),
            )?,
        }
    } else {
        Vec::new()
    };
    // The rows of subjects whose lookup fails (no id) where NCBI has no result
    // (`report::tabular::check_failed_lookups`), before any format is written.
    let parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>] = match shown {
        ShownSubjects::Fasta { parsed, .. } => parsed,
        ShownSubjects::Database { reported, .. } => &reported.parsed,
    };
    let database = matches!(shown, ShownSubjects::Database { .. });
    // A subject set's subjects are all loaded.
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let loads = match shown {
        ShownSubjects::Fasta { .. } => &every_loaded,
        ShownSubjects::Database { reported, .. } => &reported.loads,
    };
    // Each tabular format's refusal is that format's (`FormatRefusal::check`).
    for format in outputs
        .formats
        .iter_mut()
        .filter(|format| shows_ids(format))
    {
        let spec = super::tabular::tabular_spec(format.outfmt);
        if crate::algorithm::common::tabular::checks_failed_lookups(
            spec.as_ref(),
            &subject_ids,
            parsed,
            loads,
        ) {
            let checked = crate::report::tabular::check_failed_lookups(
                spec.as_ref(),
                hitlists.iter().map(|hitlist| {
                    hitlist
                        .lists()
                        .iter()
                        .filter(|list| !list.hsps.hsps.is_empty())
                        .map(|list| list.oid as usize)
                }),
                |s_idx| {
                    crate::algorithm::common::tabular::subject_lookup_failure(
                        spec.as_ref(),
                        &subject_ids,
                        parsed,
                        loads,
                        database,
                        false,
                        s_idx,
                    )
                },
                |s_idx| match shown {
                    ShownSubjects::Database { reported, oids, .. } => {
                        crate::report::tabular::failed_database_subject(
                            oids.oid(reported.ordinals[s_idx]),
                            loads.thrown(s_idx),
                        )
                    }
                    ShownSubjects::Fasta { .. } => crate::report::tabular::failed_record_title(
                        s_idx,
                        subject_records[s_idx].shown_id(),
                    ),
                },
                "TBLASTN",
            );
            format.refusal.check(checked)?;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411:
    // The alignments are rendered once from the final result, for outfmt 0 and for
    // the caller's hit records, in the order of the tabular rows.
    let anchored_formats = outputs
        .formats
        .iter()
        .any(|format| (1..=4).contains(&format_number(format.outfmt)));
    // The masks of the rows, kept for the sections of outfmt 1 to 4 only (`pairwise_hits`).
    let mut section_masks: Vec<Vec<(usize, usize)>> = Vec::new();
    let pairwise_hits = if outputs.hits.is_some()
        || outputs
            .formats
            .iter()
            .any(|format| (0..=4).contains(&format_number(format.outfmt)))
    {
        Some(std::sync::Arc::new(pairwise_hits(
            query_records,
            subject_records,
            subject_letters,
            hitlists,
            scoring,
            genetic_code,
            seg,
            mask_lowercase,
            ranges,
            pssm.is_none(),
            query_bioseq,
            anchored_formats.then_some(&mut section_masks),
        )?))
    } else {
        None
    };
    if let (Some(hits_sink), Some(pairwise_hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref())
    {
        hits_sink(std::sync::Arc::clone(pairwise_hits));
    }
    // The query-anchored formats (outfmt 1 to 4, line RG T3b), before any format of the call is
    // written: a reported database sequence that NCBI's scope finds as a query of its batch, by
    // its Seq-align id (a generated `lcl|Query_N`) or, in a description line, by its display id
    // (`common::database::refuse_query_ids_of_anchored_reports`, as BLASTN's): NCBI shows the
    // query in its place (with the protein query's letters translated as nucleotides).
    let anchored_view = outputs.formats.iter().find_map(|format| {
        crate::blastinput::app::ReportFormat::from_number(format_number(format.outfmt))
            .and_then(crate::report::pairwise::PairwiseView::of)
            .filter(|view| *view != crate::report::pairwise::PairwiseView::Pairwise)
    });
    if let (Some(view), Some(rows)) = (anchored_view, pairwise_hits.as_deref()) {
        if matches!(shown, ShownSubjects::Database { .. }) {
            let batches = crate::blastinput::query_batch::AnchoredBatchScopes::of_starts(
                query_records,
                batch_starts,
            );
            let parsed = shown.parsed_subjects();
            crate::algorithm::common::database::refuse_query_ids_of_anchored_reports(
                query_records.len(),
                |q_idx| {
                    hitlists[q_idx]
                        .lists()
                        .iter()
                        .filter(|list| !list.hsps.hsps.is_empty())
                        .filter_map(|list| u32::try_from(list.oid).ok())
                        .collect()
                },
                |q_idx| {
                    Ok(crate::report::pairwise::blastp_described_subjects(
                        view.report_sort(formatting.report_sort()),
                        None,
                        report_counts.num_descriptions,
                        q_idx,
                        rows,
                        |previous, s_idx| {
                            crate::report::db_subject::rank_joins_as_copy(
                                parsed,
                                previous as usize,
                                s_idx as usize,
                            )
                        },
                    ))
                },
                parsed,
                |rank| {
                    shown
                        .oid(rank as usize)
                        .expect("a database search names its sequences by OID")
                },
                |q_idx| batches.of_query(q_idx),
                batches.generated(),
                program,
            )?;
        }
    }
    // The sequences of the query-anchored sections of outfmt 1 to 4 (line RG T3b), made once.
    let anchored = if anchored_formats {
        let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
        Some(super::anchored::AnchoredSequences {
            queries: query_records,
            query_bioseq,
            query_placements: ranges.queries,
            subject_records,
            subject_letters,
            subject_placements: ranges.subjects,
            shown,
            show_gis,
            // The masks of outfmt 0's rows (`pairwise_hits`), half-open in the searched letters.
            query_masks: (0..query_records.len())
                .map(|q_idx| {
                    let offset = ranges.queries.offset(q_idx) as u32;
                    section_masks
                        .get(q_idx)
                        .map(Vec::as_slice)
                        .unwrap_or_default()
                        .iter()
                        .filter(|&&(left, right)| left < right)
                        .map(|&(left, right)| (left as u32 + offset, right as u32 - 1 + offset))
                        .collect()
                })
                .collect(),
            code: crate::algorithm::tblastx::report::DisplayCode::new(&code),
        })
    } else {
        None
    };
    // The warnings of each query are written before its report, once, with the first
    // format (`QueryWarnings`); the taxonomy database's messages of XML2 and JSON go to the
    // diagnostics at every format.
    let mut warning_lines = Some(query_posts);
    // The queries' result warnings that a structured first format does not write are each later
    // format's that writes them, to the caller's receiver of that format's
    // (`ReportOutputs::later_warnings`).
    // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    let keeps_later = outputs.keeps_later_query_warnings();
    let observer = &mut outputs.observer;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
    // The statistics of each query for the observer, the values of the XML reports.
    if let Some(observer) = observer.as_deref_mut() {
        let stats =
            super::blast_xml::query_stats(parameters, query_validity, scoring, ungapped_karlin);
        for (q_idx, stats) in stats.into_iter().enumerate() {
            observer.query_stats(q_idx, stats);
        }
    }
    // What the archive reads of the report's data.
    let archive_data = || super::archive::ArchiveData {
        queries: query_records,
        subject_records,
        subject_letters,
        shown,
        hitlists,
        parameters,
        ungapped_karlin,
        query_validity,
        query_batch_skipped,
        gapped: scoring.gapped.then(|| {
            lookup_protein_params(&ProteinScoringSpec {
                matrix: scoring.matrix,
                gap_open: scoring.gap_open,
                gap_extend: scoring.gap_extend,
            })
        }),
        seg,
        mask_lowercase,
        ranges,
    };
    for (format_index, format) in outputs.formats.iter_mut().enumerate() {
        // A format that failed alone writes nothing more; the next takes the warnings.
        if format.refusal.refused() {
            continue;
        }
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a call with several formats, `run_local`).
        if let Some(handle) = archive.map(|render| render.handle) {
            handle.set_muted(format_number(format.outfmt) != ARCHIVE_FORMAT);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // A tabular report with a field list or a delimiter, or outfmt 10 or 20
        // (`tabular::tabular_spec`); `None` for the default 6 and 7.
        let custom_tabular = super::tabular::tabular_spec(format.outfmt);
        // So do the warnings of the formatters of a tabular field list that shows a taxonomy
        // name (`x_CheckTaxDB`, `tabular`), and the Seq-annots' posts of the bytes that NCBI's
        // stream replaces (`seqalign::write_query_to`).
        // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
        let later = format.later_warnings(keeps_later, &mut outputs.later_warnings, format_index);
        let mut warnings = (warning_lines.is_some()
            || later.is_some()
            || crate::report::query_warnings::report_posts(format.report_format())
            || custom_tabular
                .as_ref()
                .is_some_and(crate::report::tabular::TabularSpec::looks_up_every_taxid))
        .then(|| {
            let first = warning_lines.take();
            let warnings = crate::report::query_warnings::QueryWarnings::of_posts(
                first.unwrap_or(&[]),
                archive.map(|render| render.handle),
                &mut *outputs.diagnostics,
            )
            .with_later(query_posts, later);
            // The first format's writer also writes the texts of the subjects that the reports
            // read one by one (`query_reads`).
            match first {
                Some(_) => warnings.with_subject_reads(query_reads),
                None => warnings,
            }
        });
        // The caller's receiver of the side files of outfmt 13 and 14, if it gives one.
        let side_sink = format.sink.side_files().cloned();
        // The format's report, whose error goes to its refusal (below).
        let written = (|| -> Result<()> {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            match format_number(format.outfmt) {
                // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
                // (`tabular`); never the default 6 and 7 below.
                6 | 7 | 10 | 18 | 20 if custom_tabular.is_some() => {
                    super::tabular::write_tblastn_tabular(
                        &super::tabular::TblastnTabular {
                            spec: custom_tabular.as_ref().expect("a custom tabular report"),
                            inputs: tabular,
                            queries: query_records,
                            query_bioseq,
                            letters: subject_letters,
                            subject_ids: &subject_ids,
                            subject_records,
                            shown,
                            query_batch_skipped,
                            ranges,
                            genetic_code,
                            matrix: scoring.matrix,
                            program_version: super::tabular::program_version(pssm.is_some()),
                        },
                        hitlists,
                        &mut writer,
                        warnings.as_mut(),
                        probe.as_mut(),
                    )?
                }
                6 | 7 => write_tabular(
                    &mut writer,
                    program,
                    format_number(format.outfmt) == 7,
                    query_records,
                    subject_records,
                    subject_letters,
                    &subject_ids,
                    shown,
                    hitlists,
                    query_batch_skipped,
                    ranges,
                    probe.as_mut(),
                    warnings.as_mut(),
                )?,
                number @ 0..=4 => write_pairwise(
                    &mut writer,
                    pairwise_hits
                        .as_deref()
                        .expect("pairwise hits prepared for TBLASTN outfmt 0 to 4"),
                    query_records,
                    subject_records,
                    subject_letters,
                    shown,
                    subject_set,
                    parameters,
                    ungapped_karlin,
                    query_validity,
                    query_batch_skipped,
                    scoring,
                    matrix_name,
                    report_counts,
                    formatting,
                    ranges,
                    probe.as_mut(),
                    warnings.as_mut(),
                    pssm,
                    crate::report::pairwise::PairwiseView::of(
                        crate::blastinput::app::ReportFormat::from_number(number)
                            .expect("a format number 0 to 4"),
                    )
                    .expect("a pairwise view of formats 0 to 4"),
                    anchored.as_ref(),
                )?,
                // BLAST XML: one document per run, a header before the first query and a
                // trailer after the last (`blast_xml`).
                5 => xml
                    .as_ref()
                    .expect("XML prepared for TBLASTN outfmt 5")
                    .write(
                        &mut crate::report::blast_xml::BlastXmlWriter::new(0),
                        &mut writer,
                        query_records,
                        subject_letters,
                        hitlists,
                        genetic_code,
                        seg,
                        mask_lowercase,
                        ranges,
                        warnings.as_mut(),
                    )?,
                // The XML2 and JSON documents (outfmt 16 and 15, and 14 and 13 as side files next to
                // `out` with this output as the master file): the same records, written by
                // `Xml2Writer`.
                13..=16 => {
                    let kind = xml2_kind(format.outfmt).expect("a checked XML2 or JSON format");
                    let xml = xml
                        .as_ref()
                        .expect("XML prepared for TBLASTN outfmt 13 to 16");
                    xml.write(
                        &mut crate::report::blast_xml2::Xml2Writer::for_run(
                            kind,
                            xml.xml2_run(),
                            out,
                            side_sink.as_ref(),
                        ),
                        &mut writer,
                        query_records,
                        subject_letters,
                        hitlists,
                        genetic_code,
                        seg,
                        mask_lowercase,
                        ranges,
                        warnings.as_mut(),
                    )?;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
                // (`eAsnBinary` is the same pruning and `MSerial_AsnBinary`,
                // blast_format.cpp:699-705.) The Seq-annot of each query with alignments.
                number @ (8 | 9 | 12) => super::seqalign::write_tblastn_seqalign(
                    match number {
                        8 => crate::report::serial::SerialFormat::AsnText,
                        9 => crate::report::serial::SerialFormat::AsnBinary,
                        _ => crate::report::serial::SerialFormat::Json,
                    },
                    &mut writer,
                    query_records,
                    subject_records,
                    subject_letters,
                    shown,
                    hitlists,
                    ranges,
                    warnings
                        .as_mut()
                        .expect("a Seq-annot writes its posts to the diagnostics"),
                    scoring.gapped,
                )?,
                // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:294-296
                // One BLAST archive per query batch (`archive`).
                ARCHIVE_FORMAT => {
                    let render = archive.expect("a BLAST archive request has its handle");
                    super::archive::write_tblastn_archives(
                        render.handle,
                        render,
                        &archive_data(),
                        warnings
                            .as_mut()
                            .expect("the BLAST archive writes its posts to the diagnostics"),
                        &mut writer,
                    )?
                }
                _ => bail!("unsupported TBLASTN outfmt {}", format.outfmt),
            }
            writer.flush()?;
            Ok(())
        })();
        // An error of the format's report fails the call, or this format alone (`FormatRefusal`):
        // the formats after it are written then.
        format.refusal.check(written)?;
    }
    // The messages that `search` writes after the formats (the stop messages) are the archive's
    // again, whichever format came last (as in blastn and blastp).
    if let Some(handle) = archive.map(|render| render.handle) {
        handle.set_muted(false);
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1266-1338:
// CBlastTabularInfo::PrintNumProcessed appends the final query count.
#[allow(clippy::too_many_arguments)]
fn write_tabular(
    writer: &mut impl Write,
    // The program of outfmt 7's first line (`TBLASTN`, `PSITBLASTN`).
    program: &str,
    comments: bool,
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    shown: &ShownSubjects<'_>,
    hitlists: &[KappaResultHitList],
    query_batch_skipped: &[bool],
    ranges: &ReportRanges<'_>,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-224:
    // a skipped batch has a null align set for each query it contains.
    ensure!(
        query_batch_skipped.len() == query_records.len(),
        "TBLASTN batch validity count mismatch"
    );
    // The rows are printed in the order of the final HSP list, so a running count is
    // each row's HSP index.
    let mut hsp_index: HspIndex = 0;
    // What the reports do at each subject (`report::db_subject::HitLoads`; a subject set's
    // subjects are all loaded).
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let loads = match shown {
        ShownSubjects::Fasta { .. } => &every_loaded,
        ShownSubjects::Database { reported, .. } => &reported.loads,
    };
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
            writeln!(writer, "# {program} {}", crate::wording::REPORT_VERSION)?;
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1305-1308
            // The title's bytes (`write_outfmt7_query_line`).
            crate::report::outfmt6::write_outfmt7_query_line(writer, &query.query_line())?;
            writeln!(writer, "# Database: {}", shown.name())?;
            if hit_count > 0 {
                writeln!(writer, "# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score")?;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180,204-208:
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1266-1338:
            // if (align_set) { m_Ostream << "# " << num_hits << " hits found\n"; }
            // A skipped search has a null alignment set, so the count line
            // is omitted. An ordinary searched no-hit query retains it.
            if !search_skipped {
                writeln!(writer, "# {hit_count} hits found")?;
            }
        }
        // The query's formatter is new: its rows keep no subject before the first one that
        // loads it (`TabularRows`).
        let mut rows = crate::report::tabular::TabularRows::default();
        for list in hitlist.lists() {
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            let _subject = subject_records
                .get(oid)
                .context("TBLASTN subject OID out of range")?;
            // The subject's load throws where the query's batch has its prefetch off
            // (`report::db_subject::HitLoads`): each row is written as the catch leaves it
            // (`TabularRows::write_thrown_default`).
            let thrown = loads.thrown(oid).is_some();
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
            // The subject's first row (the default fields read it) loads it where the query's
            // batch did not prefetch it (`QueryWarnings::subject_read`).
            if !list.hsps.hsps.is_empty() && !thrown {
                if let Some(warnings) = warnings.as_deref_mut() {
                    warnings.subject_read(q_idx, oid as u32, &mut *writer)?;
                }
            }
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                if thrown {
                    if let Some(probe) = probe.as_mut() {
                        writer.flush()?;
                        probe.begin(hsp_index);
                    }
                    rows.write_thrown_default(
                        writer,
                        query.shown_id(),
                        |kept| crate::algorithm::common::tabular::printed_id(subject_ids, kept),
                        linked.evalue,
                        payload.bit_score,
                        &crate::report::outfmt6::OutputConfig::ncbi_compat(),
                    )?;
                    if let Some(probe) = probe.as_mut() {
                        writer.flush()?;
                        probe.end(hsp_index);
                    }
                    hsp_index += 1;
                    continue;
                }
                let hsp = &linked.hsp;
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1087-1105:
                // Stage D's local subject coordinate conversion is 1-offset.
                let subject_len = i32::try_from(subject_letters.length(oid))?;
                let (s_start, s_end) = subject_endpoints(hsp, subject_len);
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
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108:
                // One printed row is one HSP; the probe marks it without changing it.
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-497
                // The query ID is the first word of its title, or its local ID without a title
                // (`FastaRecord::shown_id`); the subject's is `tabular_subject_ids`'s. Both are
                // written as bytes.
                writer.write_all(query.shown_id())?;
                writer.write_all(b"\t")?;
                // A failed lookup's rows are refused before (`check_failed_lookups`).
                rows.found(oid);
                writer.write_all(crate::algorithm::common::tabular::printed_id(
                    subject_ids,
                    oid,
                ))?;
                writeln!(
                    writer,
                    "\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    percent,
                    payload.align_length,
                    payload.report_mismatches,
                    payload.gap_opens,
                    q_start,
                    q_end,
                    s_start,
                    s_end,
                    evalue,
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:986-993:
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
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:342
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
/// A subject interval has the strand both, so `RemapAlignToLoc` (seq_align_util.cpp:71-78)
/// moves its row by the interval's start.
pub(super) fn record_coordinates(
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

/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`), as BLASTN's and TBLASTX's; `None`
/// where `GenerateDefline` throws (`common::tabular::SubjectId`), whose rows NCBI writes with its catch's values
/// (`report::tabular::SubjectLookup::Failed`) or stops at (`check_failed_lookups`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
fn tabular_subject_ids(
    subjects: &[FastaRecord],
) -> Vec<crate::algorithm::common::tabular::SubjectId> {
    subjects
        .iter()
        .map(|subject| {
            // A nucleotide record's ids are never refused (`FastaRecord::subject_tabular_id`).
            subject
                .subject_tabular_id(false)
                .ok()
                .flatten()
                .map(Arc::from)
        })
        .collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1087-1105: the nucleotide coordinates of a
// translated frame, descending on a minus frame.
/// The 1-based subject coordinates of the rows of `hsp` (the pairwise, tabular and BLAST XML
/// reports) in a subject of `subject_len` nucleotides: those of its frame's codons, from the end
/// of the minus strand, before the records' offsets (`record_coordinates`).
pub(super) fn subject_endpoints(
    hsp: &super::search_gapped::GappedHsp,
    subject_len: i32,
) -> (i32, i32) {
    let frame = i32::from(hsp.frame);
    if frame > 0 {
        (3 * hsp.s_start + frame, 3 * hsp.s_end + frame - 1)
    } else {
        (
            subject_len - 3 * hsp.s_start + frame + 1,
            subject_len - 3 * hsp.s_end + frame + 2,
        )
    }
}

/// The values of the outfmt 0 row of `linked` with its payload `payload` in a subject of
/// `subject_len` nucleotides (`pairwise_hits`) that the report's sorts read (`HspScores::of`),
/// without making the row, for the model of the reads (`TblastnDatabase::subject_read_texts`):
/// the payload's bit score and columns, the linked E-value, the engine's identities (the row's
/// `seqalign_num_ident`) and the row's coordinates (`subject_endpoints`). The records' offsets of
/// a ranged query move every HSP of the query alike (a database subject has none), which no sort
/// reads, and are left out.
pub(super) fn hsp_scores(
    linked: &super::stage_d_linking::LinkedHsp,
    payload: &super::stage_d_results::KappaHspPayload,
    subject_len: usize,
) -> Result<crate::report::hit_sort::HspScores> {
    let hsp = &linked.hsp;
    let (s_start, s_end) = subject_endpoints(hsp, i32::try_from(subject_len)?);
    Ok(crate::report::hit_sort::HspScores::of_values(
        payload.bit_score,
        linked.evalue,
        usize::try_from(payload.num_ident)?,
        payload.align_length,
        (
            usize::try_from(hsp.q_start + 1)?,
            usize::try_from(hsp.q_end)?,
        ),
        (usize::try_from(s_start)?, usize::try_from(s_end)?),
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1490-1589;
// c++/src/objtools/align_format/showalign.cpp:1600-1625:
// AcknowledgeBlastQuery; x_DisplayDeflines; DisplaySeqalign;
// the formatted result contains full aligned protein strings and the
// translated subject nucleotide endpoints.
/// The final HSP list with rendered alignments, in the order of the tabular rows.
#[allow(clippy::too_many_arguments)]
fn pairwise_hits(
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    hitlists: &[KappaResultHitList],
    scoring: LocalStageDScoring,
    genetic_code: u8,
    seg: Option<&SegParams>,
    mask_lowercase: bool,
    ranges: &ReportRanges<'_>,
    // Whether the alignments show the query's masks (not a PSSM query's).
    query_masks: bool,
    // The PSSM query's Bioseq letters where they are not its record's (`PsiSearch::bioseq_letters`),
    // which the rows show (`QueryLetters`); a run of such a query has no other.
    query_bioseq: Option<&[u8]>,
    // With `Some`, the masks of each query's rows (half-open, in the searched letters; none for a
    // PSSM query) are pushed in query order, for the query-anchored sections of outfmt 1 to 4,
    // which show the same masks (`render`: SEG runs once per query).
    mut masks_out: Option<&mut Vec<Vec<(usize, usize)>>>,
) -> Result<Vec<PairwiseHit>> {
    let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
    let mut hits = Vec::new();
    let mut x_window = Vec::new();
    for (q_idx, hitlist) in hitlists.iter().enumerate() {
        let query = query_records[q_idx].seq();
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1541-1557:
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:294
        // A PSSM query's search (`CPsiBlastImpl::Run`, psiblast_impl.cpp:120-170) never sets the
        // filtered regions: its alignments show no masks.
        let report_masks = if query_masks {
            encode_tblastn_lookup_query(query, seg, mask_lowercase).seg_masks
        } else {
            Vec::new()
        };
        for list in hitlist.lists() {
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            let subject = subject_records
                .get(oid)
                .context("TBLASTN subject OID out of range")?;
            ensure!(
                list.payloads.len() == list.hsps.hsps.len(),
                "TBLASTN report payload count mismatch"
            );
            // The subject's letters, or (EXPERIMENT, LOSAT_X_TBNSTORE) its packed form, of which
            // the translations and the rows below read windows.
            let subject_seq = subject_letters.nt(oid);
            let mut translation = TargetTranslation::of(subject_seq, &code);
            for (linked, payload) in list.hsps.hsps.iter().zip(&list.payloads) {
                let hsp = &linked.hsp;
                let frame = i32::from(hsp.frame);
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1087-1105:
                // translated coordinates use three nucleotide letters per
                // residue and reverse strand endpoints run in descending order.
                let (s_start, s_end) = subject_endpoints(hsp, i32::try_from(subject_seq.len())?);
                let query_letters = QueryLetters::new(query, query_bioseq);
                let (qseq, sseq) = aligned_sequences(
                    query_letters,
                    subject_seq,
                    &mut x_window,
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
                let report_positives = crate::report::row_letters::display_positives(
                    qseq.as_bytes(),
                    sseq.as_bytes(),
                    scoring.matrix,
                    query_letters.case(),
                );
                // The rows are read from the searched letters; the report prints record
                // coordinates (`record_coordinates`).
                let searched = (hsp.q_start + 1, hsp.q_end, s_start, s_end);
                let (q_start, q_end, s_start, s_end) =
                    record_coordinates(ranges, q_idx, oid, searched)?;
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:411-422
                // The report's frame comes from the record coordinate of the subject row's
                // start and the record's length (the search's frame without
                // `-subject_loc`).
                let subject_length = ranges.subjects.length(oid, subject_seq.len());
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
                    // The title after its first word, as BLASTN's and TBLASTX's; the reports take the
                    // subject's title bytes by subject index (`write_tblastn_pairwise_report`).
                    subject_title: subject.title.iter().position(|&byte| byte == b' ').map(
                        |space| String::from_utf8_lossy(&subject.title[space + 1..]).into_owned(),
                    ),
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:331-342:
                    // eDontAdjustMatrix -> 0; eCompoScaleOldMatrix -> 1;
                    // every other adjusted matrix -> 2.
                    comp_adjust_method: Some(match payload.matrix_adjust_rule {
                        EMatrixAdjustRule::DontAdjustMatrix => 0,
                        EMatrixAdjustRule::CompoScaleOldMatrix => 1,
                        _ => 2,
                    }),
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/link_hsps.c:1765-1810 sets the HSP
                    // linked count; showalign.cpp:3595-3598 prints sum_n.
                    sum_n: Some(linked.num),
                    subject_stored_gaps: Vec::new(),
                    // A PSSM query's gap residues (`row_letters::stored_gap_columns`).
                    query_stored_gaps: stored_gap_columns(
                        query,
                        b'-',
                        usize::try_from(hsp.q_start)?,
                        &payload.edit_script,
                        RowSide::Query,
                    ),
                    query_case: query_letters.case(),
                    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:519-531
                    // The Seq-align's `num_ident` (`tabular.rs`, the outfmt 8 / 11 / 12
                    // writers): the engine's identities, counted on the masked query, so that a
                    // SEG-masked letter is a mismatch there and an identity in the rows
                    // (`report_num_ident`).
                    seqalign_num_ident: usize::try_from(payload.num_ident)?,
                });
                // The model of the reads sorts by `hsp_scores` without the rows
                // (`TblastnDatabase::subject_read_texts`): the row's values where the records'
                // offsets move nothing.
                debug_assert!(
                    searched != (q_start, q_end, s_start, s_end)
                        || hsp_scores(linked, payload, subject_seq.len()).ok()
                            == hits.last().map(crate::report::hit_sort::HspScores::of),
                    "TBLASTN hsp_scores differs from its row"
                );
            }
        }
        if let Some(masks_out) = masks_out.as_deref_mut() {
            masks_out.push(report_masks);
        }
    }
    Ok(hits)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1490-1589;
// c++/src/objtools/align_format/showalign.cpp:1600-1625:
// AcknowledgeBlastQuery; x_DisplayDeflines; DisplaySeqalign;
#[allow(clippy::too_many_arguments)]
fn write_pairwise(
    writer: &mut impl Write,
    hits: &[PairwiseHit],
    query_records: &[FastaRecord],
    subject_records: &[FastaRecord],
    subject_letters: SubjectLetters<'_>,
    shown: &ShownSubjects<'_>,
    subject_set: SubjectSetNumbers,
    parameters: &LocalSubjectParameters,
    ungapped_karlin: &[crate::stats::tables::KarlinParams],
    query_validity: &[bool],
    query_batch_skipped: &[bool],
    scoring: LocalStageDScoring,
    matrix_name: &str,
    report_counts: crate::blastinput::report_counts::ReportCounts,
    formatting: crate::blastinput::report_counts::FormattingArgs,
    ranges: &ReportRanges<'_>,
    probe: Option<&mut FormatProbe<'_>>,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
    // The alignment section: outfmt 0's, or a query-anchored section (outfmt 1 to 4) of the
    // sequences `anchored`.
    view: crate::report::pairwise::PairwiseView,
    anchored: Option<&super::anchored::AnchoredSequences<'_>>,
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
    // The set's numbers, NCBI's `Int8`: a database's, or those of the subject records (the whole
    // set, each held in memory).
    let (num_sequences, total_nt) = match shown {
        ShownSubjects::Fasta { .. } => (
            subject_records.len() as i64,
            subject_letters
                .lengths()
                .map(|length| length as i64)
                .sum::<i64>(),
        ),
        // A database's lines print `pairwise` itself (NCBI's signed counts); these numbers serve
        // a FASTA set's lines only.
        ShownSubjects::Database { pairwise, .. } => {
            (pairwise.num_sequences, pairwise.total_letters)
        }
    };
    // The Gumbel block's database length takes `-dbsize` (`SubjectSetNumbers`). An ungapped
    // search has neither a Gumbel nor a gapped block (its gap costs need not be in the matrix's
    // table, `args.rs`): the footers print the ungapped block only.
    let gumbel = scoring
        .gapped
        .then(|| subject_set.gumbel_block(&spec, total_nt))
        .flatten();
    if scoring.gapped && gumbel.is_none() && !subject_set.old_fsc {
        anyhow::bail!("TBLASTN pairwise Spouge parameters are unavailable");
    }
    let report = BlastpPairwiseReport {
        version: crate::wording::REPORT_VERSION.into(),
        database_name: shown.name().into(),
        database_num_sequences: num_sequences,
        database_total_letters: total_nt,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2267: options.GetMatrixName(), as typed.
        matrix_name: Some(matrix_name.to_string()),
        gap_open: scoring.gap_open,
        gap_extend: scoring.gap_extend,
        word_threshold: scoring.threshold,
        window_size: scoring.window,
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:101-103
        // `kbp_gap`: a PSSM query's `kbp_gap_psi` (the PSSM's gapped values; the Gumbel values
        // stay the matrix's).
        gapped_karlin: match pssm {
            _ if !scoring.gapped => Default::default(),
            Some(psi) => psi.gapped_karlin(&lookup_protein_params(&spec)),
            None => lookup_protein_params(&spec),
        },
        gumbel,
        gapped: scoring.gapped,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:165-167
        // The formatter shows the size of each linked set (the `N` column) when the search is
        // ungapped and uses sum statistics (the link parameters are made with sum statistics,
        // `local_subject_parameters`).
        show_linked_set_size: !scoring.gapped && parameters.link.is_some(),
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2913-2927
        // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221);
        // `blastinput::report_counts` applies the -num_descriptions, -num_alignments and
        // -max_target_seqs values.
        num_descriptions: report_counts.num_descriptions,
        num_alignments: report_counts.num_alignments,
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
            // TBLASTN's report reads the query's title bytes by query index (`query_titles`).
            query_name: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:742-744
            // The query record's length, also with `-query_loc`.
            query_length: ranges.queries.length(q_idx, query.seq().len()),
            ungapped_karlin: karlin,
            effective_search_space: length.eff_searchsp,
        })
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
    let subject_titles: Vec<Arc<[u8]>> = subject_records
        .iter()
        .map(|record| Arc::from(record.title()))
        .collect();
    // NCBI reference (598d8ae6): c++/src/app/blast/tblastn_app.cpp:263-265
    let config = PairwiseConfig {
        program: "tblastn".into(),
        protein_matrix: scoring.matrix,
        line_length: formatting.line_length(),
        // `-sorthits` for every view, `-sorthsps` for outfmt 0 only (NCBI's warning for the others,
        // `PairwiseView::report_sort`).
        sort: view.report_sort(formatting.report_sort()),
        view,
        ..PairwiseConfig::default()
    };
    // A subject set's subjects are all loaded.
    let every_loaded = crate::report::db_subject::HitLoads::default();
    write_tblastn_pairwise_report(
        if pssm.is_some() {
            super::in_pssm::PROGRAM
        } else {
            "TBLASTN"
        },
        hits,
        writer,
        &config,
        &queries,
        &query_titles,
        query_validity,
        query_batch_skipped,
        &subject_titles,
        shown.parsed_subjects(),
        shown.subject_sources(),
        &report,
        shown.pairwise_database(),
        probe,
        warnings,
        ranges.epilog,
        anchored.map(|sequences| sequences as &dyn crate::report::query_anchored::ReportSequences),
        match shown {
            ShownSubjects::Fasta { .. } => &every_loaded,
            ShownSubjects::Database { reported, .. } => &reported.loads,
        },
    )
    // A subject whose load throws may end the run there (`x_pairwise_par::alignment_end`).
    .map_err(crate::cli::NativeError::of_writer)?;
    Ok(())
}

/// Whether the reports read other letters than the search, so that they count its HSPs again
/// from their rows (`recount_shown_counts`): the PSSM query `pssm` (`-in_pssm`) with
/// `QueryLetters::shows_other_letters` (its Bioseq's letters in lower case, or its gap residue),
/// or a search with `-db_hard_mask` (`hard_masked`, series NP, M4b PR 2), which read `B` over the
/// masked bases of the database sequences where the reports read their stored bases.
pub(super) fn recounts_shown_counts(
    pssm: Option<&crate::algorithm::pssm_search::PsiSearch>,
    hard_masked: bool,
) -> bool {
    hard_masked
        || pssm.is_some_and(|psi| {
            QueryLetters::new(psi.query.seq(), psi.bioseq_letters.as_deref()).shows_other_letters()
        })
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:996-1022 (cited at
// `blastp::blast_engine::build_pairwise_hits`) and c++/include/objtools/align_format/
// tabular.hpp:543-547
/// The report counts (`report_num_ident`, `report_mismatches`) of the HSPs
/// of `hitlists`, counted again from their rows (`aligned_sequences`) for a query whose reports
/// read other letters than its search (`QueryLetters::shows_other_letters`: a PSSM query's Bioseq
/// letters `query_bioseq`, or its gap residue), and for every query of a search with
/// `-db_hard_mask` (`hard_masked`: the search read `B` over the masked bases, the rows read the
/// stored bases; series NP, M4b PR 2). Stage D counts on the searched NCBIstdaa codes;
/// NCBI's reports compare the rows' letters, where a lower-case letter is identical to no subject
/// letter and the gap residue's `-` is identical to a gap of the subject (an identity and a gap,
/// so the mismatches can be fewer than none). Every report then reads these counts: the default
/// tabular rows, the custom rows without the alignment's rows, and outfmt 0's.
fn recount_shown_counts(
    query_records: &[FastaRecord],
    query_bioseq: Option<&[u8]>,
    subject_letters: SubjectLetters<'_>,
    hitlists: &mut [KappaResultHitList],
    genetic_code: u8,
    hard_masked: bool,
) -> Result<()> {
    let code = GeneticCode::try_from_id(genetic_code).map_err(anyhow::Error::msg)?;
    let mut x_window = Vec::new();
    for (query, hitlist) in query_records.iter().zip(hitlists.iter_mut()) {
        let letters = QueryLetters::new(query.seq(), query_bioseq);
        if !hard_masked && !letters.shows_other_letters() {
            continue;
        }
        for list in hitlist.lists_mut() {
            let oid = usize::try_from(list.oid).context("negative TBLASTN subject OID")?;
            let subject_seq = subject_letters.nt(oid);
            let mut translation = TargetTranslation::of(subject_seq, &code);
            for (linked, payload) in list.hsps.hsps.iter().zip(list.payloads.iter_mut()) {
                let hsp = &linked.hsp;
                let (qseq, sseq) = aligned_sequences(
                    letters,
                    subject_seq,
                    &mut x_window,
                    &code,
                    &mut translation,
                    hsp.frame,
                    hsp.q_start,
                    hsp.s_start,
                    hsp.s_end,
                    &payload.edit_script,
                    &[],
                )?;
                let ident = qseq
                    .bytes()
                    .zip(sseq.bytes())
                    .filter(|(q, s)| q == s)
                    .count();
                payload.report_num_ident = ident;
                payload.report_mismatches = i64::try_from(payload.align_length)?
                    - i64::try_from(ident + payload.gap_letters)?;
            }
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:971-1021;
// c++/src/objtools/align_format/showalign.cpp:2122-2149:
// the edit script consumes query and translated subject residues in order.
// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:911-918
// The subject row is translated again from the nucleotides for the display, with the
// table of `CTrans_table` (an ambiguous codon is B, Z or J where its residues are D/N, E/Q
// or I/L; TBLASTX's `display_residue`), not taken from the search's translation.
#[allow(clippy::too_many_arguments)]
fn aligned_sequences(
    // The query's letters as the reports read them (a PSSM query's Bioseq letters as they are,
    // a record's in upper case).
    query: QueryLetters<'_>,
    // EXPERIMENT (LOSAT_X_TBNSTORE): the letters, or a packed subject.
    subject: SubjectNt<'_>,
    window_buffer: &mut Vec<u8>,
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1508-1516
    // The nucleotides the row is translated from: the record's letters, or (EXPERIMENT,
    // LOSAT_X_TBNSTORE) the letters of the HSP's codons made from the store, the window NCBI
    // translates for the residues the script reads (from `s_start` on). Its codon of residue `o`
    // is the one `display_residue` reads in the whole subject: from `3 * (o - s_start) +
    // |frame| - 1` on the plus strand, and as far from the window's end as from the subject's on
    // the minus strand.
    let (nucleotides, first) = match subject {
        SubjectNt::Letters(letters) => (letters, 0),
        SubjectNt::Packed(_) => {
            let first = usize::try_from(s_start)?;
            let residues: usize = script
                .iter()
                .filter(|op| !matches!(op, GapEditOp::Ins(_)))
                .map(|op| op.num() as usize)
                .sum();
            let frame_shift = usize::from(frame.unsigned_abs());
            let length = 3 * residues + frame_shift - 1;
            let start = if frame > 0 {
                Some(3 * first)
            } else {
                (subject.len() + 1).checked_sub(3 * (first + residues) + frame_shift)
            }
            .filter(|start| start + length <= subject.len())
            .context("TBLASTN report HSP outside the subject")?;
            (subject.fetch(start, start + length, window_buffer), first)
        }
    };
    let display = crate::algorithm::tblastx::report::DisplayCode::new(code);
    let displayed_subject = |s: usize| -> u8 {
        display.residue(nucleotides, frame, (base + s as isize) as usize - first)
    };
    const NCBISTDAA_TO_AA: &[u8; 28] = b"-ABCDEFGHIKLMNPQRSTVWXYZU*OJ";
    // The rows live in a `PairwiseHit` while the reports are written: their exact size.
    let columns = script.iter().map(|op| op.num() as usize).sum();
    let mut query_string = String::with_capacity(columns);
    let mut subject_string = String::with_capacity(columns);
    for &op in script {
        for _ in 0..op.num() {
            match op {
                GapEditOp::Sub(_) => {
                    let qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    let &sa = translated
                        .get(s)
                        .context("TBLASTN report subject offset overflow")?;
                    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1541-1557 renders filter locations
                    // as lowercase in the pairwise query sequence.
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:184-187
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:945-949
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
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
                    let qa = query
                        .get(q)
                        .context("TBLASTN report query offset overflow")?;
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:184-187
                    // NCBI reference (598d8ae6): c++/src/objtools/readers/fasta.cpp:945-949
                    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2520-2521
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

/// `aligned_sequences` into two strings the caller keeps, for the outfmt 5 rows, which are
/// made for one HSP at a time. It calls `aligned_sequences` as it is (the rows of 0/6/7 are
/// made by that function alone, unchanged since the base of OF4b, DECISIONS 23) and moves its
/// strings into the caller's.
#[allow(clippy::too_many_arguments)]
pub(super) fn aligned_sequences_into(
    query: QueryLetters<'_>,
    subject: SubjectNt<'_>,
    window_buffer: &mut Vec<u8>,
    code: &GeneticCode,
    target: &mut TargetTranslation<'_>,
    frame: i8,
    q_start: i32,
    s_start: i32,
    s_end: i32,
    script: &[GapEditOp],
    report_masks: &[(usize, usize)],
    query_string: &mut String,
    subject_string: &mut String,
) -> Result<()> {
    (*query_string, *subject_string) = aligned_sequences(
        query,
        subject,
        window_buffer,
        code,
        target,
        frame,
        q_start,
        s_start,
        s_end,
        script,
        report_masks,
    )?;
    Ok(())
}

#[cfg(test)]
#[path = "stage_e_report_tests.rs"]
mod tests;
