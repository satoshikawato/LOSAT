//! The reports of a BLASTP call's final hits (`run_call`'s loop over the requested formats), one
//! function per arm of the loop's `match` over the format: the Seq-annots of outfmt 8, 9 and 12
//! (`write_seqalign`); the BLAST archive of outfmt 11 (`write_archive`); the structured reports,
//! BLAST XML (outfmt 5) and the XML2 and JSON documents (outfmt 13 to 16, `write_xml`); the
//! pairwise report of outfmt 0, the query-anchored reports of outfmt 1 to 4 and psiblast's round
//! reports (`write_pairwise`); the tabular reports with a field list or a delimiter, outfmt 10 and
//! 20 and the organism report of outfmt 18 (`write_custom_tabular`); and the default outfmt 6 and
//! 7 (`write_text`). Each function is the statements of one arm, and `run_call` calls it in that
//! arm: its error ends the loop as the arm's did.
//!
//! The NCBI references are at the functions: the Seq-annot branches of
//! `CBlastFormat::x_PrintStructuredReport` at `write_seqalign` and its BLAST XML, XML2 and JSON
//! branches at `write_xml`, `CBlastFormat::WriteArchive` at `write_archive`, the descriptions and
//! alignments of `CBlastFormat::PrintOneResultSet` at `write_pairwise`, the tabular writer of
//! `CBlastFormat::x_PrintTabularReport` at `write_custom_tabular`, and NCBI's default tabular
//! fields at `write_text`.
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`BlastpReports`, `BlastpRanges`, `BlastpReport`, `ShownSubjects`, `ShownQueries`,
//! `result_query_masks`, `set_row_masks`,
//! `write_blastp_tabular_hit_lists`) and its sibling modules `anchored`, `database`, `psi_round`,
//! `pssm_hooks`, `query_rows` and `x_soft_subject`.

use std::collections::BTreeSet;
use std::io::Write;
use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::{
    anchored, database, psi_round, pssm_hooks, query_record_masks, query_rows, result_query_masks,
    set_row_masks, write_blastp_tabular_hit_lists, x_soft_subject, BlastpRanges, BlastpReport,
    BlastpReports, ShownQueries, ShownSubjects,
};
use crate::algorithm::blastp::args::{BlastpApp, BlastpCompositionMode, ResolvedBlastpArgs};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::seqalign::EngineCounts;
use crate::algorithm::blastp::{archive, blast_xml, seqalign, tabular};
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::{FormatOutput, FormatProbe};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::query_batch::XmlFailingQuery;
use crate::report::blast_xml::XmlStats;
use crate::report::blast_xml2::SideFileSink;
use crate::report::db_subject::ParsedSubject;
use crate::report::defline::SourceTaxname;
use crate::report::pairwise::{PairwiseDatabase, PairwiseView};
use crate::report::query_warnings::QueryWarnings;
use crate::report::serial::SerialFormat;
use crate::report::tabular::TabularSpec;
use crate::report::{
    write_blastp_pairwise_window, BlastpPairwiseQuery, BlastpPairwiseReport, OutputFormat,
    PairwiseConfig, PairwiseHit, ReportContext,
};
use crate::stats::search_space::SearchSpace;
use crate::stats::{BlastGumbelBlk, KarlinParams};
use crate::utils::threading::SearchPool;

/// What the functions of this module read of `run_call`, made once before its loop over the
/// formats: each field is `run_call`'s value of the same name, by reference where it is not `Copy`
/// (`subject_records`, `subject_titles` and `shown_subjects` are the records, the titles and the
/// letters of the subjects that the reports show, a database's listed subjects' or the subject
/// set's, `report_prep::shown_of`; `subjects` are the encoded subjects). The format's output, its
/// warnings, its observer, the receiver of its side files and the subjects that the pairwise
/// report shows are the loop's own values, and `run_call`'s closures that an arm calls
/// (`first_context`, `batch_spans`, `rows_of`, `query_stats`) stay closures there: the functions
/// take each of them as a parameter.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:120-128
/// NCBI's formatter keeps the database's name, or the sequences of a subject set (bl2seq), from
/// its construction; the results of the search reach the formatter with each call
/// (blastp_app.cpp:281 and :287).
pub(super) struct ReportView<'a> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) reports: &'a BlastpReports<'a>,
    pub(super) ranges: &'a BlastpRanges<'a>,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) query_frames: &'a [Vec<QueryFrame>],
    pub(super) query_ids: &'a [Arc<[u8]>],
    pub(super) query_titles: &'a [Arc<[u8]>],
    pub(super) query_lengths: &'a [usize],
    pub(super) query_valid: &'a [bool],
    pub(super) query_batch_skipped: &'a [bool],
    pub(super) contexts: &'a [QueryContext],
    pub(super) report_contexts: &'a [QueryContext],
    pub(super) gapped_params: KarlinParams,
    pub(super) gapped_gumbel: Option<BlastGumbelBlk>,
    pub(super) search_spaces: &'a [SearchSpace],
    pub(super) exact_search_spaces: &'a [Option<i64>],
    pub(super) database: Option<&'a database::BlastpDatabase>,
    pub(super) listed_subjects: &'a Option<database::ListedSubjects>,
    pub(super) listed_titles: &'a Option<database::ListedTitles<'a>>,
    pub(super) listed_parsed: &'a [Option<Arc<ParsedSubject>>],
    pub(super) listed_sources: &'a [SourceTaxname],
    pub(super) pairwise_database: &'a Option<PairwiseDatabase>,
    pub(super) shown_database: &'a Option<PairwiseDatabase>,
    pub(super) subject_input_label: &'a [u8],
    /// NCBI's `Int8` (`report_prep::ReportLabels`).
    pub(super) shown_num_sequences: i64,
    pub(super) shown_total_letters: i64,
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) subject_titles: &'a [Arc<[u8]>],
    pub(super) shown_subjects: ShownSubjects<'a>,
    pub(super) subject_ids: &'a [crate::algorithm::common::tabular::SubjectId],
    pub(super) subject_accessions: &'a [Option<[Arc<[u8]>; 2]>],
    pub(super) subject_stitles: &'a [Arc<[u8]>],
    pub(super) hit_lists: &'a [Option<BlastpHitList>],
    /// What the XML2 and JSON reports do at each subject's `Hit`, by rank
    /// (`ListedTitles::hit_loads`, from `report_prep::prefetch_reported_subjects`).
    pub(super) hit_loads: &'a crate::report::db_subject::HitLoads,
    pub(super) engine_counts: &'a EngineCounts,
    pub(super) pairwise_hits: &'a Option<Arc<Vec<PairwiseHit>>>,
    pub(super) psiblast_hits_found: &'a Option<Vec<usize>>,
    pub(super) psi_found_before: &'a Option<psi_round::FoundBefore>,
    pub(super) rps_gapped_k: &'a Option<Vec<f64>>,
    /// Whether NCBI's result of each query has no alignment set: a query of a batch that NCBI did
    /// not search, but for rpsblast's search of several volumes, whose merge makes a set for every
    /// query (`blast_engine::no_alignment_sets`). Outfmt 7's hit count, the XML2 and JSON `hits`
    /// and message, and the organism report read it.
    pub(super) no_alignment_set: &'a [bool],
    pub(super) context: &'a ReportContext,
    pub(super) call_batches: &'a archive::CallBatches<'a>,
    pub(super) filtered_empty: bool,
    /// rpsblast over several volumes (`prelim_run::PrelimRun::rps_multi_volume`).
    pub(super) rps_multi_volume: bool,
    pub(super) xml_matrix_failure: Option<XmlFailingQuery>,
    pub(super) xml_only_xml2: bool,
    pub(super) x_rows_pool: Option<&'a SearchPool<'a>>,
    pub(super) blastp_parallel_pool: &'a SearchPool<'a>,
    /// The call ends with an exception before any query's report: the pairwise report writes its
    /// prolog alone (`empty_alignments::write_exception_prolog`).
    pub(super) prolog_only: bool,
}

/// Outfmt 8, 9 and 12: the Seq-annot of each query with alignments
/// (`seqalign::write_blastp_seqalign`), as ASN.1 text, ASN.1 binary or JSON, flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:700-706
/// The ASN.1 binary branch of `CBlastFormat::x_PrintStructuredReport`; the body quotes its ASN.1
/// text and JSON branch.
pub(super) fn write_seqalign(
    view: &ReportView<'_>,
    serial: SerialFormat,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
) -> Result<()> {
    let ReportView {
        args,
        ranges,
        query_records,
        query_lengths,
        database,
        listed_subjects,
        listed_parsed,
        pairwise_database,
        subject_records,
        hit_lists,
        engine_counts,
        ..
    } = *view;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
    // The Seq-annot of each query with alignments (`seqalign`).
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    let subjects = match (
        listed_subjects.as_ref(),
        pairwise_database.as_ref(),
        database,
    ) {
        (Some(listed), Some(pairwise), Some(database)) => {
            crate::report::record::SubjectNames::Database {
                ordinals: &listed.ordinals,
                oids: database.oid_map(),
                title: &pairwise.definition,
                parsed: listed_parsed,
            }
        }
        _ => crate::report::record::SubjectNames::Records(subject_records),
    };
    seqalign::write_blastp_seqalign(
        serial,
        // psiblast's own "Blast Type" in every round, and rpsblast's.
        args.app.program_type(),
        !args.ungapped,
        hit_lists,
        &crate::report::record::QueryNames {
            records: query_records,
            searched_lengths: query_lengths,
            placements: &ranges.query_input.placements,
        },
        subjects,
        ranges.subject_placements,
        engine_counts,
        &mut writer,
        warnings
            .as_mut()
            .expect("a Seq-annot writes its posts to the diagnostics"),
    )?;
    // Flushed when written, as the other reports are (a later batch's stop
    // messages follow it).
    writer.flush()?;
    Ok(())
}

/// Outfmt 11: one BLAST archive per query batch of the call (`archive::write_blastp_archives`),
/// with the batch's queries, their masks and the warnings of its queries, flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1067-1071
pub(super) fn write_archive(
    view: &ReportView<'_>,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    first_context: &impl Fn(usize) -> archive::FirstContext,
    batch_spans: &impl Fn() -> Vec<Range<usize>>,
) -> Result<()> {
    let ReportView {
        args,
        reports,
        ranges,
        query_frames,
        query_lengths,
        contexts,
        gapped_params,
        database,
        listed_subjects,
        listed_parsed,
        pairwise_database,
        subject_records,
        hit_lists,
        engine_counts,
        call_batches,
        ..
    } = *view;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:280-282
    // One BLAST archive per query batch of the call (`archive`).
    let handle = args
        .archive
        .as_ref()
        .expect("a BLAST archive request has its handle");
    // psiblast's archive of the call's round (`archive::PsiArchive`).
    let psi_karlin = args
        .psi_search
        .as_deref()
        .map(|psi| pssm_hooks::psi_karlin(psi, contexts, &gapped_params));
    let psi = match args.app {
        BlastpApp::Psiblast => Some(psi_archive(args, reports, database.is_some(), psi_karlin)?),
        _ => None,
    };
    // A PSSM search's blocks are its PSI blocks (`GetPsiUngappedKarlinBlk`,
    // `GetPsiGappedKarlinBlk`, build_archive.cpp:113-125), which a FASTA query's first round has as
    // its standard ones (`PsiSearch::ungapped_karlin`).
    let first_context = |first: usize| {
        let mut context = first_context(first);
        if let (Some(karlin), Some((ungapped, gapped))) = (psi_karlin, context.karlin.as_mut()) {
            *ungapped = karlin.ungapped;
            if let Some(gapped) = gapped.as_mut() {
                *gapped = karlin.gapped;
            }
        }
        context
    };
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    let subject_names = match (
        listed_subjects.as_ref(),
        pairwise_database.as_ref(),
        database,
    ) {
        (Some(listed), Some(pairwise), Some(database)) => {
            crate::report::record::SubjectNames::Database {
                ordinals: &listed.ordinals,
                oids: database.oid_map(),
                title: &pairwise.definition,
                parsed: listed_parsed,
            }
        }
        _ => crate::report::record::SubjectNames::Records(subject_records),
    };
    let masks = result_query_masks(args.app, query_frames, &ranges.query_input.placements);
    let spans = batch_spans();
    let call = archive::ArchiveCall {
        handle,
        psi: psi.as_ref(),
        gapped: !args.ungapped,
        program_type: args.app.program_type(),
        batches: archive::CallBatches {
            masks: &masks,
            spans: &spans,
            ..*call_batches
        },
        query_lengths,
        subjects: subject_names,
        subject_placements: ranges.subject_placements,
        hit_lists,
        engine_counts,
    };
    archive::write_blastp_archives(
        &call,
        warnings
            .as_mut()
            .expect("the BLAST archive writes its posts to the diagnostics"),
        &first_context,
        &mut writer,
    )?;
    // Flushed when written, as the other reports are.
    writer.flush()?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:259-271
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
// NCBI reference (598d8ae6): c++/src/algo/blast/api/search_strategy.cpp:525-561
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:150-168
/// psiblast's archive of the call's round (`WriteArchive` after the round's search, with its
/// number; `reports.psi_round`, 1 outside the round driver): the request of a PSSM search
/// (`-in_pssm`, `-in_msa`, a later round) holds the PSSM it searched (`PsiSearch::archive_pssm`)
/// with the round's `kbp_gap_psi[0]` and `kbp_psi[0]` written back (`karlin`), and the options
/// the composition-based mode of its set-up (`ArchiveRound::comp_based_stats`). The archive of a
/// PSSM search builds its request with the PSSM and the run's database, so with `-subject`
/// (`database` false) NCBI's `CExportStrategy` stops there with its engine error, after the
/// archives of the rounds before (the bl2seq archive of a FASTA query's first round is another
/// path, `x_Process_SearchDb` is not reached).
fn psi_archive(
    args: &ResolvedBlastpArgs,
    reports: &BlastpReports<'_>,
    database: bool,
    karlin: Option<pssm_hooks::PsiKarlin>,
) -> Result<archive::PsiArchive> {
    let number = reports.psi_round.map_or(1, |io| io.state.number() as i32);
    let pssm = match (args.psi_search.as_deref(), karlin) {
        (Some(psi), Some(karlin)) => {
            if !database {
                return Err(crate::blastinput::app::engine_error(
                    "Empty reference for CSearchDatabase.",
                ));
            }
            let object = psi
                .archive_pssm
                .as_ref()
                .expect("a PSSM search of a run that writes the archive keeps its PSSM");
            Some(Arc::new(object.with_final_karlin(
                crate::report::serial::pssm::FinalKarlin {
                    lambda: karlin.gapped.lambda,
                    kappa: karlin.gapped.k,
                    h: karlin.gapped.h,
                    lambda_ungapped: karlin.ungapped.lambda,
                    kappa_ungapped: karlin.ungapped.k,
                    h_ungapped: karlin.ungapped.h,
                },
            )))
        }
        _ => None,
    };
    let comp_based_stats = pssm.as_ref().map(|_| match args.comp_based_stats.mode {
        BlastpCompositionMode::NoCompositionBasedStats => 0,
        BlastpCompositionMode::CompositionBasedStats => 1,
        BlastpCompositionMode::CompositionMatrixAdjust => 2,
        BlastpCompositionMode::ForceFullMatrixAdjust => 3,
    });
    Ok(archive::PsiArchive {
        round: crate::report::archive::ArchiveRound {
            number,
            comp_based_stats,
        },
        pssm,
    })
}

/// Outfmt 5 (BLAST XML) and 13 to 16 (the XML2 and JSON documents): the structured report of the
/// queries (`blast_xml::write_blastp_xml`), the side files of outfmt 13 and 14 into the caller's
/// receiver where it gives one (`side_sink`). A matrix name that NCBI's XML does not know writes no
/// report data and fails where NCBI's run fails (`xml_matrix_failure`). Flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-710
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:745-749
/// The BLAST XML, XML2 and JSON branches of `CBlastFormat::x_PrintStructuredReport`.
pub(super) fn write_xml(
    view: &ReportView<'_>,
    outfmt: BlastpReport,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    side_sink: &Option<SideFileSink>,
    query_stats: &impl Fn() -> Vec<XmlStats>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    let ReportView {
        args,
        reports,
        ranges,
        query_records,
        query_frames,
        query_lengths,
        query_valid,
        query_batch_skipped,
        no_alignment_set,
        database,
        listed_subjects,
        listed_titles,
        listed_parsed,
        shown_database,
        subjects,
        subject_records,
        subject_titles,
        hit_lists,
        hit_loads,
        pairwise_hits,
        filtered_empty,
        rps_multi_volume,
        xml_matrix_failure,
        xml_only_xml2,
        ..
    } = *view;
    // The XML2 and JSON documents of outfmt 15 and 16 (`blast_xml2`) take the path
    // of outfmt 5: the same per-query records, written by another writer.
    let xml2_kind = match outfmt {
        BlastpReport::Xml2(kind) => Some(kind),
        _ => None,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1432-1444
    // The CBlastOutput document of the queries (`blast_xml`). NCBI makes its score
    // matrix table from the matrix name as typed, and fails for a name that is not
    // spelled as one of its own.
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218
    let positives = match crate::report::blast_xml::matrix::FormattingMatrix::new(&args.matrix_name)
    {
        Some(matrix) => matrix,
        None if query_records.is_empty() => {
            let mut writer = format.sink.open()?;
            if let Some(kind) = xml2_kind {
                // No report data is made without a query: the prolog and
                // the epilog alone.
                let mut bytes = Vec::new();
                crate::report::blast_xml2::Xml2Document::for_run(
                    kind,
                    false,
                    args.out.as_deref(),
                    side_sink.as_ref(),
                )
                .epilog(&mut bytes)?;
                writer.write_all(&bytes)?;
            }
            writer.flush()?;
            return Ok(());
        }
        None => {
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2362-2389
            // (cited at `blast_xml::matrix::XmlMatrixFailure`: the XML2 and JSON
            // report of a `-subject` or psiblast run is made when the next query's
            // results arrive, or in `PrintEpilog`)
            // The report data's constructor throws for the run's first searched query,
            // or where the report is held for its first of another Seq-id
            // (`blastp_xml_matrix_failure`, `QueryScopes::xml_failing_query`):
            // in the call of that query, after the reader's messages of the batches up
            // to its own and NCBI's posts of the bytes replaced in their search strategy
            // requests (`query_posts`); a call of the queries before it writes those
            // and no report. A run whose held report has no such query fails at
            // its end (the last call, with the epilog), unless the run's stop
            // (`QueryBatches::stop`, written by `search`) ends it first. The XML2 and
            // JSON prolog is written before the search (`PrintProlog`), once: with the
            // run's first call.
            let failure = xml_matrix_failure
                .expect("an XML format of a matrix name that NCBI's XML does not know");
            if let (Some(kind), true) = (xml2_kind, reports.x_prolog) {
                let mut writer = format.sink.open()?;
                writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                writer.flush()?;
            }
            let first = ranges.searched_base;
            let failing = failure
                .query()
                .and_then(|query| query.checked_sub(first))
                .filter(|&query| query < query_records.len());
            if let Some(warnings) = warnings.as_mut() {
                for query in 0..failing.map_or(query_records.len(), |query| query + 1) {
                    warnings.before_query(query, &mut std::io::sink())?;
                }
            }
            if failing.is_some()
                || (reports.last_call
                    && reports.epilog
                    && failure.fails_at_end(first + query_records.len()))
            {
                return Err(crate::report::blast_xml::matrix::unsupported_score_matrix());
            }
            return Ok(());
        }
    };
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    let xml_stats = query_stats();
    // The messages of each query's result: the same ones the warnings of the other
    // formats write (`query_warning_lines`).
    let xml_messages: Vec<Vec<String>> = query_records
        .iter()
        .enumerate()
        .map(|(q_idx, record)| {
            // rpsblast over several volumes keeps no warning in the merged results
            // (`query_warning_lines`).
            if rps_multi_volume {
                return Vec::new();
            }
            let mut messages = crate::report::query_warnings::query_messages(
                record.seq(),
                query_valid[q_idx],
                reports.x_queries_before + q_idx < reports.identity_warned,
                filtered_empty && !query_batch_skipped[q_idx],
            );
            // A PSSM search's (psiblast `-in_pssm`, a later round): those of
            // `PsiBlastSetupScoreBlock` too, as its warning line (`pssm_query_warning`).
            if let Some(psi) = args.psi_search.as_deref() {
                messages.extend(psi.setup_messages.iter().cloned());
            }
            messages
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:176-224
    // The result of a batch that NCBI did not search has no masks (rpsblast's merge of several
    // volumes takes volume 0's, rpsblast_local.cpp:174-175): XML2 and JSON list none for its
    // queries, which have an alignment set there (`no_alignment_set`).
    let xml_masks: Vec<Vec<crate::report::blast_xml::QueryMask>> =
        result_query_masks(args.app, query_frames, &ranges.query_input.placements)
            .into_iter()
            .zip(query_batch_skipped)
            .map(|(masks, &skipped)| {
                if skipped {
                    return Vec::new();
                }
                masks
                    .into_iter()
                    .map(|(from, to)| crate::report::blast_xml::QueryMask { from, to, frame: 0 })
                    .collect()
            })
            .collect();
    let database_length = |rank: usize| subjects[rank].aa_len;
    // The length of each record itself (`-subject_loc` searches part of it).
    let record_length = |s_idx: usize| {
        ranges
            .subject_placements
            .length(s_idx, subject_records[s_idx].seq().len())
    };
    // The taxids of the reported subjects, read for the XML2 and JSON formats only.
    let (xml_taxids, xml_taxonomy, xml_parsed_taxonomy, xml_taxdb_notes): (
        Vec<i32>,
        Vec<crate::report::blast_xml::HitTaxonomy>,
        Vec<Vec<crate::report::blast_xml::HitTaxonomy>>,
        crate::algorithm::common::database::TaxdbNotes,
    ) = match (xml2_kind, listed_titles.as_ref(), database) {
        (Some(_), Some(titles), Some(database)) => {
            titles.read_taxids(database, hit_lists, hit_loads)?
        }
        _ => (Vec::new(), Vec::new(), Vec::new(), Default::default()),
    };
    let xml_subjects = match (listed_subjects.as_ref(), listed_titles.as_ref(), database) {
        // The stored titles by rank, of the records or of LOSAT_X_BPLISTEDLEAN
        // (`subject_titles`, from `report_prep::shown_of`).
        (Some(listed), Some(_), Some(database)) => {
            crate::report::blast_xml::XmlSubjects::Database {
                ordinals: &listed.ordinals,
                oids: database.oid_map(),
                titles: crate::report::blast_xml::XmlTitles::Shared(subject_titles),
                // Only XML2 shows the taxid.
                taxids: &xml_taxids,
                taxonomy: &xml_taxonomy,
                length: &database_length,
                parsed: listed_parsed,
                parsed_taxonomy: &xml_parsed_taxonomy,
                loads: hit_loads,
            }
        }
        _ => crate::report::blast_xml::XmlSubjects::Records {
            protein: true,
            records: subject_records,
            length: &record_length,
        },
    };
    let query_names = crate::report::record::QueryNames {
        records: query_records,
        searched_lengths: query_lengths,
        placements: &ranges.query_input.placements,
    };
    let db_name = database.map(|database| database.name().to_vec());
    // `BlastOutput_program` and `_version` name the application: `blastp`, `rpsblast`
    // (`Blast_ProgramNameFromType` of its search's program, `RPSBLAST 2.17.0+`).
    let version = format!(
        "{} {}",
        args.report_program().to_ascii_uppercase(),
        crate::wording::REPORT_VERSION
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:414-417
    // `-soft_masking true` (`mask_at_hash`) ends the filter string with `m;`.
    let filter = crate::report::blast_xml::filter::filter_string(
        None,
        args.seg.filter_options(),
        None,
        None,
        args.soft_masking,
    );
    let xml = blast_xml::BlastpXml {
        run: crate::report::blast_xml::XmlRun {
            program: args.report_program(),
            app: args.app.ncbi_name(),
            version: &version,
            reference: match args.app {
                BlastpApp::Psiblast => crate::report::blast_xml::REFERENCE_COMP_BASED_STATS,
                _ => crate::report::blast_xml::REFERENCE_GAPPED_BLAST,
            },
            db_name: db_name.as_deref().unwrap_or(b""),
            matrix: Some(&args.matrix_name),
            expect: args.evalue,
            gap_open: args.scoring.gap_open,
            gap_extend: args.scoring.gap_extend,
            reward_penalty: None,
            filter: &filter,
            db_num: shown_database
                .as_ref()
                .map_or(0, |shown| shown.num_sequences),
            db_len: shown_database
                .as_ref()
                .map_or(0, |shown| shown.total_letters),
            db_len_overflow: shown_database
                .as_ref()
                .is_some_and(|shown| shown.total_letters_overflow),
            positives: Some(&positives),
            mask_char: b'X',
        },
        queries: &query_names,
        stats: xml_stats,
        messages: xml_messages,
        subjects: xml_subjects,
        taxdb_notes: xml_taxdb_notes,
        masks: xml_masks,
        // A query without NCBI's alignment set (`no_alignment_set`: of a batch that the search
        // skipped, every query in it invalid, but for the merge of several volumes): XML2 and
        // JSON write no `hits` and add to its message.
        searched: no_alignment_set.iter().map(|&none| !none).collect(),
        program: args.app.name(),
        program_type: args.app.program_type(),
        epilog: reports.epilog,
        stops: !reports.epilog && reports.last_call,
        only_xml2: xml_only_xml2,
        query_other_coding: args
            .psi_search
            .as_ref()
            .and_then(|search| search.query_other_coding),
        query_scopes: reports.query_scopes,
    };
    // The rows of the queries' HSPs (`query_rows`), one query at a time: the document
    // reads their letters query by query and formats them serially, so a window of
    // several queries would hold more rows for no gain (300 proteins x Swiss-Prot,
    // 8 threads: 46.0 MB with windows, 39.6 MB a query at a time, the same time).
    let mut xml_rows = query_rows::QueryRows::new(
        (0..hit_lists.len()).map(|q_idx| q_idx..q_idx + 1).collect(),
        pairwise_hits.as_deref().map(Vec::as_slice),
        rows_of,
    );
    match xml2_kind {
        None => blast_xml::write_blastp_xml(
            &xml,
            &mut crate::report::blast_xml::BlastXmlWriter::new(reports.x_queries_before),
            hit_lists,
            &mut xml_rows,
            &mut writer,
            warnings.as_mut(),
        )?,
        Some(kind) => {
            // gap costs only for a gapped search, `cbs` when the composition-based
            // statistics are on (blastxml2_format.cpp:531-546).
            let params = crate::report::blast_xml2::Xml2Run {
                target: if database.is_some() {
                    crate::report::blast_xml2::Xml2Target::Database
                } else {
                    crate::report::blast_xml2::Xml2Target::Subjects(subject_records.len())
                },
                subject_names: if database.is_some() {
                    None
                } else {
                    crate::report::blast_xml2::subject_names(
                        subject_records.iter().map(FastaRecord::believed),
                    )
                },
                gapped: !args.ungapped,
                cbs: match args.comp_based_stats.mode {
                    BlastpCompositionMode::NoCompositionBasedStats => 0,
                    BlastpCompositionMode::CompositionBasedStats => 1,
                    BlastpCompositionMode::CompositionMatrixAdjust => 2,
                    BlastpCompositionMode::ForceFullMatrixAdjust => 3,
                },
                query_gencode: None,
                db_gencode: None,
                // psiblast is an iterative search (`ProgramType::is_iterative`; its XML2 and
                // JSON documents are refused, `psiblast::args::writes_outfmt`).
                program_type: args.app.program_type(),
            };
            let mut document = crate::report::blast_xml2::Xml2Writer::for_run(
                kind,
                params,
                args.out.as_deref(),
                side_sink.as_ref(),
            );
            // EXPERIMENT (LOSAT_X_BPQBLOCK): a block continues the document
            // of the blocks before it and hands its state on.
            if let Some(carry) = reports.x_xml2_carry {
                document
                    .document_mut()
                    .resume(std::mem::take(&mut *carry.lock().expect("the carry")));
            }
            blast_xml::write_blastp_xml(
                &xml,
                &mut document,
                hit_lists,
                &mut xml_rows,
                &mut writer,
                warnings.as_mut(),
            )?;
            if let Some(carry) = reports.x_xml2_carry {
                *carry.lock().expect("the carry") = document.document_mut().take_carry();
            }
        }
    }
    // Flushed when written, as the other reports are.
    writer.flush()?;
    Ok(())
}

/// PHI-BLAST's outfmt 0 (`write_pairwise`): the call's searches and the rows of a query's hit
/// lists (`rows_of_lists` of `run_call`), which it makes for each occurrence.
#[derive(Clone, Copy)]
pub(super) struct PhiPairwise<'a> {
    pub(super) call: &'a super::phi_hooks::PhiCall<'a>,
    #[allow(clippy::type_complexity)]
    pub(super) rows_of_lists:
        &'a dyn Fn(&[Option<BlastpHitList>], usize) -> Result<Vec<PairwiseHit>>,
    /// Each query's messages as PHI-BLAST's report writes them before the query's report
    /// (`report_prep::phi_result_warnings`).
    pub(super) result_warnings: &'a [Vec<u8>],
}

impl<'a> PhiPairwise<'a> {
    /// The reports of PHI-BLAST's call `call` (`None` without a pattern) with its queries'
    /// messages `result_warnings` (`None` likewise).
    #[allow(clippy::type_complexity)]
    pub(super) fn of(
        call: Option<&'a super::phi_hooks::PhiCall<'a>>,
        rows_of_lists: &'a dyn Fn(&[Option<BlastpHitList>], usize) -> Result<Vec<PairwiseHit>>,
        result_warnings: Option<&'a [Vec<u8>]>,
    ) -> Option<Self> {
        call.zip(result_warnings)
            .map(|(call, result_warnings)| Self {
                call,
                rows_of_lists,
                result_warnings,
            })
    }
}

/// psiblast `-phi_pattern` (item PS8a): the report of format `outfmt` of PHI-BLAST's call
/// `phi`: outfmt 0 (`write_pairwise`), the custom tabular reports (`write_custom_tabular`), outfmt 6
/// and 7 (`write_phi_text`), and NCBI's engine error of outfmt 8, 9 and 12, whose Seq-annot writer
/// does not take the program (`ProgramNameToEnum`), at the first result with an alignment; a call
/// without one (no hit, or without NCBI's Karlin block, `PhiCall::no_statistics`) writes nothing,
/// and no error. LOSAT refuses the other formats with a pattern (`phi_hooks::phi_limits`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1658-1664
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-706
#[allow(clippy::too_many_arguments)]
pub(super) fn write_phi_report(
    view: &ReportView<'_>,
    outfmt: BlastpReport,
    phi: PhiPairwise<'_>,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
    pairwise_shown: &Option<BTreeSet<usize>>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    match outfmt {
        BlastpReport::Pairwise(pairwise) => write_pairwise(
            view,
            pairwise,
            format,
            warnings,
            probe,
            pairwise_shown,
            rows_of,
            Some(phi),
        ),
        BlastpReport::CustomTabular(spec) => {
            write_custom_tabular(view, spec, format, warnings, probe, rows_of, Some(phi))
        }
        BlastpReport::Text(text) => write_phi_text(
            view,
            phi.call,
            phi.result_warnings,
            text,
            format,
            warnings,
            probe,
        ),
        BlastpReport::Seqalign(_) => {
            let aligned = view.hit_lists.iter().flatten().any(|list| {
                list.hsplist_array
                    .iter()
                    .take(list.hsplist_count)
                    .any(|hsp_list| !hsp_list.hsps.is_empty())
            });
            match aligned {
                false => Ok(()),
                true => Err(crate::blastinput::app::engine_error(
                    "Program type 'phiblastp' not supported",
                )),
            }
        }
        BlastpReport::Xml | BlastpReport::Xml2(_) | BlastpReport::Archive => Err(anyhow::anyhow!(
            "LOSAT internal error: a PHI-BLAST report of a format that LOSAT refuses with -phi_pattern"
        )),
    }
}

/// Outfmt 0 to 4: the pairwise report of the queries (`write_blastp_pairwise_window`), or psiblast's
/// round report (`psiblast_round::write_psiblast_round_report`), a window of queries at a time
/// (`query_rows::row_windows`); outfmt 1 to 4 read the sequences and their masks
/// (`anchored::AnchoredSequences`). The subjects that it shows are the loop's `pairwise_shown`,
/// whose titles the loop has checked. `view` is the report's layout (`BlastpReport::Pairwise`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1539-1541
/// `CBlastFormat::PrintOneResultSet`'s descriptions, then its alignments
/// (`CDisplaySeqalign::DisplaySeqalign`, :1585) and the query's footer (`x_PrintOneQueryFooter`,
/// :1589).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_pairwise(
    report_view: &ReportView<'_>,
    view: PairwiseView,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
    pairwise_shown: &Option<BTreeSet<usize>>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
    phi: Option<PhiPairwise<'_>>,
) -> Result<()> {
    let ReportView {
        args,
        reports,
        ranges,
        query_records,
        query_frames,
        query_titles,
        query_lengths,
        query_valid,
        query_batch_skipped,
        contexts,
        report_contexts,
        gapped_params,
        gapped_gumbel,
        search_spaces,
        exact_search_spaces,
        database,
        listed_subjects,
        listed_parsed,
        listed_sources,
        pairwise_database,
        subject_input_label,
        shown_num_sequences,
        shown_total_letters,
        subjects,
        subject_records,
        subject_titles,
        shown_subjects,
        hit_lists,
        pairwise_hits,
        psi_found_before,
        rps_gapped_k,
        blastp_parallel_pool,
        hit_loads,
        prolog_only,
        ..
    } = *report_view;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68,71,90-92
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:187
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:246-248
    // (`args.line_length`: `-line_length`, else 60; `args.sort`: blastp's
    // `-sorthits` and `-sorthsps`, none for psiblast and rpsblast.)
    // Outfmt 1 to 4 sort by `-sorthits` as outfmt 0 and read no `-sorthsps`
    // (`PairwiseView::report_sort`).
    let pairwise_config = PairwiseConfig {
        line_length: args.line_length,
        show_frame: false,
        program: args.report_program().to_string(),
        protein_matrix: args.scoring.matrix,
        sort: view.report_sort(args.sort),
        view,
        html: args.html.as_ref().map(|links| {
            crate::report::html_ids::HtmlReport::new(
                links,
                listed_subjects.as_ref().zip(database).into_iter().flat_map(
                    |(listed, database)| {
                        listed
                            .ordinals
                            .iter()
                            .map(move |&o| database.oid_map().oid(o))
                    },
                ),
            )
        }),
        // psiblast's PSSM query of the molecule `other` (`PsiSearch::query_other_coding`).
        unreadable_query: args
            .psi_search
            .as_ref()
            .and_then(|search| search.query_other_coding)
            .map(|_| crate::report::pairwise::UnreadableQuery {
                oids: listed_subjects
                    .as_ref()
                    .zip(database)
                    .map(|(listed, database)| {
                        let oids = database.oid_map();
                        listed
                            .ordinals
                            .iter()
                            .map(|&ordinal| oids.oid(ordinal).get())
                            .collect()
                    }),
            }),
    };
    // psiblast `-in_pssm`: the footer prints the PSSM's `kbp_psi` and `kbp_gap_psi`
    // (the file's values).
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:451-467
    let psi_report_karlin = args
        .psi_search
        .as_deref()
        .map(|psi| pssm_hooks::psi_karlin(psi, contexts, &gapped_params));
    // The titles are written as bytes by query index (`query_titles`);
    // `query_name` stays empty (`BlastpPairwiseQuery` is BLASTX's too).
    let mut pairwise_queries: Vec<BlastpPairwiseQuery> = query_lengths
        .iter()
        .enumerate()
        .map(|(q_idx, &query_length)| BlastpPairwiseQuery {
            query_name: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:742-744
            // The query record's length, also with `-query_loc`.
            query_length: ranges.query_input.placements.length(q_idx, query_length),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:87,100
            ungapped_karlin: psi_report_karlin
                .map_or(contexts[q_idx].karlin_params, |karlin| karlin.ungapped),
            // NCBI's `Int8` (`exact_search_spaces`; `None` only where the
            // subjects have no letters, a space of 0).
            effective_search_space: exact_search_spaces[q_idx]
                .unwrap_or(search_spaces[q_idx].effective_space.round() as i64),
            valid: query_valid[q_idx],
            batch_skipped: query_batch_skipped[q_idx],
        })
        .collect();
    let mut pairwise_report = BlastpPairwiseReport {
        version: crate::wording::REPORT_VERSION.to_string(),
        database_name: crate::report::pairwise::scan_mode_definition_bytes(subject_input_label),
        database_num_sequences: shown_num_sequences,
        database_total_letters: shown_total_letters,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2267
        // Matrix name from `GetMatrixName`, as typed.
        matrix_name: Some(args.matrix_name.clone()),
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
        word_threshold: args.threshold,
        window_size: args.window_size,
        gapped_karlin: psi_report_karlin.map_or(gapped_params, |karlin| karlin.gapped),
        gumbel: gapped_gumbel,
        gapped: !args.ungapped,
        show_linked_set_size: false,
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2927: 500 descriptions and 250 alignments,
        // or the -num_descriptions / -num_alignments / -max_target_seqs values
        // (`ReportCounts`).
        num_descriptions: args.num_descriptions,
        num_alignments: args.num_alignments,
    };
    // PHI-BLAST: the footers' Karlin blocks and search spaces are the pattern's.
    if let Some(phi) = phi.as_ref() {
        super::phi_hooks::phi_footers(phi.call, &mut pairwise_queries, &mut pairwise_report);
    }
    // The shown subjects, whose titles the loop over the formats checks (`pairwise_shown`).
    let shown = pairwise_shown
        .as_ref()
        .expect("made for every pairwise format");
    // psiblast's two-part table looks each shown subject up in the previous round's
    // set by its handle: LOSAT refuses one that it cannot place as NCBI does.
    if let (Some(_), Some(found)) = (reports.psi_round, psi_found_before.as_ref()) {
        found.refuse_unsure_shown(shown.iter().copied())?;
    }
    // The results' query masks (`result_query_masks`: none for psiblast's rounds), which the rows
    // show in lower case.
    let query_masks = {
        let masks = result_query_masks(args.app, query_frames, &ranges.query_input.placements);
        (!masks.iter().all(Vec::is_empty)).then_some(masks)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // The windows of queries whose rows are made and written together
    // (`query_rows::row_windows`); one window of every query when the rows were
    // made up front (`rows_up_front`). Without a query, the prolog and the epilog.
    // EXPERIMENT (LOSAT_X_PAIRWISEPAR): the queries' parts formatted on the
    // search's pool after the search (`report::x_pairwise_par`); with
    // LOSAT_X_BPLKPAHEAD the next block's lookup build may still run on it.
    let x_pool = crate::report::x_pairwise_par::pool(blastp_parallel_pool);
    let mut windows = match pairwise_hits {
        _ if prolog_only => Vec::new(),
        Some(_) => vec![0..pairwise_queries.len()],
        None => query_rows::row_windows(
            hit_lists,
            pairwise_queries.len(),
            x_pool,
            pairwise_config.line_length,
        ),
    };
    if windows.is_empty() {
        windows.push(0..0);
    }
    for window in windows {
        // PHI-BLAST: each batch's report of its occurrences' rows (`phi_round`); the epilog
        // after the loop. Its results carry the query's masks (each occurrence's
        // `GetMaskedQueryRegions`, shown in lower case: `query_record_masks`), unlike psiblast's
        // other rounds (`result_query_masks`), and no subject's masks: NCBI's PHI-BLAST report
        // sets none (no `SetSubjectMasks`, so `-lcase_masking` shows a subject's letters in upper
        // case).
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2155-2187
        if let Some(phi) = phi.as_ref() {
            let phi_query_masks = query_record_masks(query_frames, &ranges.query_input.placements);
            let phi_query_masks =
                (!phi_query_masks.iter().all(Vec::is_empty)).then_some(phi_query_masks);
            use crate::report::pairwise::psiblast_round::phi_round;
            for q_idx in window {
                // A batch's report is headed by its first query (`PhiBatch`).
                let Some(batch) = phi.call.batch_at(q_idx) else {
                    continue;
                };
                // Without NCBI's Karlin block no query is searched (`PhiCall::no_statistics`):
                // the report has no pattern block and no hits.
                let search = phi.call.search_of(batch);
                let occurrences = search.map_or(&[][..], |search| &search.occurrences[..]);
                let occ_lists = phi.call.occurrence_lists(batch.searched);
                let mut occurrence_hits = Vec::with_capacity(occurrences.len());
                for occurrence in 0..occurrences.len() {
                    let mut rows = match occ_lists.get(occurrence) {
                        Some(slot @ Some(_)) => {
                            (phi.rows_of_lists)(std::slice::from_ref(slot), batch.searched)?
                        }
                        _ => Vec::new(),
                    };
                    // The rows' query masks as the other reports' (`set_row_masks`); psiblast
                    // refuses `-html`, so no spans are kept.
                    if let Some(query_masks) = phi_query_masks.as_deref() {
                        if !rows.is_empty() {
                            let no_masks = Default::default();
                            set_row_masks(
                                &mut rows,
                                Some(query_masks),
                                [&no_masks, &no_masks],
                                None,
                            );
                        }
                    }
                    occurrence_hits.push(rows);
                }
                let info = search.map(|search| phi_round::PhiQueryInfo {
                    text: &search.pattern.text,
                    probability: search.pattern.probability,
                    occurrences: occurrences
                        .iter()
                        .map(|occurrence| (occurrence.offset as u32, occurrence.length as u32))
                        .collect(),
                });
                phi_round::write_phi_round_report(
                    &occurrence_hits,
                    batch.queries.clone(),
                    &phi.call.result_warnings(batch, phi.result_warnings),
                    info.as_ref(),
                    phi.call.footer(batch),
                    &mut writer,
                    &pairwise_config,
                    &pairwise_queries,
                    query_titles,
                    subject_titles,
                    listed_parsed,
                    listed_sources,
                    hit_loads,
                    &pairwise_report,
                    pairwise_database.as_ref(),
                    probe.as_mut(),
                    warnings.as_mut(),
                )
                .map_err(crate::cli::NativeError::of_writer)?;
            }
            continue;
        }
        let mut rows: std::borrow::Cow<'_, [PairwiseHit]> = match pairwise_hits.as_deref() {
            Some(hits) => std::borrow::Cow::Borrowed(hits.as_slice()),
            None => std::borrow::Cow::Owned(rows_of(window.clone())?),
        };
        // `-db_soft_mask` and `-db_hard_mask`: the subjects' stored masks that meet the
        // query's HSPs with them, by NCBI's Seq-id, shown in lower case as `-lcase_masking`
        // shows a FASTA subject's (`db_mask_search::shown_db_subject_masks`).
        let db_subject_masks = match (
            args.db_mask_resolved.as_ref(),
            database,
            listed_subjects.as_ref(),
        ) {
            (Some(mask), Some(database), Some(listed)) => {
                use crate::algorithm::common::db_mask_search as db_mask;
                let ordinal = |rank: u32| listed.ordinals.get(rank as usize).copied();
                let stored = db_mask::read_shown_stored(
                    &rows,
                    database.seqdb(),
                    mask.id,
                    args.app.name(),
                    ordinal,
                )?;
                Some(db_mask::shown_db_subject_masks(
                    &rows,
                    |rank| stored.get(&rank).map(Vec::as_slice),
                    |rank| {
                        db_mask::subject_id_keys(
                            listed_parsed.get(rank as usize).and_then(Option::as_deref),
                            database.oid_map().oid(listed.ordinals[rank as usize]),
                        )
                    },
                ))
            }
            _ => None,
        };
        // Outfmt 1 to 4: the sections read the sequences and their masks
        // (`anchored`), not the rows' letters.
        let anchored_sequences =
            (view != crate::report::pairwise::PairwiseView::Pairwise).then(|| {
                anchored::AnchoredSequences {
                    contexts: report_contexts,
                    queries: ShownQueries::of(query_records, args, None),
                    query_placements: &ranges.query_input.placements,
                    subjects,
                    shown_subjects,
                    subject_records,
                    subject_placements: ranges.subject_placements,
                    parsed: listed_parsed,
                    listed_ordinals: listed_subjects
                        .as_ref()
                        .zip(database)
                        .map(|(listed, database)| (&listed.ordinals[..], database.oid_map())),
                    show_gis: args.show_gis,
                    query_masks: query_masks.as_deref(),
                    subject_masks: if args.lcase_masking {
                        x_soft_subject::x_subject_masks(
                            &rows,
                            subject_records,
                            ranges.subject_placements,
                        )
                    } else {
                        db_subject_masks.clone().unwrap_or_default()
                    },
                }
            });
        // Under `-html` the rows' masks of each HSP (`HtmlRowMasks`, the spans of outfmt 0).
        let mut html_masks = None;
        if anchored_sequences.is_none()
            && (query_masks.is_some() || args.lcase_masking || db_subject_masks.is_some())
        {
            // `-lcase_masking`: the subject rows show their records' lower case too; a masked
            // database search, its masks.
            let lcase_masks = if args.lcase_masking {
                x_soft_subject::x_subject_masks(&rows, subject_records, ranges.subject_placements)
            } else {
                Default::default()
            };
            let no_masks = Default::default();
            html_masks = args.html.is_some().then(Vec::new);
            set_row_masks(
                rows.to_mut(),
                query_masks.as_deref(),
                [&lcase_masks, db_subject_masks.as_ref().unwrap_or(&no_masks)],
                html_masks.as_mut(),
            );
        }
        // psiblast's rounds: the round's report (`psiblast_round`), without the
        // epilog.
        if let Some(io) = reports.psi_round {
            crate::report::pairwise::psiblast_round::write_psiblast_round_report(
                &rows,
                window,
                &mut writer,
                &pairwise_config,
                &pairwise_queries,
                query_titles,
                subject_titles,
                listed_parsed,
                listed_sources,
                hit_loads,
                &pairwise_report,
                pairwise_database.as_ref(),
                probe.as_mut(),
                warnings.as_mut(),
                crate::report::pairwise::psiblast_round::ReportRound {
                    number: io.state.number(),
                    found_before: psi_found_before
                        .as_ref()
                        .map(|found| found.found.as_slice()),
                },
                anchored_sequences.as_ref().map(|sequences| {
                    sequences as &dyn crate::report::query_anchored::ReportSequences
                }),
            )
            .map_err(crate::cli::NativeError::of_writer)?;
            continue;
        }
        write_blastp_pairwise_window(
            &rows,
            window,
            &mut writer,
            &pairwise_config,
            &pairwise_queries,
            query_titles,
            subject_titles,
            listed_parsed,
            listed_sources,
            &pairwise_report,
            // rpsblast: the gapped K of the RPS traceback (`rps_hooks::rps_search`).
            rps_gapped_k.as_deref(),
            pairwise_database.as_ref(),
            probe.as_mut(),
            warnings.as_mut(),
            reports.epilog && !prolog_only,
            x_pool,
            anchored_sequences
                .as_ref()
                .map(|sequences| sequences as &dyn crate::report::query_anchored::ReportSequences),
            html_masks.as_deref(),
            hit_loads,
        )
        // A subject whose load throws may end the run there (`x_pairwise_par::alignment_end`).
        .map_err(crate::cli::NativeError::of_writer)?;
    }
    // PHI-BLAST: psiblast's epilog after the last query (`PrintEpilog`; `-phi_pattern` refuses
    // `-html` in `phi_hooks::phi_limits`, so the epilog is plain).
    if phi.is_some() && reports.epilog && !prolog_only {
        crate::report::pairwise::psiblast_round::write_psiblast_epilog(
            &mut writer,
            &pairwise_report,
            pairwise_database.as_ref(),
            false,
        )
        .map_err(crate::cli::NativeError::of_writer)?;
    }
    Ok(())
}

/// The tabular reports with a field list or a delimiter and outfmt 10 and 20
/// (`tabular::write_blastp_tabular`), and the organism report of outfmt 18
/// (`tabular::write_blastp_tax_report`), into the format's stream: the taxonomy and the names that
/// their fields show are read before the report's first byte, and the rows of the HSPs a window of
/// queries at a time where a field shows them (`rows_of`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:782-785
///
/// The organism report of outfmt 18 is NCBI's `CBlastFormat::x_PrintTaxReport`:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:873-874
pub(super) fn write_custom_tabular(
    view: &ReportView<'_>,
    spec: &'static TabularSpec,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
    phi: Option<PhiPairwise<'_>>,
) -> Result<()> {
    let ReportView {
        args,
        reports,
        ranges,
        query_records,
        query_ids,
        query_titles,
        query_lengths,
        database,
        listed_subjects,
        listed_titles,
        listed_parsed,
        listed_sources,
        subjects,
        subject_records,
        subject_ids,
        subject_accessions,
        subject_stitles,
        hit_lists,
        engine_counts,
        pairwise_hits,
        psiblast_hits_found,
        no_alignment_set,
        context,
        x_rows_pool,
        ..
    } = *view;
    // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
    // (`tabular`); never the default 6 and 7 below. The format's stream as the sink
    // opens it (`cli::NcbiOstream::sink`): the organism report's.
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    let needs = spec.needs();
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-842
    // `slen` is the length of the subject's Bioseq: the whole record with
    // `-subject_loc`, as the outfmt 0 `Length=` line.
    let subject_length = |s_idx: usize| {
        ranges.subject_placements.length(
            s_idx,
            subjects.get(s_idx).map_or(0, |subject| subject.aa_len),
        )
    };
    // The taxids of the reported database sequences where a field shows them, and
    // the taxonomy database's names that the rows look up (read before the report's
    // first byte, so that their refusals come before it).
    let taxonomy = match (listed_titles.as_ref(), database) {
        (Some(titles), Some(database)) if needs.taxonomy => {
            titles.tabular_taxonomy(database, hit_lists, view.hit_loads)?
        }
        _ => Vec::new(),
    };
    // A `-subject` search's organism report (outfmt 18) looks names up in the run's
    // database.
    let names = match database {
        Some(database) => database.tabular_names(spec, &taxonomy, args.app.name())?,
        None => args
            .subject_taxdb
            .as_ref()
            .map(|taxdb| taxdb.report_names(args.app.name())),
    };
    let tabular_subjects = match (listed_subjects.as_ref(), database) {
        (Some(listed), Some(database)) => tabular::TabularSubjects::Database {
            ids: subject_ids,
            accessions: subject_accessions,
            ordinals: &listed.ordinals,
            oids: database.oid_map(),
            titles: subject_stitles,
            parsed: listed_parsed,
            sources: listed_sources,
            length: &subject_length,
            taxonomy: &taxonomy,
            loads: view.hit_loads,
            seen: listed_titles
                .as_ref()
                .map_or_else(Default::default, database::ListedTitles::seen),
        },
        _ => tabular::TabularSubjects::Records {
            records: subject_records,
            ids: subject_ids,
            accessions: subject_accessions,
            parsed: listed_parsed,
            length: &subject_length,
        },
    };
    let program_version = format!(
        "{} {}",
        args.report_program().to_uppercase(),
        crate::wording::REPORT_VERSION
    );
    let tab = tabular::BlastpTabular {
        spec,
        app: args.app,
        program_version: &program_version,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1458
        // psiblast's round (`psi_round`), 1 outside the rounds.
        iteration: (args.app == BlastpApp::Psiblast)
            .then(|| reports.psi_round.map_or(1, |io| io.state.number())),
        matrix: args.scoring.matrix,
        query_ids,
        query_records,
        query_titles,
        searched_lengths: query_lengths,
        range_length: ranges.query_input.range_length,
        ungapped: args.ungapped,
        database: context
            .database_name
            .as_deref()
            .or(context.subject_name.as_deref())
            .unwrap_or_default(),
        subjects: tabular_subjects,
        unsearched: no_alignment_set,
        hits_found: psiblast_hits_found.as_deref(),
        prolog: reports.x_prolog,
        epilog: reports
            .epilog
            .then_some(reports.x_queries_before + query_titles.len()),
        taxdb_found: match database {
            Some(database) => database.taxdb_found(),
            None => args.subject_taxdb_found,
        },
        names: names.as_ref(),
        engine_counts,
        query_other_coding: args
            .psi_search
            .as_ref()
            .and_then(|search| search.query_other_coding),
    };
    // Outfmt 18, the organism report: the `Query=` block shows the query record's
    // length, also with `-query_loc` (as outfmt 0's, `pairwise_queries`).
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        let record_lengths: Vec<usize> = query_lengths
            .iter()
            .enumerate()
            .map(|(q_idx, &length)| ranges.query_input.placements.length(q_idx, length))
            .collect();
        tabular::write_blastp_tax_report(
            &tab,
            &record_lengths,
            hit_lists,
            &mut writer,
            warnings.as_mut(),
            args.is_html,
        )?;
        writer.flush()?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2039-2057
    // PHI-BLAST (`-phi_pattern`): each batch's messages (`PhiCall::result_warnings`), then the
    // tabular report of each of its occurrences in turn (`write_phi_text`), each as a report of
    // its own (the values of the batch's first query repeated for each of its occurrences, the
    // rows its query with letters', `PhiBatch`; without NCBI's Karlin block one result without a
    // Seq-align set); outfmt 20's line of the field names once, and outfmt 7's last line counts 0
    // queries. The texts of the subjects' loads are not written there (an open point of PS8a:
    // they are keyed by query).
    if let Some(phi) = phi.as_ref() {
        // Outfmt 20's line of the field names is the prolog's, before the first query's messages.
        if tab.prolog {
            tabular::write_tabular_prolog(spec, &mut writer)?;
        }
        for (batch_index, batch) in phi.call.batches.iter().enumerate() {
            if let Some(warnings) = warnings.as_mut() {
                for q_idx in batch.queries.clone() {
                    warnings.before_query(q_idx, &mut writer)?;
                }
            }
            writer.write_all(&phi.call.result_warnings(batch, phi.result_warnings))?;
            let searched = batch.searched;
            let result_queries = phi.call.result_queries(batch);
            let unsearched = phi.call.search_of(batch).is_none();
            let occ_lists = if unsearched {
                &[][..]
            } else {
                phi.call.occurrence_lists(searched)
            };
            let lists: Vec<Option<BlastpHitList>> = (0..result_queries.len())
                .map(|occurrence| occ_lists.get(occurrence).cloned().flatten())
                .collect();
            let results = lists.len();
            // The rows of a result are made from its query's values, then carry the result's
            // index.
            let entry_rows = |window: Range<usize>| -> Result<Vec<PairwiseHit>> {
                let mut all = Vec::new();
                for index in window {
                    let mut rows =
                        (phi.rows_of_lists)(std::slice::from_ref(&lists[index]), searched)?;
                    for row in &mut rows {
                        row.hit.q_idx = index as u32;
                    }
                    all.extend(rows);
                }
                Ok(all)
            };
            let mut rows = tabular::needs_rows(needs.rows, spec, args.app).then(|| {
                query_rows::QueryRows::new(
                    (0..results).map(|index| index..index + 1).collect(),
                    None,
                    &entry_rows,
                )
            });
            // The HSPs of a result carry its index too, which the rows' lookup compares
            // (`tabular::hsp_rows`).
            let shown: Vec<Option<BlastpHitList>> = lists
                .iter()
                .enumerate()
                .map(|(index, list)| {
                    list.clone().map(|mut list| {
                        for hsp_list in &mut list.hsplist_array {
                            for hsp in &mut hsp_list.hsps {
                                hsp.q_idx = index as u32;
                            }
                        }
                        list
                    })
                })
                .collect();
            let ids: Vec<_> = result_queries
                .iter()
                .map(|&q_idx| query_ids[q_idx].clone())
                .collect();
            let records: Vec<_> = result_queries
                .iter()
                .map(|&q_idx| query_records[q_idx].clone())
                .collect();
            let titles: Vec<_> = result_queries
                .iter()
                .map(|&q_idx| query_titles[q_idx].clone())
                .collect();
            let lengths: Vec<_> = result_queries
                .iter()
                .map(|&q_idx| query_lengths[q_idx])
                .collect();
            let unsearched = vec![unsearched; results];
            let last = batch_index + 1 == phi.call.batches.len();
            let phi_tab = tabular::BlastpTabular {
                query_ids: &ids,
                query_records: &records,
                query_titles: &titles,
                searched_lengths: &lengths,
                unsearched: &unsearched,
                hits_found: None,
                prolog: false,
                epilog: (reports.epilog && last).then_some(0),
                ..tab
            };
            tabular::write_blastp_tabular(
                &phi_tab,
                &shown,
                &mut writer,
                None,
                probe.as_mut(),
                rows.as_mut(),
            )?;
        }
        writer.flush()?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2955-2958
    // The rows of the queries' HSPs when a field shows them (and rpsblast's
    // identities, `tabular::needs_rows`), a window at a time (`query_rows`). A tabular
    // report has no lines: NCBI reads `-line_length` only for outfmt 0 to 4
    // (blast_args.cpp:2939-2941) and keeps its default otherwise (`kDfltLineLength`,
    // blast_args.hpp:1094), the line length that the windows' estimate takes here.
    let mut rows = tabular::needs_rows(needs.rows, spec, args.app).then(|| {
        query_rows::QueryRows::new(
            query_rows::row_windows(
                hit_lists,
                hit_lists.len(),
                x_rows_pool,
                crate::report::pairwise::DEFAULT_LINE_LENGTH,
            ),
            pairwise_hits.as_deref().map(Vec::as_slice),
            rows_of,
        )
    });
    tabular::write_blastp_tabular(
        &tab,
        hit_lists,
        &mut writer,
        warnings.as_mut(),
        probe.as_mut(),
        rows.as_mut(),
    )?;
    // Flushed when written, as the default tabular report is (below).
    writer.flush()?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2034-2067
/// psiblast `-phi_pattern` (item PS8a): outfmt 6 and 7 of the `std` fields, the tabular report
/// of each query occurrence of each batch in turn (`write_blastp_tabular_hit_lists` of the
/// occurrence's part of the list of the batch's query with letters, `PhiCall::occurrence_lists`,
/// under the id and title of the batch's first query, `PhiBatch`; outfmt 7 with its whole header
/// each time), then outfmt 7's last line, whose count of the formatted queries PHI-BLAST's report
/// leaves at 0. Each batch's messages come first (`PhiCall::result_warnings`). Without NCBI's
/// Karlin block nothing is searched (`PhiCall::no_statistics`): each batch has one result without
/// a Seq-align set, whose outfmt 7 header has no line of hits found.
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2039-2057
pub(super) fn write_phi_text(
    view: &ReportView<'_>,
    phi: &super::phi_hooks::PhiCall<'_>,
    result_warnings: &[Vec<u8>],
    outfmt: OutputFormat,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
) -> Result<()> {
    let ReportView {
        reports,
        query_ids,
        query_titles,
        subject_ids,
        subject_accessions,
        context,
        ..
    } = *view;
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    for batch in &phi.batches {
        if let Some(warnings) = warnings.as_mut() {
            for q_idx in batch.queries.clone() {
                warnings.before_query(q_idx, &mut writer)?;
            }
        }
        writer.write_all(&phi.result_warnings(batch, result_warnings))?;
        let no_set = phi.search_of(batch).is_none();
        let occ_lists = if no_set {
            &[][..]
        } else {
            phi.occurrence_lists(batch.searched)
        };
        for (occurrence, q_idx) in phi.result_queries(batch).into_iter().enumerate() {
            let rewritten;
            let part: &[Option<BlastpHitList>] = match occ_lists.get(occurrence) {
                Some(Some(list)) if batch.shows_gaps() => {
                    // The rows' query is the batch's first (its id), whose letters they show.
                    let mut copy = list.clone();
                    for hsp_list in &mut copy.hsplist_array {
                        for hsp in &mut hsp_list.hsps {
                            hsp.q_idx = q_idx as u32;
                        }
                    }
                    rewritten = [Some(copy)];
                    &rewritten
                }
                Some(slot) => std::slice::from_ref(slot),
                None => &[None],
            };
            write_blastp_tabular_hit_lists(
                part,
                outfmt,
                &mut writer,
                query_ids,
                &query_titles[q_idx..=q_idx],
                subject_ids,
                subject_accessions,
                view.hit_loads,
                &[no_set],
                None,
                context,
                probe.as_mut(),
                None,
                false,
                0,
                Some(1),
                None,
            )?;
        }
    }
    if reports.epilog && outfmt == OutputFormat::TabularWithComments {
        writeln!(writer, "# BLAST processed 0 queries")?;
    }
    writer.flush()?;
    Ok(())
}

/// The default outfmt 6 and 7 (`write_blastp_tabular_hit_lists`): NCBI's default fields, with the
/// comment lines of outfmt 7, flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
pub(super) fn write_text(
    view: &ReportView<'_>,
    outfmt: OutputFormat,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
) -> Result<()> {
    let ReportView {
        args,
        reports,
        query_ids,
        query_titles,
        subject_ids,
        subject_accessions,
        hit_lists,
        psiblast_hits_found,
        no_alignment_set,
        context,
        ..
    } = *view;
    let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1102-1107
    write_blastp_tabular_hit_lists(
        hit_lists,
        outfmt,
        &mut writer,
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1306-1308
        query_ids,
        query_titles,
        subject_ids,
        subject_accessions,
        view.hit_loads,
        no_alignment_set,
        psiblast_hits_found.as_deref(),
        context,
        probe.as_mut(),
        warnings.as_mut(),
        reports.epilog,
        reports.x_queries_before,
        // psiblast's round (`psi_round`), 1 outside the rounds.
        (args.app == BlastpApp::Psiblast)
            .then(|| reports.psi_round.map_or(1, |io| io.state.number())),
        args.psi_search
            .as_ref()
            .and_then(|search| search.query_other_coding),
    )?;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:4086
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-146
    // NCBI reads the next query batch (c++/src/app/blast/blastp_app.cpp:257-259) after this
    // batch's report, and the reader's messages and error go to `cerr`, which
    // the C++ library ties to `cout`: the rows written so far reach standard
    // output before them. The report is flushed when it is written, as the
    // pairwise report is and as TBLASTN's and TBLASTX's are, so a later
    // batch's stop messages (`search`, `batches.stop`) follow the rows of the
    // batches before it (audit findings A-3 and C-1 of session SFc).
    writer.flush()?;
    Ok(())
}
