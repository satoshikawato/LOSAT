//! What a BLASTP call does with its final hits before any format is written (`run_call`, after
//! the rows of the reports, for the functions named here): the warning lines of each query
//! (`query_warning_lines`), the checks of the reported subjects and each query batch's prefetch of
//! them (`prefetch_reported_subjects`), and the refusals of what LOSAT does not report as NCBI does
//! (`refuse_unreportable`). Each function is the statements of `run_call` at its place, in their
//! order; `run_call` calls them once each, in that order, and its `?` passes their errors on.
//!
//! The NCBI reference for the prefetch before the reports (`BlastFormatter_PreFetchSequenceData`, called at
//! blastp_app.cpp:284-288) is at `prefetch_reported_subjects`. NCBI prefetches each batch's subjects and then writes each query's report, whose warnings come
//! first; LOSAT reads, before the first byte of any format, what the reports will read, and refuses
//! what it cannot write as NCBI does.
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`BlastpReport`, `BlastpReports`, `BlastpRanges`, `blastp_pairwise_view_of`,
//! `blastp_hit_list_hsp_count`) and its sibling modules `database`,
//! `psi_round` and `pssm_hooks`.
//!
//! Before the rows of the reports, `run_call` calls the functions after those three, also once
//! each and in their order: the residues that the reports read of a database's listed subjects
//! (`read_shown_residues`), NCBI's check of the prefetch environment with its stop at the first
//! batch (`refuse_prefetch_environment`), the subjects that the reports show (`shown_of`), what
//! the reports read of the final hit lists' Seq-aligns, psiblast's cut and PSSMs among them
//! (`report_values`), and the labels of the reports (`report_labels`). They read `ReportSetup`,
//! and of `blast_engine` also `ShownSubjects`, `ShownQueries`, `ShownLetters`,
//! `BlastpSubjectSource`, `psiblast_pssms`, `recount_shown_identities`, `write_pairwise_prologs`,
//! `write_request_posts`, `blastp_subject_tag` and `path_label`.

use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::{
    blastp_hit_list_hsp_count, blastp_pairwise_view_of, database, psi_round, pssm_hooks,
    BlastpRanges, BlastpReport, BlastpReports,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::hsp::{collect_hits_and_methods_from_hit_lists, BlastpHitList};
use crate::algorithm::blastp::seqalign::EngineCounts;
use crate::algorithm::common::database::BatchPrefetch;
use crate::algorithm::psiblast::pssm::PssmFiles;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::SubjectReads;

/// What the three functions of this module that run after the rows of the reports read of
/// `run_call`: each field is `run_call`'s value of the same name, by reference where it is not `Copy` (`outputs` for its formats only;
/// `subject_records` and `subject_titles` are those that the reports show, a database's listed
/// subjects' or the subject set's; `format_plans` is the requested formats). The values that
/// `run_call` makes for the prefetch (`batches_prefetch` to `batch_prefetch`) and for the PSSMs
/// (`pssm_made`, `pssm_reads`) come last.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:756-758
pub(super) struct ReportPrep<'a, 'o> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) outputs: &'a ReportOutputs<'o>,
    pub(super) reports: &'a BlastpReports<'a>,
    pub(super) ranges: &'a BlastpRanges<'a>,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) query_valid: &'a [bool],
    pub(super) query_batch_skipped: &'a [bool],
    pub(super) query_batches: &'a [RangedBatch],
    pub(super) database: Option<&'a database::BlastpDatabase>,
    pub(super) listed_titles: &'a Option<database::ListedTitles<'a>>,
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) subject_titles: &'a [Arc<[u8]>],
    pub(super) format_plans: &'a [BlastpReport],
    pub(super) hit_lists: &'a [Option<BlastpHitList>],
    pub(super) pairwise_hits: &'a Option<Arc<Vec<PairwiseHit>>>,
    pub(super) engine_counts: &'a EngineCounts,
    pub(super) psi_found_before: &'a Option<psi_round::FoundBefore>,
    pub(super) psiblast_hits_found: &'a Option<Vec<usize>>,
    pub(super) filtered_empty: bool,
    pub(super) rps_multi_volume: bool,
    pub(super) batches_prefetch: bool,
    pub(super) prefetch_taxids: bool,
    pub(super) prefetch_only: bool,
    pub(super) psiblast_prefetched: &'a Option<Vec<Vec<usize>>>,
    pub(super) batch_ranges: &'a [Range<usize>],
    pub(super) batch_prefetch: &'a BatchPrefetch,
    pub(super) pssm_made: bool,
    pub(super) pssm_reads: &'a database::PssmReads,
}

/// The warning lines of each query of the call, which the first format writes once, before the
/// query's report (`QueryWarnings`): the reader's messages of the query's batch, then, unless the
/// first format is a Seq-annot (which writes no result warnings), the query's own messages (none
/// for rpsblast over several volumes).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:283-291
///
/// Also returned: by query, the queries' own messages that a structured first format does not
/// write, kept for the later formats that do (`Posts::later`, `ReportOutputs::later_warnings`);
/// `None` where they are not kept (`ReportOutputs::keeps_later_query_warnings`).
pub(super) fn query_warning_lines(
    prep: &ReportPrep<'_, '_>,
) -> (Vec<Vec<u8>>, Option<Vec<Vec<u8>>>) {
    let ReportPrep {
        args,
        outputs,
        reports,
        ranges,
        query_records,
        query_valid,
        query_batch_skipped,
        filtered_empty,
        rps_multi_volume,
        ..
    } = *prep;
    // A Seq-annot format writes no result warnings (`ReportFormat::writes_query_warnings`);
    // the warnings go with the first format, or, where it is structured, they are kept for the
    // later formats that write them (`later`).
    let result_warnings = outputs.first_writes_query_warnings();
    let later = outputs.keeps_later_query_warnings();
    let mut later_lines = Vec::new();
    let lines = query_records
        .iter()
        .enumerate()
        .map(|(q_idx, record)| {
            let mut messages = crate::report::query_warnings::query_messages(
                record.seq(),
                query_valid[q_idx],
                reports.x_queries_before + q_idx < reports.identity_warned,
                filtered_empty && !query_batch_skipped[q_idx],
            );
            // rpsblast over several volumes: the merged results keep the errors of the volumes'
            // results and no warning (`s_CombineSearchSets`, rpsblast_local.cpp:165).
            if rps_multi_volume {
                messages.clear();
            }
            let mut lines = reports
                .reader_warnings
                .get(q_idx)
                .cloned()
                .unwrap_or_default();
            if result_warnings || later {
                let text = if let Some(psi) = args.psi_search.as_deref() {
                    // psiblast `-in_pssm`: the label is the PSSM's Seq-id, the record's local id;
                    // the messages add those of `PsiBlastSetupScoreBlock`.
                    pssm_hooks::pssm_query_warning(args.app.ncbi_name(), psi, record, &messages)
                } else {
                    crate::report::query_warnings::query_warning(
                        args.app.ncbi_name(),
                        args.first_query_ordinal
                            + ranges.x_ordinal_base
                            + ranges.query_input.ordinal(q_idx),
                        record,
                        &messages,
                    )
                };
                match result_warnings {
                    true => lines.extend(text),
                    false => later_lines.push(text),
                }
            }
            lines
        })
        .collect();
    (lines, later.then_some(later_lines))
}

/// What query `query`'s pairwise report (outfmt 0 to 4) of its hit list `shown` reads, by rank
/// (`QueryGroup::read_order` of the rows' hits as the writer reads them, `reading`; psiblast's
/// later round's two-part table, `found_before`).
fn pairwise_read_order(
    titles: &database::ListedTitles<'_>,
    shown: &[Option<BlastpHitList>],
    query: usize,
    engine_counts: &EngineCounts,
    reading: crate::report::pairwise::PairwiseReading,
    found_before: Option<&[bool]>,
) -> crate::report::db_subject::PairwiseReads {
    let hits = collect_hits_and_methods_from_hit_lists(shown, query, engine_counts);
    crate::report::x_pairwise_par::QueryGroup::read_order(
        hits.into_iter()
            .map(|(hit, _, seqalign_num_ident, _)| (hit.s_idx, (hit, seqalign_num_ident))),
        |(hit, seqalign_num_ident)| {
            crate::report::hit_sort::HspScores::of_hit(hit, *seqalign_num_ident)
        },
        |previous, s_idx| {
            crate::report::db_subject::rank_joins_as_copy(
                &titles.parsed,
                previous as usize,
                s_idx as usize,
            )
        },
        reading,
        false,
        found_before,
    )
}

/// What the prefetch of the reported subjects gives the reports (`prefetch_reported_subjects`).
pub(super) struct ReportedPrefetch {
    /// The texts of the subjects that a query's report reads one by one where its batch's
    /// prefetch is off, where the first format places them (`ReportFormat::places_subject_reads`).
    pub(super) subject_reads: SubjectReads,
    /// The posts of each batch's prefetch, by the batch's first query, which go after the
    /// reader's messages of the batch (`query_warnings::query_posts`).
    pub(super) posts: Vec<(usize, crate::report::diag_log::Posts)>,
    /// What the XML2 and JSON reports do at each subject's `Hit` (`ListedTitles::hit_loads`), for
    /// the writer.
    pub(super) hit_loads: crate::report::db_subject::HitLoads,
}

/// Check the reported subjects before any format is written and read each query batch's prefetch
/// of them: the database titles and the subjects with parsed Seq-ids that the reports of the
/// batches before a prefetch's stop show and LOSAT does not write as NCBI does are refused, and
/// the stopping batch's prefetch returns NCBI's error. The texts of a batch's prefetch go to the
/// returned posts (`ReportedPrefetch::posts`), with its first query, after the reader's messages
/// of the batch; where a batch's prefetch is off, the texts of the subjects that a query's report
/// reads go to the query's warning lines, or to the returned `ReportedPrefetch::subject_reads`
/// where the first format places them (`ReportFormat::places_subject_reads`), and the texts of
/// the loads of a query's PSSM go to the PSSM files' read messages (`psiblast_pssms`) or to the
/// round's outcome. Also returned: what the XML2 and JSON reports do at each subject's `Hit`
/// (`ListedTitles::hit_loads`), for the writer.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:767-772
pub(super) fn prefetch_reported_subjects(
    prep: &ReportPrep<'_, '_>,
    query_warning_lines: &mut [Vec<u8>],
    psiblast_pssms: &mut Option<PssmFiles>,
) -> Result<ReportedPrefetch> {
    let ReportPrep {
        args,
        outputs,
        reports,
        ranges,
        query_batches,
        database,
        listed_titles,
        format_plans,
        hit_lists,
        engine_counts,
        psiblast_hits_found,
        batches_prefetch,
        prefetch_taxids,
        prefetch_only,
        psiblast_prefetched,
        batch_ranges,
        batch_prefetch,
        pssm_made,
        pssm_reads,
        query_records,
        psi_found_before,
        ..
    } = *prep;
    // Whether the reports write the same where every load of a subject with a positive taxid throws
    // (`db_subject::unchanged_by_load_throws`, of psiblast's formats as of blastp's, as
    // `query_warning_lines` reads them); psiblast's PSSM reads no letters of such a subject
    // (`ListedTitles::pssm_reads`).
    let reports_unchanged_by_load_throws =
        crate::report::db_subject::unchanged_by_load_throws_outfmts(
            outputs.formats.iter().map(|format| format.outfmt),
            crate::blastinput::app::Program::Blastp,
            args.app.request_kind(),
        );
    // How the call's pairwise reports (outfmt 0 to 4) read each query's subjects, as their
    // writers read them (`PairwiseView::reading`; the formats of a run whose loads they make share
    // one view's order, `FormatThrow::and`), and psiblast's later round's previous set (its
    // two-part table, `psiblast_round`).
    let pairwise_reading = format_plans.iter().find_map(|outfmt| match outfmt {
        BlastpReport::Pairwise(view) => Some(view.reading(
            args.sort,
            crate::report::db_subject::Shown {
                descriptions: args.num_descriptions,
                alignments: args.num_alignments,
            },
        )),
        _ => None,
    });
    let found_before = psi_found_before
        .as_ref()
        .map(|found| found.found.as_slice());
    // A database title that LOSAT does not report as NCBI does is an error before any
    // format is written.
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
    // The first batch whose prefetch looks up a Seq-align id whose lookup fails on a file stops
    // the run in its prefetch (`shared::batch_prefetch`, below): NCBI reads none of that batch's
    // reports, nor any of the batches after it, so the checks of the reports are of the batches
    // before it.
    // The XML2 and JSON formats write the taxonomy database's texts at each `Hit`, so of their
    // batches' prefetch only its stop is asked (`shared::batch_prefetch_stop`).
    let xml2_prefetch = listed_titles.is_some()
        && format_plans
            .first()
            .is_some_and(|outfmt| matches!(outfmt, BlastpReport::Xml2(_)));
    // The ranks that the prefetch of the batch of the searched queries `span` reads: psiblast's of
    // its traceback hit lists (`psiblast_prefetched`).
    let prefetched_ranks = |span: std::ops::Range<usize>| -> Vec<usize> {
        match &psiblast_prefetched {
            Some(ranks) => ranks[span].concat(),
            None if args.app.prunes_after_traceback() => hit_lists[span]
                .iter()
                .flat_map(|hit_list| {
                    database::ListedTitles::reported_ranks(std::slice::from_ref(hit_list))
                })
                .collect(),
            None => database::ListedTitles::reported_ranks(&hit_lists[span]).collect(),
        }
    };
    let prefetch_stop = listed_titles
        .as_ref()
        .filter(|_| batches_prefetch || xml2_prefetch)
        .and_then(|titles| titles.prefetch_stop(batch_ranges, batch_prefetch, prefetched_ranks));
    let reported_lists = &hit_lists[..prefetch_stop
        .as_ref()
        .map_or(hit_lists.len(), |span| span.start)];
    // What the XML2 and JSON reports of the batches before a stop (all of them where none stops)
    // do at each subject's `Hit` (`ListedTitles::hit_loads`, the one answer: a `Hit` whose lookup
    // throws where a batch's prefetch is off, and the subjects that no report loads), which the
    // refusals, the stop's record of the loads before it, the names and the taxonomy reads take
    // (review 1 of line DBID pfoff-1, H1 and M1; psiblast writes none of these formats,
    // `psiblast::args::PsiblastArgs::check_outfmt`).
    let hit_loads = listed_titles
        .as_ref()
        .map_or_else(Default::default, |titles| {
            titles.hit_loads(
                reported_lists,
                batch_prefetch,
                crate::report::db_subject::FormatThrow::of_outfmts(
                    outputs.formats.iter().map(|format| format.outfmt),
                    crate::blastinput::app::Program::Blastp,
                    args.app.request_kind(),
                )
                // The call's refusal, over the formats of its other searches too.
                .of_search(
                    outputs.other_outfmts,
                    crate::blastinput::app::Program::Blastp,
                    args.app.request_kind(),
                ),
                |query| {
                    pairwise_reading.map_or_else(Default::default, |reading| {
                        pairwise_read_order(
                            titles,
                            std::slice::from_ref(&reported_lists[query]),
                            query,
                            engine_counts,
                            reading,
                            found_before,
                        )
                    })
                },
            )
        });
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:185-203 (cited at
    // `db_subject::RunLoads`)
    // A report whose load ends the run ends NCBI's application there: no batch after it is
    // prefetched or read one by one, so a stop of a later batch's prefetch is never reached
    // (`HitLoads::reaches`; review 3 of line DBID pfoff-3, N1). The writer ends the run at that
    // report with its error.
    let prefetch_stop = prefetch_stop.filter(|span| hit_loads.reaches(span.start));
    if let Some(titles) = listed_titles.as_ref() {
        titles.check_reported_titles(reported_lists)?;
        // The subjects with parsed Seq-ids that the reports asked for do not write as NCBI
        // (the formats that read a subject's Bioseq need its sequence found by its Seq-align
        // id, which a table of fields that read nothing of the subject does not; outfmt 0, 5 and
        // 13 to 16 join the hits of matching Seq-align ids, which LOSAT does for copies of one
        // record; 13 to 16 and the tables of taxonomy fields write their taxids).
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450
        // psiblast's PSSM (`-out_pssm`, `-out_ascii_pssm`, and each round's PSSM for the next
        // round) reads the letters of each included subject by its Seq-align id, as those formats
        // do (review 1 of line DBID tabfields, H1). `ListedTitles::pssm_reads` reads those, before
        // the cut, where a PSSM is made; a reported subject that no PSSM reads and no format reads is
        // written as NCBI writes it (a `giim` id in rounds with `6 qseqid`, line DBID pssmrows).
        titles.check_parsed_reports(
            reported_lists,
            batch_prefetch,
            format_plans.iter().any(|outfmt| outfmt.reads_subjects()),
            format_plans.iter().any(|outfmt| match outfmt {
                BlastpReport::Pairwise(_) | BlastpReport::Xml | BlastpReport::Xml2(_) => true,
                BlastpReport::CustomTabular(spec) => spec.needs().joins_subjects(),
                _ => false,
            }),
            format_plans.iter().any(|outfmt| match outfmt {
                BlastpReport::Xml2(_) => true,
                BlastpReport::CustomTabular(spec) => spec.needs().taxonomy,
                _ => false,
            }),
            &hit_loads,
            format_plans.iter().fold(
                crate::report::tabular::PrintedIds::default(),
                |words, outfmt| words.or(outfmt.printed_ids()),
            ),
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
    // Before the reports of a batch NCBI fetches every subject that the batch reports, and with
    // `taxdb.bti` on the search path the taxonomy database writes its messages there
    // (`database::prefetch_messages`: the texts of opening it at the run's first lookup, a
    // volume's lookup of a taxid once per run): after the reader's messages of the batch and
    // before the warnings of its first query. They go with the warnings, where the first format
    // is outfmt 0, 5, 6 or 7, or with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or 12, which write
    // the reader's messages only (`BlastpReport::prefetches_subjects`). They are read here, after
    // the subjects' refusals and before any format is written, so a taxonomy database that LOSAT
    // cannot read is refused before any byte; without `taxdb.bti` nothing is read.
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:262-264
    // psiblast's batches, of one query each (`read_psiblast_query_batches`), prefetch the same
    // way, the subjects of the traceback hit lists (`psiblast_prefetched`).
    // Where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off (`s_PreFetchSeqs`, cited at
    // `shared::prefetch_may_be_off`), each query's report fetches the subjects it shows, after
    // the query's own warnings: those of its final hit list, psiblast's after the cut, where the
    // first format's rows read their Bioseqs (`TabularSpec::subject_fetch`; `nident` alone on
    // rows with identities, `ListedTitles::fetched_ranks`).
    // The texts of the subjects that the reports of the queries whose batch's prefetch is off
    // read one by one, by query, where the first format places them
    // (`ReportFormat::places_subject_reads`).
    let mut subject_reads = crate::report::query_warnings::SubjectReads::default();
    // The posts of each batch's prefetch, by the batch's first query.
    let mut prefetch_posts: Vec<(usize, crate::report::diag_log::Posts)> = Vec::new();
    let first_report = outputs
        .formats
        .first()
        .and_then(crate::api::local_blast::FormatOutput::report_format);
    let places_reads =
        first_report.is_some_and(crate::blastinput::app::ReportFormat::places_subject_reads);
    // Outfmt 0 to 4 (`blastp_pairwise_view_of`): NCBI's `PrintOneResultSet` sorts the hits by
    // `-sorthits` and writes the description lines before the alignments in all five
    // (blast_format.cpp:1523-1541), so they read the subjects in one order.
    let pairwise_first = first_report
        .is_some_and(|format| blastp_pairwise_view_of(args.app, format as i32).is_some());
    // The query records that NCBI's scope holds at the prefetch of the batch of the searched
    // queries `span` (`shared::PrefetchQueries`, made only where the batch stops): psiblast's round
    // of a query gives its query's batch (`BlastpReports::prefetch_scope`); the PSSM's query of
    // `-in_pssm` (`PsiSearch::query_scope`); else the input records that NCBI's reader added for
    // the batch (`RangedBatch::input`: the records that `-query_loc` drops included, numbered in
    // the run, `FastaRecord::local_id`), and in outfmt 5, whose scope the formatter keeps
    // (`ResetScopeHistory`, cited at `shared::stop_lookup_unsure`), those of the batches before.
    let xml_scope = crate::blastinput::query_batch::outfmts_have_xml(
        outputs.formats.iter().map(|format| format.outfmt),
    );
    let batch_scope = |span: &std::ops::Range<usize>| {
        if let Some(scope) = reports.prefetch_scope {
            return scope();
        }
        if let Some(psi) = args.psi_search.as_deref() {
            return psi.query_scope();
        }
        let input = query_batches
            .iter()
            .find(|batch| batch.searched == *span)
            .map_or(0..0, |batch| batch.input.clone());
        let input = if xml_scope { 0..input.end } else { input };
        crate::blastinput::query_batch::QueryScope::of_shown_records(&ranges.input_records[input])
    };
    // A call of the library's batches after the first (`api::run_batches`): the calls before
    // searched with a database of their own, whose loads nothing here records.
    let unrecorded_calls = args.first_query_ordinal > 0;
    // Which rows of the first format's report read the subjects' Bioseqs, which NCBI's
    // formatters fetch one by one where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off: a
    // custom tabular report's only for fields that read them (`TabularSpec::subject_fetch`).
    let rows = match format_plans.first() {
        Some(BlastpReport::CustomTabular(spec)) => {
            spec.subject_fetch(crate::blastinput::app::Program::Blastp)
        }
        _ => crate::report::tabular::SubjectFetch::Every,
    };
    // The Seq-aligns of the batch of the searched queries `span`, one for each HSP of its results
    // (psiblast's before the cut, `psiblast_hits_found`), which `PRE_FETCH_SEQS_LIMIT` counts
    // (`shared::prefetch_off`; as `batch_prefetch` counts them).
    let batch_seq_aligns = |span: std::ops::Range<usize>| -> usize {
        match psiblast_hits_found.as_ref() {
            Some(found) => found[span].iter().sum(),
            None => blastp_hit_list_hsp_count(&hit_lists[span]),
        }
    };
    // Outfmt 5 first: whether NCBI's XML stream cuts a query's report before its hits
    // (`shared::XmlCuts`), the one answer for the run's record of the batches' loads (`loaded`:
    // a cut report loads nothing) and for a later batch's prefetch stop, which asks for the texts
    // that the cut reports left (`shared::XmlReportsBefore`, from the reads of the queries'
    // reports modelled here).
    let xml_first = first_report == Some(crate::blastinput::app::ReportFormat::Xml);
    let xml_cut = |query: usize| {
        query_records.get(query).is_some_and(|record| {
            crate::report::blast_xml::names::iteration_cut_before_hits(
                record.local_id.as_bytes(),
                record.warning_id(),
                record.believed(),
                &record.title,
                true,
            )
        })
    };
    let xml_cuts = crate::algorithm::common::database::XmlCuts::of(xml_first, &xml_cut);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:764-766 (cited at
    // `shared::batch_loads`)
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:811-813 (cited at
    // `shared::BatchLoads`: the tabular report reads the Seq-aligns after the cut to the hitlist
    // size)
    // The subjects, by rank, that NCBI's BLAST database data loader loaded for the call's batches
    // before the searched query `start` (`shared::LoadedBefore::subjects`, review 3 of line DBID
    // loadord, H2): a batch whose prefetch is on loads the subjects of its Seq-aligns (psiblast's
    // of its traceback hit lists, `prefetched_ranks`, beyond `-max_target_seqs`); one whose
    // prefetch is off, those that the rows of the reports NCBI writes read, of the hit lists
    // they show (after psiblast's cut, `ListedTitles::fetched_ranks`; review 6, H1; none of an
    // outfmt 5 report that NCBI's XML stream cut, `xml_cuts`: review 6, H2), and the subjects
    // that psiblast's PSSM includes (`PssmReads::one_by_one`, where a PSSM is made). Made only
    // where the batch at `start` stops, and for the run's record of psiblast's rounds where a
    // prefetch of the run can stop (`ListedTitles::note_loaded`).
    let hit_loads_ref = &hit_loads;
    let loaded = |start: usize| {
        crate::algorithm::common::database::LoadedBefore {
            hits: &hit_lists[..start],
            spans: batch_ranges,
            rows,
            cuts: xml_cuts,
        }
        .subjects(
            |hits, span| hits.get(span),
            move |span, _| batch_seq_aligns(span),
            move |span, hits, loads| match loads {
                crate::algorithm::common::database::BatchLoads::Prefetched => {
                    prefetched_ranks(span)
                }
                crate::algorithm::common::database::BatchLoads::Read(rows) => {
                    // Not a subject whose load throws, nor one that no report reaches after
                    // such a subject (`HitLoads::loaded`).
                    let mut ranks: Vec<usize> = database::ListedTitles::fetched_ranks(hits, rows)
                        .filter(|&rank| hit_loads_ref.loaded(rank))
                        .collect();
                    if pssm_made {
                        for query in span {
                            ranks.extend_from_slice(pssm_reads.one_by_one(query));
                        }
                    }
                    ranks
                }
            },
        )
    };
    if let (Some(titles), Some(database)) = (listed_titles.as_ref(), database) {
        // Without `taxdb.bti` the prefetch of a table that reads the subjects writes nothing
        // before it stops: the stopping batch's prefetch alone (`ListedTitles::prefetch_stop`).
        // The XML2 and JSON formats: the stop alone (`shared::held_batch_prefetch_stop`, after the
        // lookups of the reports of the batches before), whose other subjects' loads write their
        // texts at each `Hit`.
        if let Some(span) = prefetch_stop.as_ref().filter(|_| xml2_prefetch) {
            let scope = || batch_scope(span);
            // The loads of the batches before (psiblast's traceback beyond the cut, which no
            // report reads), beside their reports' lookups, which the stop records as every stop
            // records the call's loads before it (without the record's gate, f8 decision b of
            // line DBID loadord; review 5 L3).
            titles.batch_prefetch_stop(
                database,
                &hit_lists[..span.start],
                loaded(span.start),
                &hit_loads,
                prefetched_ranks(span.clone()),
                &crate::algorithm::common::database::PrefetchQueries::new(span.start, &scope)
                    .after_unrecorded_calls(unrecorded_calls),
            )?;
        } else if let Some(span) = prefetch_stop
            .as_ref()
            .filter(|_| !(prefetch_taxids || psiblast_prefetched.is_some() || prefetch_only))
        {
            let scope = || batch_scope(span);
            // No report's reads are modelled before it (`prefetch_taxids` is off).
            titles.batch_prefetch(
                database,
                loaded(span.start),
                database::ListedTitles::reported_ranks(&hit_lists[span.clone()]),
                &crate::algorithm::common::database::PrefetchQueries::new(span.start, &scope)
                    .after_unrecorded_calls(unrecorded_calls)
                    .where_reports_unchanged_by_load_throws(reports_unchanged_by_load_throws),
            )?;
        }
        if prefetch_taxids || psiblast_prefetched.is_some() || prefetch_only {
            // A custom tabular report with a plural name looks every taxid of its rows up, taxid 0
            // too, which the prefetch does not: the rows of a batch with hits may open the
            // taxonomy database before a later batch's prefetch (`TaxdbRun::open_for_rows`).
            let plural_names = matches!(
                format_plans.first(),
                Some(BlastpReport::CustomTabular(spec)) if spec.looks_up_every_taxid()
            );
            // A row whose subject's load throws looks nothing up (`HitLoads::loaded`).
            let open_for_rows = |lists: &[Option<BlastpHitList>]| -> Result<()> {
                if plural_names
                    && database::ListedTitles::reported_ranks(lists)
                        .any(|rank| hit_loads.loaded(rank))
                {
                    database.open_for_tabular_rows(args.app.name())?;
                }
                Ok(())
            };
            for batch in query_batches {
                // No batch after the run's end is prefetched, nor read one by one.
                if !hit_loads.reaches(batch.searched.start) {
                    break;
                }
                let lists = &hit_lists[batch.searched.clone()];
                let seq_aligns = batch_seq_aligns(batch.searched.clone());
                if crate::algorithm::common::database::prefetch_off(seq_aligns) {
                    if !prefetch_taxids {
                        continue;
                    }
                    for query in batch.searched.clone() {
                        if !hit_loads.reaches(query) {
                            break;
                        }
                        if rows != crate::report::tabular::SubjectFetch::Never {
                            let shown = std::slice::from_ref(&hit_lists[query]);
                            // Rows that read a subject only for its identities (`nident` alone)
                            // read none of a query whose HSPs have none, which looks no taxid up.
                            // A subject whose load throws loads nothing, nor does one that no
                            // report reaches after it (`HitLoads::loaded`).
                            let loaded = |rank: &usize| hit_loads.loaded(*rank);
                            let mut fetched = database::ListedTitles::fetched_ranks(shown, rows)
                                .filter(loaded)
                                .peekable();
                            if rows == crate::report::tabular::SubjectFetch::Every
                                || fetched.peek().is_some()
                            {
                                // The first format writes the texts where it reads each
                                // subject (`places_subject_reads`), else before the report.
                                if places_reads {
                                    // The pairwise report reads the subjects that it shows as
                                    // its writer reads them (`QueryGroup::read_order`,
                                    // `PairwiseReads::subjects`), the other formats in their
                                    // rows' order.
                                    let ranks: Vec<usize> =
                                        match pairwise_reading.filter(|_| pairwise_first) {
                                            Some(reading) => pairwise_read_order(
                                                titles,
                                                shown,
                                                query,
                                                engine_counts,
                                                reading,
                                                found_before,
                                            )
                                            .subjects()
                                            .into_iter()
                                            .filter(loaded)
                                            .collect(),
                                            None => fetched.collect(),
                                        };
                                    subject_reads.set(
                                        query,
                                        titles.subject_read_texts(database, ranks, plural_names)?,
                                    );
                                } else {
                                    let messages = titles.prefetch_messages(database, fetched)?;
                                    query_warning_lines[query].extend(messages);
                                }
                            }
                            open_for_rows(shown)?;
                        }
                        // The query's PSSM reads after its report (psiblast_app.cpp:283-308):
                        // the blobs it loads first look their taxids up then, and NCBI writes
                        // their texts among the PSSM's warnings (`ListedTitles::pssm_load_texts`).
                        // A report that ends the run makes no PSSM (`HitLoads::reaches`).
                        if pssm_made && hit_loads.reaches(query + 1) {
                            let texts = titles.pssm_load_texts(
                                database,
                                pssm_reads,
                                query,
                                hit_lists[query].as_ref(),
                                args.inclusion_ethresh,
                            )?;
                            if texts.iter().any(|(_, text)| !text.is_empty()) {
                                if let Some(files) = psiblast_pssms.as_mut() {
                                    files.read_messages(query).add_load_texts(texts);
                                } else if let Some(io) = reports.psi_round {
                                    let mut outcome =
                                        io.outcome.lock().expect("the round's outcome");
                                    if let Some(outcome) = outcome.as_mut() {
                                        outcome.next_warnings.add_load_texts(texts);
                                    }
                                }
                            }
                        }
                    }
                    continue;
                }
                // The taxonomy texts, and NCBI's warning where the prefetch fails at a subject
                // whose lookup fails (`ListedTitles::batch_prefetch`), also without `taxdb.bti`.
                let prefetched = |query: usize| {
                    prefetch_posts
                        .iter()
                        .find(|(first, _)| *first == query)
                        .map(|(_, posts)| posts)
                };
                let xml_before = xml_cuts.reports_before(subject_reads.queries(), &prefetched);
                let scope = || batch_scope(&batch.searched);
                let queries = crate::algorithm::common::database::PrefetchQueries::new(
                    batch.searched.start,
                    &scope,
                )
                .after_unrecorded_calls(unrecorded_calls)
                .after_xml_reports(xml_before.as_ref())
                .where_reports_unchanged_by_load_throws(reports_unchanged_by_load_throws);
                let messages = match psiblast_prefetched.as_ref() {
                    Some(ranks) => {
                        let prefetched = titles.prefetched_ranks(
                            database,
                            &ranks[batch.searched.clone()].concat(),
                            prefetch_taxids,
                            !prefetch_only,
                        )?;
                        titles.batch_prefetch(
                            database,
                            loaded(batch.searched.start),
                            prefetched,
                            &queries,
                        )?
                    }
                    None => {
                        if !prefetch_taxids {
                            titles.check_prefetched(
                                database,
                                database::ListedTitles::reported_ranks(lists),
                            )?;
                        }
                        titles.batch_prefetch(
                            database,
                            loaded(batch.searched.start),
                            database::ListedTitles::reported_ranks(lists),
                            &queries,
                        )?
                    }
                };
                open_for_rows(lists)?;
                // After the reader's messages of the batch, in its first query
                // (`query_warnings::query_posts`).
                if !messages.is_empty() {
                    prefetch_posts.push((batch.searched.start, messages));
                }
            }
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:246-264
        // Every round prefetches into the one scope of the query's batch, which the formatter
        // resets after the batch (`ResetScopeHistory`, psiblast_app.cpp:513), and into the one
        // BLAST database data loader of the run, which keeps every Seq-id it registered
        // (`RegisterIds`, cited at `shared::stop_lookup_unsure`).
        // psiblast's round of a query is a call of its own: the run's record of what its prefetch
        // (the traceback hit lists, `psiblast_prefetched`), or where it is off its reports and its
        // PSSM's reads, loaded (`loaded`), for the stop of a later round's or query's prefetch
        // (`shared::TaxdbRun::note_loaded`, on the database that the rounds share). The call
        // reaches here where its prefetch did not stop.
        if reports.psi_round.is_some() {
            titles.note_loaded(database, || loaded(hit_lists.len()));
        }
    }
    Ok(ReportedPrefetch {
        subject_reads,
        posts: prefetch_posts,
        hit_loads,
    })
}

/// Refuse, before any format of the call is written, what LOSAT does not report as NCBI does
/// (with outfmt 0 and `-parse_deflines`, a query that is also a subject and has masks: outfmt 0's
/// refusal alone, which `run_call` makes after this one and puts first, `refuse_shared_masks_in`):
/// with outfmt 1 to 4 of a `-subject` search, a query whose own record is one of the subjects and
/// whose query-anchored section NCBI aborts on (the error is that abort's); and with outfmt 1 to
/// 4 of a database search, or outfmt 0 to 4 with the queries' ids, a reported sequence that shows
/// a query's id of the batch that reports it. `rows_of` makes the rows of a window of queries,
/// where the call has not made every row before (`pairwise_hits`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1412
pub(super) fn refuse_unreportable(
    prep: &ReportPrep<'_, '_>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    let ReportPrep {
        args,
        reports,
        ranges,
        query_records,
        query_batches,
        database,
        listed_titles,
        subject_records,
        subject_titles,
        format_plans,
        hit_lists,
        pairwise_hits,
        psi_found_before,
        engine_counts,
        ..
    } = *prep;
    // Outfmt 1 to 4 of a `-subject` search: a query whose own record is one of the subjects
    // (`-parse_deflines`) and whose query-anchored section NCBI BLAST+ aborts on is refused
    // before any format of the call is written (`pairwise::blastp_section_aborts`; line RG T3a,
    // PD-LOSAT-NCBI-DEFECTS rule 3): in a block search (`x_query_block_plan`), the call of the
    // query's block, after the bytes of the blocks before, as outfmt 0's refusal of the masks of
    // such a record (`run_call`, `believed::refuse_shared_masks_in`; NCBI's own bytes before its
    // abort depend on its buffers, so no prefix is NCBI's). The formats 1 to 4 merge one section (their views differ in the letters they
    // show).
    let anchored_view = format_plans.iter().find_map(|outfmt| match outfmt {
        BlastpReport::Pairwise(view)
            if *view != crate::report::pairwise::PairwiseView::Pairwise =>
        {
            Some(*view)
        }
        _ => None,
    });
    // Outfmt 1 to 4 of an `-ungapped` search that shows an HSP of length 0: the merge of NCBI's
    // query-anchored section can crash on it (`empty_alignments::refuse_merge_frees`), unless
    // the search ends with an exception first (`-subject`'s mapper).
    if let Some(view) = anchored_view {
        let parsed = listed_titles
            .as_ref()
            .map_or(&[][..], |titles| titles.parsed.as_slice());
        super::empty_alignments::refuse_merge_frees(
            view,
            args.sort,
            args.num_alignments,
            hit_lists,
            pairwise_hits.as_deref().map(Vec::as_slice),
            rows_of,
            subject_titles,
            &|previous, s_idx| {
                crate::report::db_subject::rank_joins_as_copy(
                    parsed,
                    previous as usize,
                    s_idx as usize,
                )
            },
            super::empty_alignments::subject_start_exception(
                hit_lists,
                ranges.subject_placements,
                database.is_none(),
            )
            .is_some(),
        )?;
    }
    if let (Some(view), true, None) = (anchored_view, args.parse_deflines, database) {
        for (q_idx, record) in query_records.iter().enumerate() {
            if record.believed().is_none_or(|ids| ids.shared.is_none()) {
                continue;
            }
            let rows: std::borrow::Cow<'_, [PairwiseHit]> = match pairwise_hits.as_deref() {
                Some(hits) => std::borrow::Cow::Borrowed(hits.as_slice()),
                None => std::borrow::Cow::Owned(rows_of(q_idx..q_idx + 1)?),
            };
            let aborts = crate::report::pairwise::blastp_section_aborts(
                view,
                args.sort,
                args.num_alignments,
                q_idx,
                &rows,
                subject_titles,
                |s_idx| {
                    subject_records.get(s_idx as usize).is_some_and(|subject| {
                        crate::report::query_anchored::subject_is_query_record(record, subject)
                    })
                },
            )?;
            if aborts {
                return Err(crate::report::query_anchored::aborted_section(
                    record,
                    args.app.name(),
                ));
            }
        }
    }
    // The subjects whose description lines the reports of query `q_idx` show, by rank (outfmt 0
    // to 4; `anchored`: of outfmt 1 to 4 only): the union over the call's pairwise formats of
    // `pairwise::blastp_described_subjects` of the rows that the reports read (sorted by the
    // format's `-sorthits`, none in psiblast and rpsblast; psiblast's round after the first shows
    // the HSPs of its two-part table, `psiblast_round`). Each such line looks its title up again
    // by the subject's display id (`ParsedSubject::display_id`), which the refusals below ask
    // for, only for a query that reports a subject whose display id its batch has.
    let pairwise_views: Vec<crate::report::pairwise::PairwiseView> = format_plans
        .iter()
        .filter_map(|outfmt| match outfmt {
            BlastpReport::Pairwise(view) => Some(*view),
            _ => None,
        })
        .collect();
    let described = |q_idx: usize, anchored: bool| -> Result<Vec<u32>> {
        let rows: std::borrow::Cow<'_, [PairwiseHit]> = match pairwise_hits.as_deref() {
            Some(hits) => std::borrow::Cow::Borrowed(hits.as_slice()),
            None => std::borrow::Cow::Owned(rows_of(q_idx..q_idx + 1)?),
        };
        let parsed = listed_titles
            .as_ref()
            .map_or(&[][..], |titles| titles.parsed.as_slice());
        let mut described: Vec<u32> = Vec::new();
        for &view in &pairwise_views {
            if anchored && view == crate::report::pairwise::PairwiseView::Pairwise {
                continue;
            }
            for s_idx in crate::report::pairwise::blastp_described_subjects(
                view.report_sort(args.sort),
                psi_found_before
                    .as_ref()
                    .map(|found| found.found.as_slice()),
                args.num_descriptions,
                q_idx,
                &rows,
                |previous, s_idx| {
                    crate::report::db_subject::rank_joins_as_copy(
                        parsed,
                        previous as usize,
                        s_idx as usize,
                    )
                },
            ) {
                if !described.contains(&s_idx) {
                    described.push(s_idx);
                }
            }
        }
        Ok(described)
    };
    // Outfmt 1 to 4 of a database search whose queries are read without `-parse_deflines`: a
    // reported database sequence whose Seq-align id is the generated id (`lcl|Query_N`) of a query
    // of the batch that reports it is that query's own Bioseq for NCBI's scope, which shows the
    // query's title and letters in its place and can abort on its merge (review 2 of line RG
    // T3a, D2), and one whose description line shows such an id as its display id gets the
    // query's title there (review 3 of T3a, D1): refused here, with the abort above
    // (`ListedTitles::refuse_generated_query_ids`). Under `-parse_deflines` every query has a
    // believed id (generated for a record without one) and the refusal of a reported sequence with
    // a query's id is the scopes' of the batches, for every format by its Seq-align id
    // (`ListedTitles::refuse_ids_of_queries` above) and for outfmt 0 to 4 by its display id
    // (below). A `-subject` search cannot meet it: without `-parse_deflines` the subjects are
    // `lcl|Subject_N`, and with it a subject with a query's id is the reader's refusal
    // (`id_scope`) unless it is the query's own record (the abort above). NCBI's display of such a
    // sequence in outfmt 0, 6 and 7 (the query's title and letters) is line RG item X7's, the
    // owner of this behaviour for every format. The refusal runs in the search of all queries,
    // before any byte of the call's reports (in a psiblast round after the first, after the
    // earlier rounds' output): such a sequence has parsed Seq-ids, and the block
    // search never takes a database with them (`database::x_titles_allow_blocks`). Were a block
    // search to reach it, it would refuse in the call of the query's block as the abort above,
    // after the reports of the blocks before (NCBI's bytes for their queries).
    // rpsblast's BLAST XML, XML2 and JSON (outfmt 5, 13 to 16): the same collision of an entry's
    // generated id with a query of its batch, where NCBI's documents take the query's title for
    // the entry (the archive, outfmt 11, and the Seq-annots, 8, 9 and 12, are NCBI's bytes with
    // it): refused here before any byte, for the reported entries only
    // (`ListedTitles::refuse_generated_ids`; PS15a-1 and PS15a-2, followup
    // `rps_lcl_query_title.md`). Outfmt 0, 6, 7 and 10 are line RG item X7's.
    let rps_xml = args.rps_search.is_some()
        && format_plans
            .iter()
            .any(|outfmt| matches!(outfmt, BlastpReport::Xml | BlastpReport::Xml2(_)));
    // rpsblast's organism report (outfmt 18) and comma-separated values with a header line
    // (outfmt 20) where NCBI's rows show what its formatter finds by such an id, the query
    // (`TabularSpec::refuses_generated_id`, PS15a-2): refused as the documents above; a table of
    // other fields, and the organism report of a profile without a taxid, are NCBI's bytes.
    // The reports that can refuse (an HSP with identities, a profile with a taxid).
    let rps_tables: Vec<&crate::report::tabular::TabularSpec> = if args.rps_search.is_some() {
        format_plans
            .iter()
            .filter_map(|outfmt| match outfmt {
                BlastpReport::CustomTabular(spec)
                    if spec.refuses_generated_id(
                        crate::blastinput::app::Program::Blastp,
                        1,
                        true,
                    ) =>
                {
                    Some(*spec)
                }
                _ => None,
            })
            .collect()
    } else {
        Vec::new()
    };
    if let (true, false, Some(titles)) = (
        anchored_view.is_some() || rps_xml || !rps_tables.is_empty(),
        args.parse_deflines,
        listed_titles.as_ref(),
    ) {
        let scopes: Vec<(usize, crate::blastinput::query_batch::QueryScope)> = query_batches
            .iter()
            .map(|batch| {
                (
                    batch.searched.end,
                    crate::blastinput::query_batch::QueryScope::of_generated_ids(
                        ranges.input_records[batch.input.clone()]
                            .iter()
                            .map(|record| record.local_id.as_str()),
                    ),
                )
            })
            .collect();
        let none = crate::blastinput::query_batch::QueryScope::default();
        let scope_of = |query: usize| {
            scopes
                .get(scopes.partition_point(|&(end, _)| end <= query))
                .map_or(&none, |(_, scope)| scope)
        };
        if rps_xml {
            titles.refuse_generated_ids(
                hit_lists,
                &scope_of,
                |_, _, _| true,
                crate::algorithm::common::database::XML_EFFECT,
            )?;
        }
        if !rps_tables.is_empty() {
            titles.refuse_generated_ids(
                hit_lists,
                &scope_of,
                |query, list_index, subject| {
                    let taxid = subject.first_source_taxid().is_some();
                    let Some(Some(hit_list)) = hit_lists.get(query) else {
                        return false;
                    };
                    hit_list.hsplist_array[list_index]
                        .hsps
                        .iter()
                        .enumerate()
                        .any(|(index_in_list, hsp)| {
                            // The Seq-align's `num_ident` (`EngineCounts`).
                            let num_ident = engine_counts.num_ident(
                                (query as u32, list_index as u32, index_in_list as u32),
                                hsp.num_ident,
                            );
                            rps_tables.iter().any(|spec| {
                                spec.refuses_generated_id(
                                    crate::blastinput::app::Program::Blastp,
                                    num_ident,
                                    taxid,
                                )
                            })
                        })
                },
                crate::algorithm::common::database::TABLE_EFFECT,
            )?;
        }
        if anchored_view.is_some() {
            titles
                .refuse_generated_query_ids(hit_lists, scope_of, |q_idx| described(q_idx, true))?;
        }
    }
    // Outfmt 0 to 4 of a database search with queries whose ids the reports show
    // (`-parse_deflines`, a PSSM query with its ids): a sequence whose description line shows a
    // query's id of the batch that reports it as its display id (`gi|779|lcl|abc` and a query
    // `lcl|abc`), where NCBI shows the query's title (review 3 of line RG T3a, D1): refused here
    // with the refusals above, before any format is written
    // (`ListedTitles::refuse_display_ids_of_queries`; by its Seq-align id in every format,
    // `ListedTitles::refuse_ids_of_queries`). It runs in the search of all queries, as the
    // refusal above: before any byte of the call, in a psiblast round after the first after the
    // earlier rounds' output.
    if let (false, Some(titles), Some(scopes)) = (
        pairwise_views.is_empty(),
        listed_titles.as_ref(),
        reports.query_scopes,
    ) {
        titles.refuse_display_ids_of_queries(hit_lists, scopes, |q_idx| described(q_idx, false))?;
    }
    Ok(())
}

// The steps from the final hit lists to the labels of the reports (`ReportSetup`), before the
// rows of the reports. `hit_list_to_seqalign_order` is `blast_engine`'s import (from `hsp`),
// which only this module reads.
use super::{
    blastp_subject_tag, hit_list_to_seqalign_order, path_label, psiblast_pssms,
    recount_shown_identities, write_pairwise_prologs, write_request_posts, BlastpSubjectSource,
    ShownLetters, ShownQueries, ShownSubjects,
};
use crate::algorithm::blastp::args::BlastpApp;
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::report::{OutputFormat, ReportContext};
use crate::stats::KarlinParams;

/// What the functions from `read_shown_residues` to `report_labels` read of `run_call`: each
/// field is `run_call`'s value of the same name, by reference where it is not `Copy`
/// (`requested_formats` is the requested formats, `run_call`'s `format_plans` after them;
/// `total_db_len` is the letters of the searched subjects, NCBI's `Int8`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:283-285
/// NCBI makes its result set of the Seq-aligns of the final hit lists and the queries' ids, which
/// its formatter reads; LOSAT reads what its reports need of the final hit lists with these values.
pub(super) struct ReportSetup<'a> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) reports: &'a BlastpReports<'a>,
    pub(super) ranges: &'a BlastpRanges<'a>,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) query_batches: &'a [RangedBatch],
    pub(super) database: Option<&'a database::BlastpDatabase>,
    pub(super) requested_formats: &'a [BlastpReport],
    pub(super) contexts: &'a [QueryContext],
    pub(super) gapped_params: &'a KarlinParams,
    pub(super) subject_source: BlastpSubjectSource<'a>,
    pub(super) subject_label: &'a str,
    pub(super) total_db_len: i64,
    /// rpsblast: the HSPs of each query's profiles that left its merged list past the
    /// formatter's cut (`rps_hooks::RpsSearched::unlisted`), which outfmt 7's count and the
    /// prefetch's limit still count; empty otherwise.
    pub(super) unlisted: &'a [usize],
}

/// The residues that the reports read of a database's listed subjects, read from the database
/// again after the redo (`ListedSubjects::x_read_shown`); nothing for a subject set.
/// `x_listed_lean` is LOSAT_X_BPLISTEDLEAN's switch (`x_bp_listed_lean`), which `run_call` reads
/// once for this read and the read of the titles after it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1110-1113
/// NCBI loads the subjects of the Seq-aligns into the formatter's scope, from which its reports
/// read them. It loads them there when it prefetches (blast_app_util.cpp:764-772), after its read
/// of `PRE_FETCH_SEQS_LIMIT`, and otherwise when a report reads each subject; LOSAT reads them
/// here for every format.
pub(super) fn read_shown_residues(
    setup: &ReportSetup<'_>,
    outputs: &ReportOutputs<'_>,
    listed_subjects: &mut Option<database::ListedSubjects>,
    hit_lists: &[Option<BlastpHitList>],
    x_listed_lean: bool,
) -> Result<()> {
    let ReportSetup {
        args,
        reports,
        query_records,
        database,
        requested_formats,
        ..
    } = *setup;
    // The residues that the reports read, read from the
    // database again after the redo (`ListedSubjects::x_read_shown`): every listed subject's for
    // the records; where every report is tabular rows of the default fields or Seq-aligns, which
    // read the letters only to recount the HSPs whose query or subject has an O, those of the
    // subjects with an O and of the HSPs whose query has one; every final subject's otherwise.
    if let (Some(listed), Some(database)) = (listed_subjects.as_mut(), database) {
        let residues = if !x_listed_lean {
            database::XShownResidues::Listed
        } else if outputs.hits.is_none()
            && args.psi_pssm.is_none()
            && args.rps_search.is_none()
            && reports.psi_round.is_none()
            && requested_formats.iter().all(|outfmt| {
                matches!(
                    outfmt,
                    BlastpReport::Seqalign(_)
                        | BlastpReport::Archive
                        | BlastpReport::Text(
                            OutputFormat::Tabular | OutputFormat::TabularWithComments
                        )
                )
            })
        {
            database::XShownResidues::WithO
        } else {
            database::XShownResidues::Final
        };
        // The query test of `recount_shown_identities`.
        let query_recounted: Vec<bool> = if residues == database::XShownResidues::WithO {
            let queries = ShownQueries::of(query_records, args);
            (0..query_records.len())
                .map(|q_idx| queries.letters(q_idx).is_some_and(ShownLetters::recounted))
                .collect()
        } else {
            Vec::new()
        };
        listed.x_read_shown(database, hit_lists, residues, &query_recounted)?;
    }
    Ok(())
}

/// NCBI's read of `PRE_FETCH_SEQS_LIMIT` before the first batch's reports
/// (`app::check_prefetch_environment`): where NCBI cannot convert its value, what NCBI writes
/// before it stops is written (the first batch's search strategy request to the strategy file,
/// `CallRequests::end_at`; the outfmt 0 prolog, the batch's reader messages and the posts of that
/// request to `outputs`), and NCBI's error is returned.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:732-734
pub(super) fn refuse_prefetch_environment(
    setup: &ReportSetup<'_>,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let ReportSetup {
        args,
        reports,
        subject_source,
        subject_label,
        ..
    } = *setup;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:280-285
    // The prefetch before a batch's reports reads `PRE_FETCH_SEQS_LIMIT`, after the batch's search
    // and the read of its definition lines (in `run_call`): a value that NCBI cannot convert stops the
    // first batch here, after the outfmt 0 prolog and the messages of its records, before any of
    // its reports (`app::check_prefetch_environment`; the first call of a run that writes reports
    // meets it: not a chunk's search, which returned earlier in `run_call`). psiblast's iterations fetch ahead
    // before they choose the archive, rpsblast's and blastp's write the archive without it.
    let archive_alone = outputs.only_archive() && args.app != BlastpApp::Psiblast;
    if let Err(error) = crate::blastinput::app::check_prefetch_environment(archive_alone) {
        // NCBI wrote the first batch's search strategy request before its search (`CallRequests`).
        if let Some(call) = reports.call_requests {
            call.end_at(call.first_batch_end())?;
        }
        if reports.x_prolog {
            write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
        }
        if let Some(first_batch_messages) = reports.reader_warnings.first() {
            outputs.diagnostics.write_all(first_batch_messages)?;
        }
        // The posts of the bytes replaced in that request follow the batch's reader messages
        // (`SaveSearchStrategy` after `GetNextSeqBatch`, blastp_app.cpp:259-262).
        write_request_posts(args, &mut *outputs.diagnostics)?;
        return Err(error);
    }
    Ok(())
}

/// The records, the titles and the letters of the subjects that the reports show
/// (`ShownSubjects`): a database's listed subjects' (`listed_titles`; with LOSAT_X_BPLISTEDLEAN
/// their titles by rank and their stored residues, `subjects`), or else the subject set's
/// (`set_records`, `set_titles`: `run_call`'s `subject_records` and `subject_titles` before the
/// returned ones shadow them).
///
/// No NCBI counterpart: how LOSAT's reports reach the subjects' letters and titles; NCBI's
/// formatter reads them from its scope.
#[allow(clippy::type_complexity)]
pub(super) fn shown_of<'a>(
    listed_titles: &'a Option<database::ListedTitles<'_>>,
    subjects: &'a [EncodedProtein],
    set_records: &'a [FastaRecord],
    set_titles: &'a [Arc<[u8]>],
) -> (
    &'a [FastaRecord],
    std::borrow::Cow<'a, [Arc<[u8]>]>,
    ShownSubjects<'a>,
) {
    match listed_titles.as_ref().map(|titles| &titles.shown) {
        // A database's subjects have no records with the switch: every reader of a listed
        // subject takes `shown_subjects`, the titles by rank or the stored residues.
        Some(database::ListedShown::Lean(titles)) => (
            &[],
            std::borrow::Cow::Borrowed(&titles[..]),
            ShownSubjects::Listed { subjects, titles },
        ),
        Some(database::ListedShown::Records(records)) => (
            records,
            std::borrow::Cow::Owned(
                records
                    .iter()
                    .map(|record| Arc::from(record.title()))
                    .collect(),
            ),
            ShownSubjects::Records(records),
        ),
        None => (
            set_records,
            std::borrow::Cow::Borrowed(set_titles),
            ShownSubjects::Records(set_records),
        ),
    }
}

/// What `report_values` makes for the reports, each `run_call`'s value of the same name: the
/// prefetch of each batch (`batches_prefetch` to `batch_prefetch`), the PSSMs' reads, files and
/// round outcome (`pssm_reads` to `pssm_made`), psiblast's numbers of HSPs before its cut by query
/// (`psiblast_hits_found`), and the engine's identity counts of the HSPs that the recount changed
/// (`engine_counts`).
///
/// No NCBI counterpart: LOSAT's values between its Seq-aligns and its reports; NCBI's result set
/// and formatter hold theirs (`report_values`).
pub(super) struct ReportValues {
    pub(super) batches_prefetch: bool,
    pub(super) prefetch_taxids: bool,
    pub(super) prefetch_only: bool,
    pub(super) psiblast_prefetched: Option<Vec<Vec<usize>>>,
    pub(super) batch_ranges: Vec<Range<usize>>,
    pub(super) batch_prefetch: BatchPrefetch,
    pub(super) pssm_reads: database::PssmReads,
    pub(super) psiblast_pssms: Option<PssmFiles>,
    pub(super) psi_found_before: Option<psi_round::FoundBefore>,
    pub(super) pssm_made: bool,
    pub(super) psiblast_hits_found: Option<Vec<usize>>,
    pub(super) engine_counts: EngineCounts,
}

/// What the reports read of the Seq-aligns of the final hit lists, made from them before any
/// report: the lists in the Seq-aligns' order, whether and how each batch prefetches the
/// subjects, the PSSMs' reads, `-out_pssm`'s PSSM files (staged, `PssmFiles`) and a psiblast
/// round's outcome, psiblast's cut of the lists to the formatter's (after the prefetch's ranks and
/// the PSSMs, which read the lists before it), and the identities that the reports show
/// (`recount_shown_identities`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1549-1551
pub(super) fn report_values(
    setup: &ReportSetup<'_>,
    listed_subjects: &Option<database::ListedSubjects>,
    listed_titles: &Option<database::ListedTitles<'_>>,
    shown_subjects: ShownSubjects<'_>,
    hit_lists: &mut [Option<BlastpHitList>],
) -> Result<ReportValues> {
    let ReportSetup {
        args,
        reports,
        ranges,
        query_records,
        query_batches,
        database,
        requested_formats,
        contexts,
        gapped_params,
        ..
    } = *setup;
    // The reports are made from the Seq-aligns of the final hit lists, in the order and with
    // the HSPs that BlastHitList2SeqAlign_OMF gives them (`hit_list_to_seqalign_order`).
    for hit_list in hit_lists.iter_mut().flatten() {
        hit_list_to_seqalign_order(hit_list);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:804-813
    // psiblast's Seq-aligns are those of its larger traceback hit lists (`traceback_hitlist_size`):
    // the `# N hits found` line of -outfmt 7 counts them all, and the formatter keeps the first
    // `max_target_seqs` subjects that have one (the pairwise report shows at most as many, its
    // numbers of descriptions and alignments).
    // NCBI's cut counts the Seq-aligns of subjects whose Seq-align ids match as those of one
    // subject (`PruneSeqalign`, cited at `joined_cut`): a copy of one record that follows
    // another is kept with it; where it keeps another subject for others, the run is LOSAT's
    // refusal before any format is written.
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1066-1133
    // The archive holds the results as the search made them: no formatter cuts rpsblast's merged
    // volumes or psiblast's lists before `WriteArchive` (the other formats' formatters cut them).
    let prunes_lists = args.app.prunes_after_traceback()
        && !matches!(requested_formats.first(), Some(BlastpReport::Archive));
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:707-711
    // psiblast's XML formatter cuts the round's results themselves (`TrimSeqAlign`, the subjects
    // that `PruneSeqalign` keeps) before the round's PSSM and inclusion set are made from them
    // (psiblast_app.cpp:273-309); the other formats cut a copy after them.
    let cuts_before_pssm = prunes_lists
        && args.app == crate::algorithm::blastp::args::BlastpApp::Psiblast
        && matches!(requested_formats.first(), Some(BlastpReport::Xml));
    let psiblast_kept: Option<Vec<usize>> = match (args.app, listed_titles.as_ref()) {
        (_, Some(titles)) if prunes_lists => {
            Some(titles.joined_cut(hit_lists, args.max_target_seqs)?)
        }
        _ => None,
    };
    // Whether the batches prefetch the subjects of a database, and whether that looks their
    // taxids up (`prefetch_reported_subjects`'s `prefetch_messages`).
    let batches_prefetch = listed_titles.is_some()
        && requested_formats
            .first()
            .is_some_and(|outfmt| outfmt.prefetches_subjects());
    let prefetch_taxids =
        batches_prefetch && database.is_some_and(|database| database.taxdb_found());
    // Whether that prefetch reads the subjects though the first format reads nothing of them
    // (the Seq-annots with `PRE_FETCH_SEQS_LIMIT=2147483647`, a table of fields that read nothing
    // of the subject): a subject with parsed Seq-ids that it cannot read as NCBI does is refused
    // also without `taxdb.bti` (`ListedTitles::check_prefetched`).
    let prefetch_only = batches_prefetch
        && requested_formats
            .first()
            .is_some_and(|outfmt| !outfmt.reads_subjects());
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:673-678
    // NCBI's prefetch reads the subject of every Seq-align of a batch's results, which are
    // psiblast's before its formatter's cut: the subjects of the traceback hit lists, also those
    // beyond `max_target_seqs` that no report shows. Their ranks, by query, for the prefetch
    // (`ListedTitles::prefetched_ranks`).
    let psiblast_prefetched: Option<Vec<Vec<usize>>> =
        (args.app.prunes_after_traceback() && batches_prefetch).then(|| {
            hit_lists
                .iter()
                .map(|hit_list| {
                    let mut ranks: Vec<usize> =
                        database::ListedTitles::reported_ranks(std::slice::from_ref(hit_list))
                            .collect();
                    ranks.dedup();
                    ranks
                })
                .collect()
        });
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:291-309
    // psiblast's PSSM is made from the Seq-aligns of the round's results, before the
    // formatter's cut; its files are staged here (`psiblast_pssms`) and written after the
    // reports (`run_call`'s `output_result`). A refusal between here and the reports
    // (`refuse_unreportable` in `run_call`) or a report write error drops them: no PSSM file has
    // changed (`PssmFiles`).
    // Whether each query's batch prefetches (`shared::BatchPrefetch`), by the Seq-aligns of its
    // results, one for each HSP (psiblast's before the cut, as `psiblast_hits_found` below): a
    // batch whose prefetch is off reads its subjects one by one, for the PSSM too. The rounds of
    // one psiblast query share NCBI's scope, but every subject that a round whose prefetch is off
    // reads one by one is refused there (its reports and its PSSM, below), so no weak match
    // reaches a later round (review 2 of line DBID weakbatch, L1).
    let batch_ranges: Vec<std::ops::Range<usize>> = query_batches
        .iter()
        .map(|batch| batch.searched.clone())
        .collect();
    let batch_prefetch = crate::algorithm::common::database::BatchPrefetch::of(
        hit_lists.len(),
        0,
        &batch_ranges,
        |span| {
            blastp_hit_list_hsp_count(&hit_lists[span.clone()])
                + setup
                    .unlisted
                    .get(span)
                    .map_or(0, |counts| counts.iter().sum())
        },
    );
    // The numbers of HSPs of psiblast's lists before the formatter's cut (`# N hits found` of
    // outfmt 7), and the cut: here for psiblast's XML report, after the PSSMs otherwise.
    let mut psiblast_hits_found: Option<Vec<usize>> = None;
    let mut cut_lists = |hit_lists: &mut [Option<BlastpHitList>]| {
        psiblast_hits_found = Some(
            hit_lists
                .iter_mut()
                .enumerate()
                .map(|(query, hit_list_opt)| {
                    let found = blastp_hit_list_hsp_count(std::slice::from_ref(hit_list_opt))
                        + setup.unlisted.get(query).copied().unwrap_or(0);
                    if let Some(hit_list) = hit_list_opt {
                        hit_list.prune_by_size(
                            psiblast_kept
                                .as_ref()
                                .map_or(args.max_target_seqs, |kept| kept[query]),
                        );
                    }
                    found
                })
                .collect(),
        );
    };
    if cuts_before_pssm {
        cut_lists(hit_lists);
    }
    // `-out_pssm`'s PSSM, before it is made (a later query's PSSM error reruns the search without
    // the files, review 2 of line DBID weakbatch H2).
    // What the PSSMs read by the subjects' Seq-align ids (`-out_pssm`'s and each round's next
    // one, `ListedTitles::pssm_reads`), refused where one is made: `-out_pssm`'s here, a
    // round's below.
    let pssm_reads = match (
        args.psi_pssm.is_some() || reports.psi_round.is_some(),
        listed_titles.as_ref(),
        database,
    ) {
        // The length of each query's hit list after the formatter's cut (`prune_by_size`, below).
        (true, Some(titles), Some(database)) => titles.pssm_reads(
            database,
            hit_lists,
            &batch_prefetch,
            args.inclusion_ethresh,
            prefetch_taxids,
            &|query| match psiblast_kept.as_ref() {
                _ if !prunes_lists => usize::MAX,
                Some(kept) => kept[query],
                None => args.max_target_seqs,
            },
        )?,
        _ => database::PssmReads::default(),
    };
    if args.psi_pssm.is_some() {
        pssm_reads.check()?;
    }
    let subject_keys = psi_round::SubjectKeys {
        listed_ordinals: listed_subjects
            .as_ref()
            .zip(database)
            .map(|(listed, database)| (listed.ordinals.as_slice(), database.oid_map())),
        parsed: listed_titles
            .as_ref()
            .map_or(&[][..], |titles| titles.parsed.as_slice()),
        reads: &pssm_reads,
    };
    // The search again of the queries up to a PSSM's stop takes the PSSMs made before it
    // (`BlastpReports::stopped_pssms`).
    let stopped_pssms = reports
        .stopped_pssms
        .and_then(|stopped| stopped.lock().expect("the stopped PSSMs").take());
    let psiblast_pssms: Option<crate::algorithm::psiblast::pssm::PssmFiles> =
        match args.psi_pssm.as_ref() {
            Some(_) if stopped_pssms.is_some() => stopped_pssms,
            Some(request) => Some(psiblast_pssms(
                args,
                request,
                hit_lists,
                ranges,
                shown_subjects,
                contexts,
                gapped_params,
                subject_keys,
            )?),
            None => None,
        };
    // psiblast's rounds: the round's outcome (its inclusion set and the next round's PSSM, from
    // the same Seq-aligns), and the subjects of the previous round's set for its two-part table
    // (`psi_round::round_outcome`).
    let psi_found_before: Option<psi_round::FoundBefore> = match reports.psi_round {
        Some(io) => psi_round::round_outcome(
            args,
            io,
            hit_lists,
            ranges,
            shown_subjects,
            contexts,
            gapped_params,
            subject_keys,
        )?,
        None => None,
    };
    // The PSSM of each round that makes the next round's PSSM (review 1 of line DBID weakbatch,
    // H1: psiblast_app.cpp:291-309 `ComputePssmForNextIteration`) reads the letters of the
    // subjects before the cut by their Seq-align ids.
    if let Some(io) = reports.psi_round {
        if io.makes_next() {
            pssm_reads.check()?;
        }
    }
    // Whether this call makes a PSSM: `-out_pssm`'s, or a round's for the next round.
    let pssm_made = args.psi_pssm.is_some() || reports.psi_round.is_some_and(|io| io.makes_next());
    if prunes_lists && !cuts_before_pssm {
        cut_lists(hit_lists);
    }
    let mut engine_counts = EngineCounts::default();
    recount_shown_identities(
        hit_lists,
        ShownQueries::of(query_records, args),
        shown_subjects,
        args.scoring.matrix,
        // `-db_hard_mask`: the search read `X` over the masks of the database sequences, and the
        // reports count their stored letters (series NP, M4b).
        psi_round::searched_other_letters(args)
            || args
                .db_mask_resolved
                .as_ref()
                .is_some_and(|mask| mask.kind == crate::blastinput::db_mask::MaskKind::Hard),
        &mut engine_counts,
    );
    Ok(ReportValues {
        batches_prefetch,
        prefetch_taxids,
        prefetch_only,
        psiblast_prefetched,
        batch_ranges,
        batch_prefetch,
        pssm_reads,
        psiblast_pssms,
        psi_found_before,
        pssm_made,
        psiblast_hits_found,
        engine_counts,
    })
}

/// The labels that the reports show of the subjects (`report_labels`), each `run_call`'s value of
/// the same name: the subject tag (`subject_input_label`), the numbers of sequences and letters of
/// the `Database:` lines (NCBI's `Int8`), the database's sums as NCBI's `Int8` values (`shown_database`, `None` for
/// a subject set: the XML report's `db-num`, `db-len` and their overflow), and the context of the
/// tabular reports.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:138-139
pub(super) struct ReportLabels {
    pub(super) subject_input_label: String,
    pub(super) shown_num_sequences: i64,
    pub(super) shown_total_letters: i64,
    pub(super) shown_database: Option<crate::report::pairwise::PairwiseDatabase>,
    pub(super) context: ReportContext,
}

/// The labels that the reports show of the subjects: the subject tag, the numbers of sequences and
/// letters of the `Database:` lines (`shown_records`: the records that the reports show,
/// `shown_of`'s), and the context of the tabular reports, with their `# Database:` line.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:120-122
pub(super) fn report_labels(
    setup: &ReportSetup<'_>,
    shown_records: &[FastaRecord],
) -> ReportLabels {
    let ReportSetup {
        args,
        database,
        subject_label,
        total_db_len,
        ..
    } = *setup;
    let subject_name_label = if subject_label.is_empty() {
        path_label(&args.subject)
    } else {
        subject_label.to_string()
    };
    let subject_input_label = blastp_subject_tag(args, subject_label);
    // The number of sequences that the `Database:` lines show: the sums over the names of the
    // `-db` list (`PairwiseDatabase`, a name given twice twice), or the records of the subject
    // set. With a taxid list the reports show the whole database's totals and the statistics
    // read the list's (line DI-8).
    let shown_database = database.map(|database| database.pairwise_database());
    let shown_num_sequences = match &shown_database {
        Some(shown) => shown.num_sequences,
        None => shown_records.len() as i64,
    };
    let shown_total_letters = match &shown_database {
        Some(shown) => shown.total_letters,
        None => total_db_len,
    };
    let context = ReportContext {
        subject_name: Some(subject_name_label.clone()),
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:797-801
        // The `# Database:` line of outfmt 7 shows the `-db` value as NCBI keeps it.
        database_name: Some(match database {
            Some(database) => database.name().to_string(),
            None => format!("User specified sequence set (Input: {subject_input_label})"),
        }),
        program: args.app.ncbi_name().to_string(),
    };
    ReportLabels {
        subject_input_label,
        shown_num_sequences,
        shown_total_letters,
        shown_database,
        context,
    }
}
