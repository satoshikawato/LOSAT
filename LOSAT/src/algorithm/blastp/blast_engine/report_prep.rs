//! What a BLASTP call does with its final hits before any format is written (`run_call`, after
//! the rows of the reports): the warning lines of each query (`query_warning_lines`), the checks
//! of the reported subjects and each query batch's prefetch of them
//! (`prefetch_reported_subjects`), and the refusals of what LOSAT does not report as NCBI does
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
//! `blastp_hit_list_hsp_count`, `query_record_masks`) and its sibling modules `database`,
//! `psi_round` and `pssm_hooks`.

use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::{
    blastp_hit_list_hsp_count, blastp_pairwise_view_of, database, psi_round, pssm_hooks,
    query_record_masks, BlastpRanges, BlastpReport, BlastpReports,
};
use crate::algorithm::blastp::args::ResolvedBlastpArgs;
use crate::algorithm::blastp::hsp::{collect_hits_and_methods_from_hit_lists, BlastpHitList};
use crate::algorithm::blastp::seqalign::EngineCounts;
use crate::algorithm::common::database::BatchPrefetch;
use crate::algorithm::psiblast::pssm::PssmFiles;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::report::pairwise::PairwiseHit;
use crate::report::query_warnings::SubjectReads;

/// What the three functions of this module read of `run_call`: each field is `run_call`'s value of
/// the same name, by reference where it is not `Copy` (`outputs` for its formats only;
/// `subject_records` and `subject_titles` are those that the reports show, a database's listed
/// subjects' or the subject set's; `format_plans` is the requested formats). The values that
/// `run_call` makes for the prefetch (`batches_prefetch` to `batch_prefetch`) and for the PSSMs
/// (`pssm_made`, `pssm_reads`) come last.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:756-758
/// ```c++
/// void BlastFormatter_PreFetchSequenceData(const blast::CSearchResultSet& results,
/// 		                                 CRef<CScope> scope,
/// 		                                 blast::CFormattingArgs::EOutputFormat format_type)
/// ```
pub(super) struct ReportPrep<'a, 'o> {
    pub(super) args: &'a ResolvedBlastpArgs,
    pub(super) outputs: &'a ReportOutputs<'o>,
    pub(super) reports: &'a BlastpReports<'a>,
    pub(super) ranges: &'a BlastpRanges<'a>,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) query_frames: &'a [Vec<QueryFrame>],
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
/// ```c++
///     string retval(m_Errors.GetQueryId());
///     if ( !retval.empty() ) {    // in case the query id is not known
///         retval += ": ";
///     }
///     ITERATE(TQueryMessages, iter, m_Errors) {
///         if ((**iter).GetSeverity() == eBlastSevWarning) {
///             retval += (*iter)->GetMessage(false) + " ";
///         }
///     }
/// ```
pub(super) fn query_warning_lines(prep: &ReportPrep<'_, '_>) -> Vec<Vec<u8>> {
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
    // the warnings go with the first format.
    let result_warnings = outputs
        .formats
        .first()
        .and_then(crate::api::local_blast::FormatOutput::report_format)
        .is_none_or(crate::blastinput::app::ReportFormat::writes_query_warnings);
    query_records
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
            if result_warnings {
                lines.extend(if let Some(psi) = args.psi_search.as_deref() {
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
                });
            }
            lines
        })
        .collect()
}

/// Check the reported subjects before any format is written and read each query batch's prefetch
/// of them: the database titles and the subjects with parsed Seq-ids that the reports of the
/// batches before a prefetch's stop show and LOSAT does not write as NCBI does are refused, and
/// the stopping batch's prefetch returns NCBI's error. The texts of a batch's prefetch go to the
/// warning lines of its first query, after the reader's messages of the batch; where a batch's
/// prefetch is off, the texts of the subjects that a query's report reads go to the query's
/// warning lines, or to the returned `SubjectReads` where the first format places them
/// (`ReportFormat::places_subject_reads`), and the texts of the loads of a query's PSSM go to the
/// PSSM files' read messages (`psiblast_pssms`) or to the round's outcome.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:767-772
/// ```c++
///     try {
///        CScope::TIds ids;
///        vector<TSeqRange> ranges;
///        s_ExtractSeqidsAndRanges(results, ids, ranges);
///        _TRACE("Prefetching " << ids.size() << " sequence lengths");
///        LoadSequencesToScope(ids, ranges, scope);
/// ```
pub(super) fn prefetch_reported_subjects(
    prep: &ReportPrep<'_, '_>,
    query_warning_lines: &mut [Vec<u8>],
    psiblast_pssms: &mut Option<PssmFiles>,
) -> Result<SubjectReads> {
    let ReportPrep {
        args,
        outputs,
        reports,
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
        ..
    } = *prep;
    // A database title that LOSAT does not report as NCBI does is an error before any
    // format is written.
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288 (quoted below)
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
    if let Some(titles) = listed_titles.as_ref() {
        titles.check_reported_titles(reported_lists)?;
        // The subjects with parsed Seq-ids that the reports asked for do not write as NCBI
        // (the formats that read a subject's Bioseq need its sequence found by its Seq-align
        // id, which a table of fields that read nothing of the subject does not; outfmt 0, 5 and
        // 13 to 16 join the hits of matching Seq-align ids, which LOSAT does for copies of one
        // record; 13 to 16 and the tables of taxonomy fields write their taxids).
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/psi_pssm_input.cpp:440-450
        // ```c++
        //     CSeq_loc seqloc(const_cast<CSeq_id&>(*ds.GetIds().back()), subj_start,
        //                     subj_start+subjlen-1);
        //
        //     try {
        //         CSeqVector sv(seqloc, scope);
        //         sv.SetCoding(CSeq_data::e_Ncbistdaa);
        //         sv.GetSeqData(0, kInvalidSeqPos, sequence_data);
        //     } catch (const CException&) {
        //         sequence_data.erase();
        //         ERR_POST(Warning << "Failed to retrieve sequence " <<
        //                  seqloc.GetInt().GetId().AsFastaString());
        // ```
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
            format_plans.iter().fold(
                crate::report::tabular::PrintedIds::default(),
                |words, outfmt| words.or(outfmt.printed_ids()),
            ),
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:284-288
    // ```c++
    //                 BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		                            fmt_args->GetFormattedOutputChoice());
    //                 ITERATE(CSearchResultSet, result, *results) {
    //                     formatter.PrintOneResultSet(**result, query_batch);
    //                 }
    // ```
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
    // ```c++
    //                    results = psiblast->Run();
    //                    BlastFormatter_PreFetchSequenceData(*results, scope,
    //                 		   m_CmdLineArgs->GetFormattingArgs()->GetFormattedOutputChoice());
    // ```
    // psiblast's batches, of one query each (`read_psiblast_query_batches`), prefetch the same
    // way, the subjects of the traceback hit lists (`psiblast_prefetched`).
    // Where `PRE_FETCH_SEQS_LIMIT` turns a batch's prefetch off (`s_PreFetchSeqs`, quoted at
    // `shared::prefetch_may_be_off`), each query's report fetches the subjects it shows, after
    // the query's own warnings: those of its final hit list, psiblast's after the cut, where the
    // first format's rows read their Bioseqs (`TabularSpec::subject_fetch`; `nident` alone on
    // rows with identities, `ListedTitles::fetched_ranks`).
    // The texts of the subjects that the reports of the queries whose batch's prefetch is off
    // read one by one, by query, where the first format places them
    // (`ReportFormat::places_subject_reads`).
    let mut subject_reads = crate::report::query_warnings::SubjectReads::default();
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
    if let (Some(titles), Some(database)) = (listed_titles.as_ref(), database) {
        // Without `taxdb.bti` the prefetch of a table that reads the subjects writes nothing
        // before it stops: the stopping batch's prefetch alone (`ListedTitles::prefetch_stop`).
        // The XML2 and JSON formats: the stop alone (`shared::batch_prefetch_stop`), whose
        // other subjects' loads write their texts at each `Hit`.
        if let Some(span) = prefetch_stop.as_ref().filter(|_| xml2_prefetch) {
            titles.batch_prefetch_stop(database, prefetched_ranks(span.clone()), span.start)?;
        } else if let Some(span) = prefetch_stop
            .as_ref()
            .filter(|_| !(prefetch_taxids || psiblast_prefetched.is_some() || prefetch_only))
        {
            titles.batch_prefetch(
                database,
                database::ListedTitles::reported_ranks(&hit_lists[span.clone()]),
                span.start,
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
            let rows = match format_plans.first() {
                Some(BlastpReport::CustomTabular(spec)) => {
                    spec.subject_fetch(crate::blastinput::app::Program::Blastp)
                }
                _ => crate::report::tabular::SubjectFetch::Every,
            };
            let open_for_rows = |lists: &[Option<BlastpHitList>]| -> Result<()> {
                if plural_names
                    && lists
                        .iter()
                        .flatten()
                        .any(|hit_list| hit_list.hsplist_count > 0)
                {
                    database.open_for_tabular_rows(args.app.name())?;
                }
                Ok(())
            };
            for batch in query_batches {
                let lists = &hit_lists[batch.searched.clone()];
                // The batch's Seq-aligns, one for each HSP of its results.
                let seq_aligns = match psiblast_hits_found.as_ref() {
                    Some(found) => found[batch.searched.clone()].iter().sum(),
                    None => blastp_hit_list_hsp_count(lists),
                };
                if crate::algorithm::common::database::prefetch_off(seq_aligns) {
                    if !prefetch_taxids {
                        continue;
                    }
                    for query in batch.searched.clone() {
                        if rows != crate::report::tabular::SubjectFetch::Never {
                            let shown = std::slice::from_ref(&hit_lists[query]);
                            // Rows that read a subject only for its identities (`nident` alone)
                            // read none of a query whose HSPs have none, which looks no taxid up.
                            let mut fetched =
                                database::ListedTitles::fetched_ranks(shown, rows).peekable();
                            if rows == crate::report::tabular::SubjectFetch::Every
                                || fetched.peek().is_some()
                            {
                                // The first format writes the texts where it reads each
                                // subject (`places_subject_reads`), else before the report.
                                if places_reads {
                                    // The pairwise report reads its subjects in the order of
                                    // its sorted groups (`QueryGroup::read_order` of the rows'
                                    // hits, `-sorthits`), the other formats in their rows' order.
                                    let ranks: Vec<usize> = if pairwise_first {
                                        let hits = collect_hits_and_methods_from_hit_lists(
                                            shown,
                                            query,
                                            engine_counts,
                                        );
                                        crate::report::x_pairwise_par::QueryGroup::read_order(
                                            hits.into_iter().map(
                                                |(hit, _, seqalign_num_ident, _)| {
                                                    (hit.s_idx, (hit, seqalign_num_ident))
                                                },
                                            ),
                                            |(hit, seqalign_num_ident)| {
                                                crate::report::hit_sort::HspScores::of_hit(
                                                    hit,
                                                    *seqalign_num_ident,
                                                )
                                            },
                                            |previous, s_idx| {
                                                crate::report::db_subject::rank_joins_as_copy(
                                                    &titles.parsed,
                                                    previous as usize,
                                                    s_idx as usize,
                                                )
                                            },
                                            args.sort,
                                            false,
                                        )
                                        .into_iter()
                                        .map(|rank| rank as usize)
                                        .collect()
                                    } else {
                                        fetched.collect()
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
                        if pssm_made {
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
                let messages = match psiblast_prefetched.as_ref() {
                    Some(ranks) => {
                        let prefetched = titles.prefetched_ranks(
                            database,
                            &ranks[batch.searched.clone()].concat(),
                            prefetch_taxids,
                            !prefetch_only,
                        )?;
                        titles.batch_prefetch(database, prefetched, batch.searched.start)?
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
                            database::ListedTitles::reported_ranks(lists),
                            batch.searched.start,
                        )?
                    }
                };
                open_for_rows(lists)?;
                if messages.is_empty() {
                    continue;
                }
                let first = batch.searched.start;
                let reader = reports.reader_warnings.get(first).map_or(0, Vec::len);
                query_warning_lines[first].splice(reader..reader, messages);
            }
        }
    }
    Ok(subject_reads)
}

/// Refuse, before any format of the call is written, what LOSAT does not report as NCBI does:
/// with outfmt 0 and `-parse_deflines`, a query that is also a subject and has masks; with
/// outfmt 1 to 4 of a `-subject` search, a query whose own record is one of the subjects and
/// whose query-anchored section NCBI aborts on (the error is that abort's); and with outfmt 1 to
/// 4 of a database search, or outfmt 0 to 4 with the queries' ids, a reported sequence that shows
/// a query's id of the batch that reports it. `rows_of` makes the rows of a window of queries,
/// where the call has not made every row before (`pairwise_hits`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1412
/// ```c++
/// CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
///                         CConstRef<blast::CBlastQueryVector> queries,
/// ```
pub(super) fn refuse_unreportable(
    prep: &ReportPrep<'_, '_>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
) -> Result<()> {
    let ReportPrep {
        args,
        reports,
        ranges,
        query_records,
        query_frames,
        query_batches,
        database,
        listed_titles,
        subject_records,
        subject_titles,
        format_plans,
        hit_lists,
        pairwise_hits,
        psi_found_before,
        ..
    } = *prep;
    // Outfmt 0 only: the query-anchored sections draw these masks as NCBI does
    // (`query_anchored::fill_loc_list`).
    if args.parse_deflines
        && format_plans.iter().any(|outfmt| {
            matches!(
                outfmt,
                BlastpReport::Pairwise(crate::report::pairwise::PairwiseView::Pairwise)
            )
        })
    {
        let masks = query_record_masks(query_frames, &ranges.query_input.placements);
        crate::blastinput::fasta_reader::believed::refuse_shared_masks(
            query_records,
            |query| masks.get(query).is_some_and(|m| !m.is_empty()),
            args.lcase_masking,
            args.app.name(),
        )?;
    }
    // Outfmt 1 to 4 of a `-subject` search: a query whose own record is one of the subjects
    // (`-parse_deflines`) and whose query-anchored section NCBI BLAST+ aborts on is refused
    // before any format of the call is written (`pairwise::blastp_section_aborts`; line RG T3a,
    // PD-LOSAT-NCBI-DEFECTS rule 3): in a block search (`x_query_block_plan`), the call of the
    // query's block, after the bytes of the blocks before, as outfmt 0's refusal of the masks of
    // such a record above (NCBI's own bytes before its abort depend on its buffers, so no prefix
    // is NCBI's). The formats 1 to 4 merge one section (their views differ in the letters they
    // show).
    let anchored_view = format_plans.iter().find_map(|outfmt| match outfmt {
        BlastpReport::Pairwise(view)
            if *view != crate::report::pairwise::PairwiseView::Pairwise =>
        {
            Some(*view)
        }
        _ => None,
    });
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
                psi_found_before.is_some(),
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
    if let (Some(_), false, Some(titles)) =
        (anchored_view, args.parse_deflines, listed_titles.as_ref())
    {
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
        titles.refuse_generated_query_ids(
            hit_lists,
            |query| {
                scopes
                    .get(scopes.partition_point(|&(end, _)| end <= query))
                    .map_or(&none, |(_, scope)| scope)
            },
            |q_idx| described(q_idx, true),
        )?;
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
