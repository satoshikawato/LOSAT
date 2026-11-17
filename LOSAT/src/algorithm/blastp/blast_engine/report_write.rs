//! The BLAST archive, the tabular reports and the organism report of a BLASTP call's final hits
//! (`run_call`'s loop over the requested formats): outfmt 11 (`write_archive`); the tabular
//! reports with a field list or a delimiter, outfmt 10 and 20 and the organism report of
//! outfmt 18 (`write_custom_tabular`); and the default outfmt 6 and 7 (`write_text`). Each
//! function is the statements of one arm of the loop's `match` over the format, and `run_call`
//! calls it in that arm: its error ends the loop as the arm's did.
//!
//! The NCBI references are at the functions: `CBlastFormat::WriteArchive` at `write_archive`, the
//! tabular writer of `CBlastFormat::x_PrintTabularReport` at `write_custom_tabular`, and NCBI's
//! default tabular fields at `write_text`.
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`BlastpReports`, `BlastpRanges`, `query_record_masks`, `write_blastp_tabular_hit_lists`) and
//! its sibling modules `database` and `query_rows`.

use std::io::Write;
use std::ops::Range;
use std::sync::Arc;

use anyhow::Result;

use super::{
    database, query_record_masks, query_rows, write_blastp_tabular_hit_lists, BlastpRanges,
    BlastpReports,
};
use crate::algorithm::blastp::args::{BlastpApp, ResolvedBlastpArgs};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::BlastpHitList;
use crate::algorithm::blastp::seqalign::EngineCounts;
use crate::algorithm::blastp::{archive, tabular};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::{FormatOutput, FormatProbe};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::db_subject::ParsedSubject;
use crate::report::defline::SourceTaxname;
use crate::report::pairwise::PairwiseDatabase;
use crate::report::query_warnings::QueryWarnings;
use crate::report::tabular::TabularSpec;
use crate::report::{OutputFormat, PairwiseHit, ReportContext};
use crate::utils::threading::SearchPool;

/// What the three functions of this module read of `run_call`, made once before its loop over the
/// formats: each field is `run_call`'s value of the same name, by reference where it is not `Copy`
/// (`subject_records` are the records that the reports show, a database's listed subjects' or the
/// subject set's; `subjects` are the encoded subjects). The format's output, its warnings and its
/// observer are the loop's own values, and `run_call`'s closures that an arm calls
/// (`first_context`, `batch_spans`, `rows_of`) stay closures there: the functions take each of
/// them as a parameter.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:120-128
/// ```c++
///     m_DbName = db_adapter.GetDatabaseName();
///     m_IsBl2Seq = (m_DbName == kEmptyStr ? true : false);
///     m_IsDbScan = db_adapter.IsDbScanMode();
///     if (m_IsBl2Seq) {
///         m_SeqInfoSrc.Reset(db_adapter.MakeSeqInfoSrc());
///     }
///     else {
///     	m_SearchDb = db_adapter.GetSearchDatabase();
///     }
/// ```
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
    pub(super) database: Option<&'a database::BlastpDatabase>,
    pub(super) listed_subjects: &'a Option<database::ListedSubjects>,
    pub(super) listed_titles: &'a Option<database::ListedTitles<'a>>,
    pub(super) listed_parsed: &'a [Option<Arc<ParsedSubject>>],
    pub(super) listed_sources: &'a [SourceTaxname],
    pub(super) pairwise_database: &'a Option<PairwiseDatabase>,
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) subject_ids: &'a [crate::algorithm::common::tabular::SubjectId],
    pub(super) subject_accessions: &'a [Option<[Arc<[u8]>; 2]>],
    pub(super) subject_stitles: &'a [Arc<[u8]>],
    pub(super) hit_lists: &'a [Option<BlastpHitList>],
    pub(super) engine_counts: &'a EngineCounts,
    pub(super) pairwise_hits: &'a Option<Arc<Vec<PairwiseHit>>>,
    pub(super) psiblast_hits_found: &'a Option<Vec<usize>>,
    pub(super) tabular_unsearched: &'a [bool],
    pub(super) context: &'a ReportContext,
    pub(super) call_batches: &'a archive::CallBatches<'a>,
    pub(super) x_rows_pool: Option<&'a SearchPool<'a>>,
}

/// Outfmt 11: one BLAST archive per query batch of the call (`archive::write_blastp_archives`),
/// with the batch's queries, their masks and the warnings of its queries, flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1067-1071
/// ```c++
/// CBlastFormat::WriteArchive(blast::IQueryFactory& queries,
///                            blast::CBlastOptionsHandle& options_handle,
///                            const CSearchResultSet& results,
///                            unsigned int num_iters,
///                            const list<CRef<CBlast4_error> > & msg)
/// ```
pub(super) fn write_archive(
    view: &ReportView<'_>,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    first_context: &impl Fn(usize) -> archive::FirstContext,
    batch_spans: &impl Fn() -> Vec<Range<usize>>,
) -> Result<()> {
    let ReportView {
        args,
        ranges,
        query_frames,
        query_lengths,
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
    // ```c++
    //             if (fmt_args->ArchiveFormatRequested(args)) {
    //                 formatter.WriteArchive(*queries, *m_OptsHndl, *results,  0, m_Bah.GetMessages());
    //                 m_Bah.ResetMessages();
    // ```
    // One BLAST archive per query batch of the call (`archive`).
    let handle = args
        .archive
        .as_ref()
        .expect("a BLAST archive request has its handle");
    let mut writer = format.sink.open()?;
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
    let masks = query_record_masks(query_frames, &ranges.query_input.placements);
    let spans = batch_spans();
    let call = archive::ArchiveCall {
        handle,
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
        first_context,
        &mut writer,
    )?;
    // Flushed when written, as the other reports are.
    writer.flush()?;
    Ok(())
}

/// The tabular reports with a field list or a delimiter and outfmt 10 and 20
/// (`tabular::write_blastp_tabular`), and the organism report of outfmt 18
/// (`tabular::write_blastp_tax_report`), into the format's stream: the taxonomy and the names that
/// their fields show are read before the report's first byte, and the rows of the HSPs a window of
/// queries at a time where a field shows them (`rows_of`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:782-785
/// ```c++
///         CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
///         if(!m_CustomDelim.empty()) {
///             tabinfo.SetCustomDelim(m_CustomDelim);
///         }
/// ```
///
/// The organism report of outfmt 18 is NCBI's `CBlastFormat::x_PrintTaxReport`:
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:873-874
/// ```c++
/// void
/// CBlastFormat::x_PrintTaxReport(const blast::CSearchResults& results)
/// ```
pub(super) fn write_custom_tabular(
    view: &ReportView<'_>,
    spec: &'static TabularSpec,
    format: &mut FormatOutput<'_>,
    warnings: &mut Option<QueryWarnings<'_>>,
    probe: &mut Option<FormatProbe<'_>>,
    rows_of: &impl Fn(Range<usize>) -> Result<Vec<PairwiseHit>>,
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
        tabular_unsearched,
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
    // ```c++
    //        		const CBioseq_Handle& subject_bh =
    //                 scope.GetBioseqHandle(align.GetSeq_id(1));
    //             ...
    //             m_SubjectLength = subject_bh.GetBioseqLength();
    // ```
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
            titles.tabular_taxonomy(database, hit_lists)?
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
        args.app.ncbi_name().to_uppercase(),
        crate::NCBI_BLAST_VERSION
    );
    let tab = tabular::BlastpTabular {
        spec,
        app: args.app,
        program_version: &program_version,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1458
        // ```c++
        //         x_PrintTabularReport(results, itr_num);
        // ```
        // psiblast's round (`psi_round`), 1 outside the rounds.
        iteration: (args.app == BlastpApp::Psiblast)
            .then(|| reports.psi_round.map_or(1, |io| io.state.number())),
        matrix: args.scoring.matrix,
        query_ids,
        query_records,
        query_titles,
        searched_lengths: query_lengths,
        range_length: ranges.query_input.range_length,
        database: context
            .database_name
            .as_deref()
            .or(context.subject_name.as_deref())
            .unwrap_or_default(),
        subjects: tabular_subjects,
        unsearched: tabular_unsearched,
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
        )?;
        writer.flush()?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2955-2958
    // ```c++
    //     	if (args[kArgLineLength]) {
    //    		 ERR_POST(Warning << "The parameter -line_length is not applicable for "
    //    				    		 "output formats > 4 .");
    //     	}
    // ```
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

/// The default outfmt 6 and 7 (`write_blastp_tabular_hit_lists`): NCBI's default fields, with the
/// comment lines of outfmt 7, flushed when written.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/format_flags.cpp:38-40
/// ```c++
/// const char* kDfltArgTabularOutputFmt =
///     "qaccver saccver pident length mismatch gapopen qstart qend sstart send "
///     "evalue bitscore";
/// ```
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
        tabular_unsearched,
        context,
        ..
    } = *view;
    let mut writer = format.sink.open()?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1102-1107
    // ```c++
    //     ITERATE(list<ETabularField>, iter, m_FieldsToShow) {
    //         // Add tab in front of field, except for the first field.
    //         if (iter != m_FieldsToShow.begin())
    //             m_Ostream << m_FieldDelimiter;
    //         x_PrintField(*iter);
    //     }
    // ```
    write_blastp_tabular_hit_lists(
        hit_lists,
        outfmt,
        &mut writer,
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1306-1308
        // ```c++
        //     CAlignFormatUtil::AcknowledgeBlastQuery(bioseq, kLineLength, m_Ostream,
        //                                             m_ParseLocalIds, kHtmlFormat,
        //                                             kTabularFormat, rid);
        // ```
        query_ids,
        query_titles,
        subject_ids,
        subject_accessions,
        tabular_unsearched,
        psiblast_hits_found.as_deref(),
        context,
        probe.as_mut(),
        warnings.as_mut(),
        reports.epilog,
        reports.x_queries_before,
        // psiblast's round (`psi_round`), 1 outside the rounds.
        (args.app == BlastpApp::Psiblast)
            .then(|| reports.psi_round.map_or(1, |io| io.state.number())),
    )?;
    // NCBI reference (598d8ae6): c++/src/corelib/ncbidiag.cpp:4086
    // ```c++
    //         CDiagHandler* handler = new CStreamDiagHandler(&NcbiCerr, true, kLogName_Stderr);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:134-146
    // ```c++
    // CRef<CBlastQueryVector>
    // CBlastInput::GetNextSeqBatch(CScope& scope)
    // {
    //     CRef<CBlastQueryVector> retval(new CBlastQueryVector);
    //     TSeqPos size_read = 0;
    //
    //     while (size_read < GetBatchSize()) {
    //
    //         if (End())
    //             break;
    //
    //         CRef<CBlastSearchQuery> q;
    //         try { q.Reset(m_Source->GetNextSequence(scope)); }
    // ```
    // NCBI reads the next query batch (blastp_app.cpp:257-259) after this
    // batch's report, and the reader's messages and error go to `cerr`, which
    // the C++ library ties to `cout`: the rows written so far reach standard
    // output before them. The report is flushed when it is written, as the
    // pairwise report is and as TBLASTN's and TBLASTX's are, so a later
    // batch's stop messages (`search`, `batches.stop`) follow the rows of the
    // batches before it (audit findings A-3 and C-1 of session SFc).
    writer.flush()?;
    Ok(())
}
