//! TBLASTX's reports: the order of the subjects, the prologs, and every output format of a run.
//!
//! This file is a child module of `blast_engine/run_impl.rs` (declared there with `#[path]`, as
//! `x_seed_release.rs` is), so that it reaches the run's types without changing them.
//!
//! NCBI reference: `c++/src/algo/blast/format/blast_format.cpp:68-1411` (`CBlastFormat`), the order
//! of the hits of `c++/src/algo/blast/core/blast_hits.c:3562-3567`, the prologs of
//! `c++/src/app/blast/tblastx_app.cpp:173-209`, and the subject ids of
//! `c++/src/objtools/align_format/tabular.cpp:863-870`.

use super::{QueryRuns, ShownSubjects};
use crate::algorithm::tblastx::args::TblastxArgs;
use crate::algorithm::tblastx::blast_engine::{
    report, FormatProbe, OutputSink, ReportOutputs, TblastxHsp, TblastxQueryStats,
};
use crate::algorithm::tblastx::database::DatabaseTitles;
use crate::algorithm::tblastx::subjects::SubjectSetNumbers;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::utils::genetic_code::GeneticCode;
use anyhow::Result;
use std::sync::Arc;

/// ABI v1's rejection (`V1Checks`) of the outfmt 0 titles that it did not write as NCBI,
/// of the subjects that the reports show: those with hits in the final hit lists, which
/// NCBI's description tables list (blast_format.cpp:1540, `report/defline.rs`). The check
/// comes after the outfmt 0 prolog, where NCBI writes the titles of the first query's
/// report. ABI v1's titles are its `bio` deflines (`web_api::v1_bio::from_bio`), which
/// `check_report_titles` keeps ASCII.
fn check_shown_subject_titles(hits: &[TblastxHsp], subjects: &[FastaRecord]) -> Result<()> {
    let shown: std::collections::BTreeSet<usize> =
        hits.iter().map(|hsp| hsp.hit.s_idx as usize).collect();
    for index in shown {
        crate::report::defline::check_shown_subject_title(
            &String::from_utf8_lossy(&subjects[index].title),
            index + 1,
            "TBLASTX",
        )?;
    }
    Ok(())
}

/// What every report of a run needs besides its HSPs.
pub(crate) struct TblastxReportRun<'a> {
    pub args: &'a TblastxArgs,
    /// The queries reported (those of the batches searched so far) and the subjects.
    pub queries: &'a [FastaRecord],
    pub subjects: &'a [FastaRecord],
    /// The subjects' letters that the displayed rows read (`-db`: the stored codes of the
    /// reported sequences; their records hold none).
    pub subject_letters: report::SubjectLetters<'a>,
    /// `-db`: the database and its reported sequences (`subjects` are their records, by rank).
    pub database: Option<ReportDatabase<'a>>,
    /// The tabular query IDs (`FastaRecord::shown_id`), by query.
    pub query_ids: &'a [Arc<[u8]>],
    /// The per-query statistics and whether the query's batch was searched.
    pub query_stats: &'a [TblastxQueryStats],
    pub searched: &'a [bool],
    /// The warnings written before each query's report (`QueryWarnings`): the reader's
    /// messages of the query's batch (before its first query) and the query's own.
    pub warnings: &'a [Vec<u8>],
    /// The queries' input and ranges.
    pub ranges: &'a QueryRuns<'a>,
    /// Whether the reports end with the epilog (not when the search stops at a batch).
    pub epilog: bool,
    /// The subject set's numbers (the BLAST XML computes a query's length adjustment again
    /// from them).
    pub numbers: &'a SubjectSetNumbers,
    /// `-db`, XML2 and JSON: the taxid of each reported sequence by rank and what the taxonomy
    /// database says about it (empty for the other formats).
    pub xml_taxa: &'a crate::algorithm::tblastx::database::XmlTaxa,
    /// Whether ABI v1's rejection of outfmt 0 titles applies (`V1Checks`).
    pub v1_titles: bool,
    /// EXPERIMENT (LOSAT_X_TBXRPTPAR): the pool the report stage runs on (tblastx/x_rpt_par.rs),
    /// `None` for the serial stage.
    pub x_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool that formats the outfmt 0 text
    /// (`report::x_pairwise_par`), `None` for the serial writer.
    pub x_pairwise_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
}

/// The database of a `-db` search as its reports show it.
#[derive(Clone, Copy)]
pub(crate) struct ReportDatabase<'a> {
    /// The `-db` value (outfmt 7's `# Database:` line).
    pub name: &'a str,
    /// The `Database:` lines of outfmt 0.
    pub pairwise: &'a crate::report::pairwise::PairwiseDatabase,
    /// The reported sequences by rank (their tabular ids).
    pub reported: &'a crate::algorithm::tblastx::database::ReportedSubjects,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3562-3567
// ```c
//    if (!results->hitlist_array[hsp_list->query_index]) {
//        results->hitlist_array[hsp_list->query_index] =
//            Blast_HitListNew(hitlist_size);
//    }
//    Blast_HitListUpdate(results->hitlist_array[hsp_list->query_index],
//                        hsp_list);
// ```
/// The HSPs of `queries` in report order: query by query, the hit lists of `-max_target_seqs`
/// (default 500) in NCBI's order (`report::final_hit_order`), after `-culling_limit`
/// (`report::culled_hit_order`). Each query's order depends on its own HSPs only, so the order of
/// a query batch is that query's order in the whole run.
pub(super) fn report_order(
    args: &TblastxArgs,
    queries: &[FastaRecord],
    hits: Vec<TblastxHsp>,
) -> Vec<TblastxHsp> {
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:2960-2962,2978
    // ```c
    //     	if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //     		hitlist_size = args[kArgMaxTargetSequences].AsInteger();
    //     	}
    //     ...
    //     opt.SetHitlistSize(hitlist_size);
    // ```
    // The hit list size is `args.report_counts.hitlist_size` (`blastinput::report_counts`).
    // NCBI reference: c++/src/algo/blast/api/blast_options_local_priv.hpp:1320-1324
    // ```c
    // CBlastOptionsLocal::SetCullingLimit(int s)
    // {
    //     if (s <= 0) {
    //         return;
    //     }
    // ```
    // A culling limit above 0 installs NCBI's culling writer and pipe
    // (`report::culled_hit_order`).
    let hitlist_size = args.report_counts.hitlist_size;
    // NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:313-319
    // ```c
    //     if (args.Exist(kArgQueryCovHspPerc) && args[kArgQueryCovHspPerc]) {
    //         opt.SetQueryCovHspPerc(args[kArgQueryCovHspPerc].AsDouble());
    //     }
    //
    //     if (args.Exist(kArgMaxHSPsPerSubject) && args[kArgMaxHSPsPerSubject]) {
    //         opt.SetMaxHspsPerSubject(args[kArgMaxHSPsPerSubject].AsInteger());
    //     }
    // ```
    let query_cov_hsp_perc = args.qcov_hsp_perc.unwrap_or(0.0);
    let query_lengths: Vec<usize> = if args.culling_limit > 0 || query_cov_hsp_perc > 0.0 {
        queries.iter().map(|query| query.seq().len()).collect()
    } else {
        Vec::new()
    };
    let filters = report::HitFilters {
        max_hsps: args.max_hsps.unwrap_or(0),
        query_cov_hsp_perc,
        query_lengths: &query_lengths,
    };
    if args.culling_limit > 0 {
        report::culled_hit_order(
            hits,
            hitlist_size,
            args.culling_limit,
            &query_lengths,
            &filters,
        )
    } else {
        report::final_hit_order(hits, hitlist_size, &filters)
    }
}

/// The `Database:` text of the reports of a local `-subject` search.
///
/// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:797-799
/// ```c++
/// 	    string dbname;
/// 	    if (m_IsDbScan)
/// 		dbname = string("User specified sequence set (Input: ") + m_SubjectTag + string(")");
/// ```
fn tblastx_database_name(args: &TblastxArgs) -> String {
    format!(
        "User specified sequence set (Input: {})",
        args.subject_label()
    )
}

/// Whether the outfmt 0 prolog of `sink` is written before the search, as NCBI writes it
/// (`write_pairwise_prologs`). A file sink is created when its report is written
/// (`OutputSink::open`), so its prolog stays at the start of the report.
fn prolog_before_search(sink: &OutputSink<'_>) -> bool {
    !matches!(sink, OutputSink::File(_))
}

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:173-178
// ```c
//         formatter.PrintProlog();
//
//         /*** Process the input ***/
//         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//
//             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:118-119
// ```c
// {
//     m_Outfile.exceptions(NcbiBadbit);
// ```
/// Writes and flushes the outfmt 0 prolog of every outfmt 0 sink whose prolog comes before
/// the search (`before_search`), before NCBI reads the first query batch (a failed write
/// stops the run there), or, for a run that stops before any report, of the others.
pub(super) fn write_pairwise_prologs(
    outputs: &mut ReportOutputs<'_>,
    args: &TblastxArgs,
    subjects: &ShownSubjects<'_>,
    before_search: bool,
) -> Result<()> {
    for format in outputs.formats.iter_mut() {
        if report::output_format(format.outfmt) == report::TblastxOutputFormat::Pairwise
            && prolog_before_search(&format.sink) == before_search
        {
            let mut writer = format.sink.open()?;
            match subjects {
                ShownSubjects::Fasta(records) => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        report::NCBI_TBLASTX_VERSION,
                        &tblastx_database_name(args),
                        records.len(),
                        records.iter().map(|r| r.seq().len()).sum(),
                        None,
                    )?
                }
                ShownSubjects::Database { pairwise, .. } => {
                    crate::report::pairwise::write_tblastx_pairwise_prolog(
                        &mut writer,
                        report::NCBI_TBLASTX_VERSION,
                        "",
                        0,
                        0,
                        Some(pairwise),
                    )?
                }
            }
            writer.flush()?;
        }
    }
    Ok(())
}

// NCBI reference: ncbi-blast/c++/src/app/blast/tblastx_app.cpp:203-209
// ```c
//                 ITERATE(CSearchResultSet, result, *results) {
//                     formatter.PrintOneResultSet(**result, query_batch);
//                 }
//             }
//         }
//
//         formatter.PrintEpilog(opt);
// ```
// NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1574-1577
// ```c
// // Sort HSPs with e-values as first priority and scores as
// // tie-breakers, since that is the order we want to see them in
// // in Seq-aligns.
// Blast_HSPListSortByEvalue(hsp_list);
// ```
/// Writes the final hits (already filtered by e-value) to every requested format, and
/// gives the caller the same final list (`ReportOutputs::hits`).
fn write_tblastx_outputs(
    mut hits: Vec<TblastxHsp>,
    outputs: &mut ReportOutputs<'_>,
    run: &TblastxReportRun<'_>,
) -> Result<()> {
    let args = run.args;
    let query_code = GeneticCode::from_id(args.query_gencode);
    let db_code = GeneticCode::from_id(args.db_gencode);
    let output_formats: Vec<report::TblastxOutputFormat> = outputs
        .formats
        .iter()
        .map(|format| report::output_format(format.outfmt))
        .collect();
    if run.v1_titles && output_formats.contains(&report::TblastxOutputFormat::Pairwise) {
        check_shown_subject_titles(&hits, run.subjects)?;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1199-1220
    // ```c++
    //     // Set the identity score
    //     if (hsp->num_ident >= 0) {
    //         static const string kNumIdent("num_ident");
    //         scores.push_back(s_MakeScore(kNumIdent, 0.0, hsp->num_ident, true));
    //     }
    //     ...
    //     if (hsp->num_positives > 0) {
    //             static const string kNumPositives("num_positives");
    //             scores.push_back(s_MakeScore(kNumPositives, 0.0, hsp->num_positives, true));
    //     }
    // ```
    // The Seq-aligns (outfmt 8, 12) show the engine's identities and positives
    // (`report::engine_counts`), read before the display counts replace the identities.
    let seqalign = output_formats
        .iter()
        .any(|format| matches!(format, report::TblastxOutputFormat::Seqalign(_)));
    let engine_counts = if seqalign {
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the counts over chunks of the list on the pool.
        match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::engine_counts(
                pool,
                &hits,
                run.queries,
                run.subject_letters,
                &query_code,
                &db_code,
            ),
            None => report::engine_counts(
                &hits,
                run.queries,
                run.subject_letters,
                &query_code,
                &db_code,
            ),
        }
    } else {
        Vec::new()
    };
    // The display counts are read by the outfmt 0, 6 and 7 reports and the caller's hit
    // records only.
    if outputs.hits.is_some()
        || output_formats
            .iter()
            .any(|format| format.reads_displayed_counts())
    {
        // EXPERIMENT (LOSAT_X_TBXRPTPAR): the identities over chunks of the list on the pool.
        match run.x_pool {
            Some(pool) => crate::algorithm::tblastx::x_rpt_par::set_displayed_identities(
                pool,
                &mut hits,
                run.queries,
                run.subject_letters,
                &query_code,
                &db_code,
            ),
            None => report::set_displayed_identities(
                &mut hits,
                run.queries,
                run.subject_letters,
                &query_code,
                &db_code,
            ),
        }
    }
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The caller receives the same final result that every formatter prints, with the
    // displayed rows of outfmt 0.
    let pairwise = (outputs.hits.is_some()
        || output_formats.contains(&report::TblastxOutputFormat::Pairwise))
    // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows as an indexed map on the pool.
    .then(|| match run.x_pool {
        Some(pool) => crate::algorithm::tblastx::x_rpt_par::pairwise_hits(
            pool,
            &hits,
            &report::PairwiseRows::new(
                run.queries,
                run.subjects,
                run.subject_letters,
                run.query_stats,
                &query_code,
                &db_code,
                &run.ranges.input.placements,
                run.ranges.subject_placements,
            ),
        ),
        None => report::pairwise_hits(
            &hits,
            run.queries,
            run.subjects,
            run.subject_letters,
            run.query_stats,
            &query_code,
            &db_code,
            &run.ranges.input.placements,
            run.ranges.subject_placements,
        ),
    });
    // The tabular rows print record coordinates (`report::shift_to_records`); the outfmt 0
    // rows above were read from the searched letters first.
    report::shift_to_records(
        &mut hits,
        &run.ranges.input.placements,
        run.ranges.subject_placements,
    );
    if let (Some(hits_sink), Some(pairwise)) = (outputs.hits.as_mut(), pairwise.as_ref()) {
        hits_sink(pairwise);
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output (`tabular_subject_ids`).
    let subject_ids = if output_formats.iter().any(|format| format.is_tabular()) {
        match &run.database {
            Some(database) => database.reported.tabular_subject_ids()?,
            None => tabular_subject_ids(&hits, run.subjects)?,
        }
    } else {
        Vec::new()
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    //
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // ```c++
    // 	alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    // ```
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index, as BLASTN's.
    let subject_titles: Vec<Arc<[u8]>> =
        if output_formats.contains(&report::TblastxOutputFormat::Pairwise) {
            run.subjects
                .iter()
                .map(|subject| Arc::from(subject.title.as_slice()))
                .collect()
        } else {
            Vec::new()
        };
    let database = match &run.database {
        Some(database) => database.name.to_string(),
        None => tblastx_database_name(args),
    };
    // The `Query=` and `# Query:` text of each query: its title bytes.
    let query_titles: Vec<&[u8]> = run.queries.iter().map(|query| query.title()).collect();
    let unsearched: Vec<bool> = run.searched.iter().map(|&searched| !searched).collect();
    // The warnings of each query are written before its report, once, with the first format
    // (`QueryWarnings`); the taxonomy database's messages of XML2 and JSON go to the
    // diagnostics at every format.
    let mut warning_lines = Some(run.warnings);
    let observer = &mut outputs.observer;
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(&output_formats).enumerate()
    {
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        let mut warnings = (warning_lines.is_some()
            || matches!(output_format, report::TblastxOutputFormat::Xml2(_)))
        .then(|| crate::report::query_warnings::QueryWarnings {
            before: warning_lines.take().unwrap_or(&[]),
            sink: &mut *outputs.diagnostics,
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let prolog = !prolog_before_search(&format.sink);
        let mut writer = format.sink.open()?;
        match output_format {
            report::TblastxOutputFormat::Pairwise => {
                let pairwise = pairwise
                    .as_deref()
                    .expect("pairwise hits are built for outfmt 0");
                let queries: Vec<crate::report::pairwise::TblastxPairwiseQuery> = run
                    .queries
                    .iter()
                    .zip(run.query_stats)
                    .zip(&query_titles)
                    .enumerate()
                    .map(|(q_idx, ((query, stats), title))| {
                        crate::report::pairwise::TblastxPairwiseQuery {
                            query_name: title.to_vec(),
                            // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:742-744
                            // ```c++
                            //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
                            //             out << "\nLength=";
                            //             out << cbs.GetInst().GetLength() <<"\n";
                            // ```
                            // The query record's length, also with `-query_loc`.
                            query_length: run
                                .ranges
                                .input
                                .placements
                                .length(q_idx, query.seq().len()),
                            karlin: stats.karlin,
                            effective_search_space: stats.eff_searchsp,
                        }
                    })
                    .collect();
                let (num_descriptions, num_alignments) = (
                    args.report_counts.num_descriptions,
                    args.report_counts.num_alignments,
                );
                let pairwise_report = crate::report::pairwise::TblastxPairwiseReport {
                    version: report::NCBI_TBLASTX_VERSION.to_string(),
                    database_name: database.clone(),
                    // `-db`: the `Database:` lines are the database's (`database` below).
                    database_num_sequences: match &run.database {
                        Some(database) => database.pairwise.num_sequences as usize,
                        None => run.subjects.len(),
                    },
                    database_total_letters: match &run.database {
                        Some(database) => database.pairwise.total_letters as usize,
                        None => run.subjects.iter().map(|r| r.seq().len()).sum(),
                    },
                    word_threshold: args.threshold,
                    window_size: args.window_size as usize,
                    num_descriptions,
                    num_alignments,
                    unsearched: unsearched.clone(),
                    // A search that stops at a query batch (`Empty CBlastQueryVector` with
                    // `-query_loc`) writes no epilog (tblastx_app.cpp:209, `PrintEpilog`).
                    epilog: run.epilog,
                    prolog,
                };
                crate::report::pairwise::write_tblastx_pairwise_report(
                    pairwise,
                    &mut writer,
                    &crate::report::pairwise::PairwiseConfig {
                        program: "tblastx".to_string(),
                        ..crate::report::pairwise::PairwiseConfig::default()
                    },
                    &queries,
                    &subject_titles,
                    &pairwise_report,
                    run.database.map(|database| database.pairwise),
                    probe.as_mut(),
                    warnings.as_mut(),
                    run.x_pairwise_pool,
                )?;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText || m_FormatType == CFormattingArgs::eJsonSeqalign) {
            //         if (results.HasAlignments()) {
            //     	    CRef<CSeq_align_set> aln_set (new CSeq_align_set);
            //             CBlastFormatUtil::PruneSeqalign(*(results.GetSeqAlign()), *aln_set, m_HitlistSize);
            //             if(m_FormatType == CFormattingArgs::eAsnText)
            //             	m_Outfile << MSerial_AsnText << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //             else
            //             	m_Outfile << MSerial_Json << *x_WrapAlignmentInSeqAnnot(aln_set, db_title);
            //         }
            // ```
            // The Seq-annot of each query with alignments.
            report::TblastxOutputFormat::Seqalign(serial) => {
                crate::algorithm::tblastx::seqalign::write_tblastx_seqalign(
                    serial,
                    &mut writer,
                    &hits,
                    &engine_counts,
                    run.queries,
                    &run.ranges.input.placements,
                    &match &run.database {
                        Some(database) => {
                            crate::algorithm::tblastx::seqalign::SeqalignSubjects::Database {
                                ordinals: &database.reported.ordinals,
                                title: &database.pairwise.definition,
                            }
                        }
                        None => crate::algorithm::tblastx::seqalign::SeqalignSubjects::Records {
                            records: run.subjects,
                            placements: run.ranges.subject_placements,
                        },
                    },
                    run.subject_letters,
                    warnings.as_mut(),
                )?
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-744
            // ```c++
            //     } else if (m_FormatType == CFormattingArgs::eXml) {
            //         ...
            //         	BlastXML_FormatReport(xml_output, &report_data, &m_Outfile,
            //                                   m_BlastXMLIncremental.GetPointer());
            // ```
            // The document of every query; the hits are the final list after `shift_to_records`
            // (record coordinates), the rows are read from the searched letters.
            report::TblastxOutputFormat::Xml => {
                crate::algorithm::tblastx::blast_xml::write_tblastx_xml(
                    &mut crate::report::blast_xml::BlastXmlWriter::new(0),
                    false,
                    &mut writer,
                    &hits,
                    run,
                    warnings.as_mut(),
                )?
            }
            // The XML2 and JSON documents (outfmt 16 and 15, and 14 and 13 as side files next to
            // `-out` with this output as the master file): the same records, written by
            // `Xml2Writer`.
            report::TblastxOutputFormat::Xml2(kind) => {
                crate::algorithm::tblastx::blast_xml::write_tblastx_xml(
                    &mut crate::report::blast_xml2::Xml2Writer::for_run(
                        kind,
                        crate::algorithm::tblastx::blast_xml::xml2_run(run),
                        run.args.out.as_deref(),
                    ),
                    true,
                    &mut writer,
                    &hits,
                    run,
                    warnings.as_mut(),
                )?
            }
            report::TblastxOutputFormat::Tabular
            | report::TblastxOutputFormat::TabularWithComments => match run.x_pool {
                // EXPERIMENT (LOSAT_X_TBXRPTPAR): the rows formatted on the pool; with the format
                // observer (which flushes around every row) the serial writer.
                Some(pool) if probe.is_none() => {
                    crate::algorithm::tblastx::x_rpt_par::write_tabular(
                        pool,
                        &hits,
                        &mut writer,
                        output_format == report::TblastxOutputFormat::TabularWithComments,
                        &query_titles,
                        &database,
                        run.query_ids,
                        &subject_ids,
                        &unsearched,
                        run.epilog,
                        warnings.as_mut(),
                    )?
                }
                _ => report::write_tabular(
                    &hits,
                    &mut writer,
                    output_format == report::TblastxOutputFormat::TabularWithComments,
                    &query_titles,
                    &database,
                    run.query_ids,
                    &subject_ids,
                    &unsearched,
                    run.epilog,
                    probe.as_mut(),
                    warnings.as_mut(),
                )?,
            },
        }
        writer.flush()?;
    }
    // Without an output format the warnings are written all at once.
    if let Some(lines) = warning_lines {
        for query_warnings in lines {
            outputs.diagnostics.write_all(query_warnings)?;
        }
    }
    Ok(())
}

/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`), as BLASTN's.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
/// ```c++
///         } catch (const CException&) {
///             list<CRef<CSeq_id> > subject_ids;
///             CRef<CSeq_id> id(new CSeq_id());
///             id->Assign(align.GetSeq_id(1));
///             subject_ids.push_back(id);
///             SetSubjectId(subject_ids);
///             bioseqs_found = false;
///         }
/// ```
/// A subject with a row whose `GenerateDefline` throws (a `Subject_` first word and a title
/// in no encoding that `CUtf8::GuessEncoding` recognises) makes NCBI 2.17.0 write part of the
/// row and stop with a `CCoreException` that names its build's source files (exit 255;
/// `AUTHORITY.md` §J-6, §K-1, RP-20). No result can be defined for it, so LOSAT rejects such
/// a subject explicitly before it writes any format; a subject without rows is never
/// named.
fn tabular_subject_ids(hits: &[TblastxHsp], subjects: &[FastaRecord]) -> Result<Vec<Arc<[u8]>>> {
    let ids: Vec<Option<Arc<[u8]>>> = subjects
        .iter()
        .map(|subject| {
            crate::report::defline::tabular_subject_id(
                subject.local_id.as_bytes(),
                &subject.title,
                false,
            )
            .ok()
            .map(Arc::from)
        })
        .collect();
    for hsp in hits {
        let s_idx = hsp.hit.s_idx as usize;
        if ids.get(s_idx).is_some_and(Option::is_none) {
            anyhow::bail!(
                "subject record {} ({}) has a title whose first word starts with 'Subject_' and whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises; NCBI BLAST+ writes part of a tabular row for such a subject and stops with an exception that names its build's source files (exit 255), which is not supported by LOSAT's TBLASTX",
                s_idx + 1,
                String::from_utf8_lossy(subjects[s_idx].shown_id())
            );
        }
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| Arc::from(&b""[..])))
        .collect())
}

/// The results of the query batches searched so far, by query.
#[derive(Default)]
pub(super) struct SearchedBatches {
    pub(super) hits: Vec<TblastxHsp>,
    pub(super) query_stats: Vec<TblastxQueryStats>,
    pub(super) searched: Vec<bool>,
}

/// What the reports of the batches searched so far need besides their results.
pub(super) struct BatchReports<'a> {
    pub(super) args: &'a TblastxArgs,
    pub(super) query_runs: &'a QueryRuns<'a>,
    pub(super) query_ids: &'a [Arc<[u8]>],
    pub(super) subjects: &'a ShownSubjects<'a>,
    pub(super) numbers: &'a SubjectSetNumbers,
    /// `-db`: the stored titles of the batches' reported sequences (`run_in_pool`).
    pub(super) titles: std::cell::RefCell<DatabaseTitles>,
    pub(super) v1_titles: bool,
    /// EXPERIMENT (LOSAT_X_TBXRPTPAR): the search's pool, for the report stage
    /// (`TblastxReportRun::x_pool`).
    pub(super) x_search_pool: &'a crate::utils::threading::SearchPool<'a>,
}

impl BatchReports<'_> {
    /// Writes the reports of the queries of the batches searched (`done`), with NCBI's
    /// epilog (`PrintEpilog`) when `epilog`. `warnings` are by query (at least those of
    /// `done`).
    pub(super) fn write(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
        epilog: bool,
    ) -> Result<()> {
        // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); written
        // after the batches are searched. No NCBI counterpart.
        let t_output = crate::utils::xtiming::clock();
        let reported = done.searched.len();
        let queries = &self.query_runs.records[..reported];
        let mut hits = report_order(self.args, queries, done.hits);
        // `-db`: the reported sequences by rank (their records hold the stored titles, their
        // letters are read from the stored codes).
        let reported_subjects;
        let (subjects, subject_letters, database) = match self.subjects {
            ShownSubjects::Fasta(records) => {
                (*records, report::SubjectLetters::Records(records), None)
            }
            ShownSubjects::Database {
                database,
                pairwise,
                stores,
            } => {
                reported_subjects = crate::algorithm::tblastx::database::reported_subjects(
                    &self.titles.borrow(),
                    &mut hits,
                    stores,
                )?;
                (
                    reported_subjects.records.as_slice(),
                    report::SubjectLetters::Stored {
                        stores,
                        ordinals: &reported_subjects.ordinals,
                    },
                    Some(ReportDatabase {
                        name: database.name(),
                        pairwise,
                        reported: &reported_subjects,
                    }),
                )
            }
        };
        // `-db`, XML2 and JSON: the taxid of each reported sequence, and what the taxonomy
        // database says about it (the refusals of both come before any byte).
        let xml_taxa = match (&self.subjects, &database) {
            (ShownSubjects::Database { database, .. }, Some(shown))
                if outputs.formats.iter().any(|format| {
                    matches!(
                        report::output_format(format.outfmt),
                        report::TblastxOutputFormat::Xml2(_)
                    )
                }) =>
            {
                database.read_taxids(&shown.reported.ordinals)?
            }
            _ => Default::default(),
        };
        let written = write_tblastx_outputs(
            hits,
            outputs,
            &TblastxReportRun {
                args: self.args,
                queries,
                subjects,
                subject_letters,
                database,
                query_ids: &self.query_ids[..reported],
                query_stats: &done.query_stats,
                searched: &done.searched,
                warnings: &warnings[..reported],
                ranges: self.query_runs,
                epilog,
                numbers: self.numbers,
                xml_taxa: &xml_taxa,
                v1_titles: self.v1_titles,
                x_pool: crate::algorithm::tblastx::x_rpt_par::pool(self.x_search_pool),
                x_pairwise_pool: crate::report::x_pairwise_par::pool(self.x_search_pool),
            },
        );
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
        written
    }

    /// Writes what NCBI has written when a query batch fails, before its error: the reports
    /// of the batches before (`done`; without the epilog), or, for the first batch, the
    /// outfmt 0 prolog of the sinks that write it with the reports (the others wrote it
    /// before the search, `write_pairwise_prologs`).
    pub(super) fn write_before_failed_batch(
        &self,
        outputs: &mut ReportOutputs<'_>,
        done: SearchedBatches,
        warnings: &[Vec<u8>],
    ) -> Result<()> {
        if done.searched.is_empty() {
            write_pairwise_prologs(outputs, self.args, self.subjects, false)?;
            // The XML2 and JSON documents start with their prolog too.
            for format in outputs.formats.iter_mut() {
                if let report::TblastxOutputFormat::Xml2(kind) =
                    report::output_format(format.outfmt)
                {
                    let mut writer = format.sink.open()?;
                    writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
                    writer.flush()?;
                }
            }
            Ok(())
        } else {
            self.write(outputs, done, warnings, false)
        }
    }
}
