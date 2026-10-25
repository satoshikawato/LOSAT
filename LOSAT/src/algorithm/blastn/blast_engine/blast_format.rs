//! BLASTN's reports: the pairwise, tabular, XML, SAM and archive output of a query batch.
//!
//! This file is a child module of `blast_engine/run.rs` (declared there with `#[path]`, as
//! `database.rs` is), so that it reaches the run's types without changing them.
//!
//! NCBI reference: `c++/src/algo/blast/format/blast_format.cpp:68-2272` (`CBlastFormat`), the
//! Seq-aligns of `c++/src/algo/blast/api/blast_seqalign.cpp:1520-1637`, and the subject ids of
//! `c++/src/objtools/align_format/tabular.cpp:482-869`.

use super::blastn_app::{ReportRanges, NO_DATA_MESSAGE};
use super::{
    context_lengths, database, drop_hsps_of_score_zero, finish_query_hit_list, record_title,
    x_qbatch_par, x_rpt_stream, BlastnTiming, ShownSubjects, SubjectSource,
};
use crate::algorithm::blastn::archive;
use crate::algorithm::blastn::args::{BlastnArgs, QueryStrand};
use crate::algorithm::blastn::coordination::configure_task;
use crate::algorithm::blastn::hsp::{BlastnHitList, BlastnOutputFormat, NCBI_BLASTN_VERSION};
use crate::algorithm::blastn::pairwise::{pairwise_hits, shown_query_masks, DisplayMasks};
use crate::algorithm::blastn::scoring::ContextKarlin;
use crate::algorithm::blastn::seqalign::write_blastn_seqalign;
use crate::algorithm::blastn::x_subject_store::SubjectLetters;
use crate::api::local_blast::{FormatOutput, FormatProbe, OutputSink, ReportOutputs};
use crate::blastinput::fasta_reader::FastaRecord;
use crate::report::pairwise::{
    write_blastn_pairwise_prolog, write_blastn_pairwise_report, BlastnPairwiseQuery,
    BlastnPairwiseReport, PairwiseConfig, PairwiseDatabase,
};
use crate::report::record::{QueryNames, SubjectNames};
use crate::utils::dust::MaskedInterval;
use anyhow::Result;
use std::sync::Arc;

/// The SAM document that a report call adds its queries to (outfmt 17): the stream's, which
/// the last call writes, or the document of a call of its own (the off path), which the call
/// writes at its end. NCBI writes the document when its formatter is destroyed.
struct SamTarget<'a> {
    document: &'a mut crate::report::sam::SamDocument,
    /// Whether the call writes the document after adding its queries.
    finish: bool,
}

/// What the pairwise report, the hit records and the query warnings need besides the
/// final hit list.
struct BlastnReportInputs<'a> {
    queries: &'a [FastaRecord],
    /// The subject records and their letters (theirs, or with LOSAT_X_BNSTORE their packed
    /// subjects').
    subjects: SubjectLetters<'a>,
    /// Where the searched letters lie in the records (`-query_loc`, `-subject_loc`).
    ranges: &'a ReportRanges<'a>,
    /// Per query, plus strand: DUST and, with `-lcase_masking`, the input lowercase.
    query_masks: &'a [Vec<MaskedInterval>],
    /// The DUST-only masks of every query (outfmt 11's `LCaseMask`); empty without it.
    archive_dust: &'a [Vec<MaskedInterval>],
    lcase_masking: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`).
    query_eff_searchsp: &'a [i64],
    /// The first valid context of a query, as its offset in the query's two contexts: 1 for
    /// `-strand minus`, else 0 (`context 2 * q + first_context`).
    first_context: usize,
    megablast: bool,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    /// The task is megablast or blastn (`BlastnPairwiseReport::zero_gap_extension_formula`).
    zero_gap_extension_formula: bool,
    /// The two-hit window of the task, printed in the epilog when not 0.
    window_size: usize,
    /// Per query context (the plus strand of query `q` is context `2 * q`); `None` for
    /// the contexts of an invalid query (not searched by NCBI) and of the strand that
    /// `-strand` does not search.
    context_karlin: &'a [Option<ContextKarlin>],
    /// The numbers of descriptions and alignments of the pairwise report
    /// (`BlastnArgs::report_counts`).
    num_descriptions: usize,
    num_alignments: usize,
    /// `-culling_limit` (0: no culling), read when the hit lists are finished.
    culling_limit: i32,
    db_num_seqs: usize,
    db_len_total: usize,
    /// The `Database:` lines of a database search's outfmt 0 report.
    database: Option<&'a PairwiseDatabase>,
    /// The tabular subject ids of a database search's reported sequences, by rank
    /// (`database::ReportedSubjects::tabular_subject_ids`); `None` for a subject set.
    database_subject_ids: Option<&'a [Arc<[u8]>]>,
    /// The ordinal ID of each reported sequence of a database search, by rank
    /// (`database::ReportedSubjects::ordinals`); `None` for a subject set.
    database_ordinals: Option<&'a [u32]>,
    /// The opened database of a database search, which the XML2 and JSON formats read the taxids
    /// of the reported sequences from; `None` for a subject set.
    database_reader: Option<&'a database::BlastnDatabase>,
    /// Whether the hit lists were finished after each query batch's traceback (a database
    /// search, `run_in_pool`), not by `post_process_hits_and_write`.
    finished_hit_lists: bool,
    /// Whether NCBI's epilog (`PrintEpilog`: the outfmt 0 footer, outfmt 7's last line)
    /// follows the reports; an error in a later batch skips it.
    epilog: bool,
    /// Whether the outfmt 0 report starts with the prolog, which is otherwise written
    /// before the search (`write_pairwise_prologs`).
    prolog: bool,
    /// The queries whose reports earlier calls wrote, which outfmt 7's epilog counts too
    /// (NCBI's `m_QueriesFormatted`): with LOSAT_X_RPTSTREAM the last call writes only the
    /// epilogs, after the reports of every batch (`x_rpt_stream`); 0 otherwise.
    queries_before: usize,
    /// EXPERIMENT (LOSAT_X_RPTSTREAM): the report values of a subject set's subjects, made once
    /// for the reports of all batches (`x_rpt_stream::SubjectCache`); `None`: made here.
    subject_cache: Option<&'a x_rpt_stream::SubjectCache>,
    /// EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool that formats the outfmt 0 rows and text
    /// (`report::x_pairwise_par`); `None`: the serial writer.
    x_pairwise_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// EXPERIMENT (LOSAT_X_BNTABPAR): the pool that formats the outfmt 6 / 7 rows
    /// (`algorithm::blastn::x_tab_par`); `None`: the serial writer.
    x_tab_pool: Option<&'a crate::utils::threading::SearchPool<'a>>,
    /// `-out`: the base of the side files of outfmt 13 and 14 and of the names their master
    /// file lists.
    out: Option<&'a std::path::Path>,
    /// The e-value threshold of the search (BLAST XML's `Parameters_expect`).
    expect: f64,
    /// The DUST options as NCBI's filter string prints them (BLAST XML's `Parameters_filter`);
    /// `None`: no DUST.
    dust: Option<(i32, i32, i32)>,
    /// The database length behind the effective lengths (`-dbsize`, or the total of the
    /// subjects), for the length adjustment that BLAST XML's `Statistics_hsp-len` shows.
    eff_db_length: i64,
    /// The run's BLAST archive (outfmt 11), and the queries of each batch of this call
    /// (`BatchResults::batches`).
    archive: Option<&'a crate::report::archive::ArchiveHandle>,
    batch_spans: &'a [std::ops::Range<usize>],
    /// Whether both strands of each query are searched (`-strand both`).
    both_strands: bool,
    minus_strand: bool,
}

/// The per-query and run-level data of the pairwise report.
// Kept out of line: this runs only for outfmt 0 or hit records, and inlining it into the
// shared post-processing makes every BLASTN run compile it in Wasm hosts.
#[inline(never)]
fn blastn_pairwise_report(
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
) -> Result<(Vec<BlastnPairwiseQuery>, BlastnPairwiseReport)> {
    let context_karlin = report.context_karlin;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:82-104
    // ```c
    //     // find the first valid context corresponding to this query
    // ```
    // The first valid context of a query is its minus context with `-strand minus`.
    let first_context = report.first_context;
    let queries = report
        .queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| BlastnPairwiseQuery {
            query_name: query_titles[q_idx].to_vec(),
            // NCBI reference: c++/src/objtools/align_format/align_format_util.cpp:742-744
            // ```c++
            //         if(cbs.IsSetInst() && cbs.GetInst().CanGetLength()){
            //             out << "\nLength=";
            //             out << cbs.GetInst().GetLength() <<"\n";
            // ```
            // The length of the query's record, also where `-query_loc` searches part of it.
            query_length: report
                .ranges
                .queries
                .placements
                .length(q_idx, query.seq().len()),
            karlin: context_karlin[2 * q_idx + first_context]
                .map(|blocks| (blocks.ungapped, blocks.gapped)),
            // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_results.cpp:82-104
            // ```c
            //     // find the first valid context corresponding to this query
            //     ...
            //     m_SearchSpace = ctx->eff_searchsp;
            // ```
            // An invalid query has no valid context, and its search space stays 0.
            effective_search_space: if context_karlin[2 * q_idx + first_context].is_some() {
                report.query_eff_searchsp[2 * q_idx + first_context]
            } else {
                0
            },
        })
        .collect();
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
    // ```c
    //     m_NumDescriptions = m_DfltNumDescriptions;
    //     m_NumAlignments = m_DfltNumAlignments;
    //     ...
    //     if (args.Exist(kArgMaxTargetSequences) && args[kArgMaxTargetSequences]) {
    //         m_NumDescriptions = args[kArgMaxTargetSequences].AsInteger();
    //         m_NumAlignments = args[kArgMaxTargetSequences].AsInteger();
    // ```
    // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221).
    let (num_descriptions, num_alignments) = (report.num_descriptions, report.num_alignments);
    let pairwise_report = BlastnPairwiseReport {
        version: NCBI_BLASTN_VERSION.to_string(),
        megablast: report.megablast,
        database_name: subject_title.to_string(),
        database_num_sequences: report.db_num_seqs,
        database_total_letters: report.db_len_total,
        database: report.database.cloned(),
        reward: report.reward,
        penalty: report.penalty,
        gap_open: report.gap_open,
        gap_extend: report.gap_extend,
        zero_gap_extension_formula: report.zero_gap_extension_formula,
        window_size: report.window_size,
        num_descriptions,
        num_alignments,
        unsearched,
        epilog: report.epilog,
        prolog: report.prolog,
    };
    Ok((queries, pairwise_report))
}

/// The values of the first query of a batch (`first`, an index of the call) that the archive
/// takes (`archive::FirstContext`): the blocks and search statistics of its first valid
/// context (`CBlastAncillaryData`, blast_results.cpp:72-115; the plus strand's, or the minus
/// strand's with `-strand minus`), and the search space of its context 0 (the option
/// `SplitQuery_SetEffectiveSearchSpace` sets).
fn archive_first_context(
    report: &BlastnReportInputs<'_>,
    query_lengths: &[usize],
    first: usize,
) -> Result<archive::FirstContext> {
    let context0_space = report.query_eff_searchsp[2 * first];
    let Some(blocks) = report.context_karlin[2 * first + report.first_context] else {
        return Ok(archive::FirstContext {
            karlin: None,
            search_space: 0,
            length_adjustment: 0,
            context0_space,
        });
    };
    let search_space = report.query_eff_searchsp[2 * first + report.first_context];
    let (length_adjustment, _) = context_lengths(
        query_lengths[first] as i64,
        report.eff_db_length,
        report.db_num_seqs as i64,
        &blocks.gapped,
    )?;
    Ok(archive::FirstContext {
        karlin: Some(blocks),
        search_space,
        length_adjustment,
        context0_space,
    })
}

/// The `Statistics` of each query of a BLAST XML report (outfmt 5): those of its first valid
/// context, which NCBI's `CBlastAncillaryData` reads, from the gapped Karlin-Altschul block
/// (blastxml_format.cpp:628-650, data4xmlformat.cpp:236-298): the length adjustment, which the
/// search keeps no field for, from `context_lengths` with the search's own inputs, and the
/// effective search space the search kept. A query without a valid context has `XmlStats::NONE`.
fn blastn_xml_stats(
    report: &BlastnReportInputs<'_>,
    query_lengths: &[usize],
) -> Result<Vec<crate::report::blast_xml::XmlStats>> {
    (0..report.queries.len())
        .map(|q_idx| {
            let Some(blocks) = report.context_karlin[2 * q_idx + report.first_context] else {
                return Ok(crate::report::blast_xml::XmlStats::NONE);
            };
            let eff_space = report.query_eff_searchsp[2 * q_idx + report.first_context];
            let (length_adjustment, searched_space) = context_lengths(
                query_lengths[q_idx] as i64,
                report.eff_db_length,
                report.db_num_seqs as i64,
                &blocks.gapped,
            )?;
            debug_assert_eq!(
                searched_space, eff_space,
                "the search's own effective space"
            );
            Ok(crate::report::blast_xml::XmlStats {
                hsp_len: length_adjustment as i32,
                eff_space,
                kappa: blocks.gapped.k,
                lambda: blocks.gapped.lambda,
                entropy: blocks.gapped.h,
            })
        })
        .collect()
}

/// Moves the reported coordinates of every HSP from the searched letters to the records:
/// by the start of the query's interval (`RemapToQueryLoc`) and of the subject's interval
/// (`s_RemapToSubjectLoc`), on both strands. Without ranges nothing moves.
///
/// NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
/// ```c++
/// void RemapToQueryLoc(CRef<CSeq_align> sar, const CSeq_loc & query)
/// {
///     _ASSERT(sar);
///     const int query_row = 0;
///
///     TSeqPos q_shift = 0;
///
///     if (query.IsInt()) {
///         q_shift = query.GetInt().GetFrom();
///     }
///
///     if (q_shift > 0) {
///         sar->OffsetRow(query_row, q_shift);
///     }
/// }
/// ```
/// NCBI reference: c++/src/objmgr/util/seq_align_util.cpp:71-78
/// ```c++
///     // Create source seq-loc
///     CSeq_loc src_loc(*id, 0, GetLength(loc, scope) - 1);
///     ENa_strand strand = loc.GetStrand();
///     if (strand != eNa_strand_unknown) {
///         src_loc.SetStrand(strand);
///     }
///     CSeq_loc_Mapper mapper(src_loc, loc, scope);
///     return mapper.Map(align, row);
/// ```
/// A subject interval has the strand both, so its mapping is the shift by its start.
fn shift_to_records(hit_lists: &mut [Option<BlastnHitList>], ranges: &ReportRanges<'_>) {
    let (queries, subjects) = (&ranges.queries.placements, ranges.subjects.placements);
    if !queries.ranged() && !subjects.ranged() {
        return;
    }
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            for hsp in &mut hsp_list.hsps {
                let q_shift = queries.offset(hsp.q_idx as usize);
                let s_shift = subjects.offset(hsp.s_idx as usize);
                hsp.q_start += q_shift;
                hsp.q_end += q_shift;
                hsp.s_start += s_shift;
                hsp.s_end += s_shift;
            }
        }
    }
}

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) of every subject, by
/// subject index (`report::defline::tabular_subject_id`: the title's first word, or the
/// local ID, and for `lcl|Subject_...` the first word of `GenerateDefline`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-841
/// ```c++
///     bool setSubjectId = (x_IsFieldRequested(eSubjectSeqId) ||
///     		 	 	 	 x_IsFieldRequested(eSubjectGi) ||
///     		 	 	 	 x_IsFieldRequested(eSubjectAccession) ||
///     		 	 	 	 x_IsFieldRequested(eSubjAccessionVersion));
///
///     if(setSubjectIds || setSubjectTaxInfo || setSubjectTaxInfoAll || setSubjectTitle ||
///        x_IsFieldRequested(eSubjectStrand) || setSubjectId)
///     {
///         try {
///        		const CBioseq_Handle& subject_bh =
///                 scope.GetBioseqHandle(align.GetSeq_id(1));
///             if(setSubjectId) {
///             	SetSubjectId(subject_bh);
///             }
/// ```
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:862-869
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
/// a subject explicitly before it writes any format, pending the PD-LOSAT-NCBI-DEFECTS
/// question; a subject without rows is never named.
fn tabular_subject_ids(
    hit_lists: &[Option<BlastnHitList>],
    subjects: &[FastaRecord],
) -> Result<Vec<Arc<[u8]>>> {
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
    for hit_list in hit_lists.iter().flatten() {
        for hsp_list in &hit_list.hsplist_array {
            let Some(hsp) = hsp_list.hsps.first() else {
                continue;
            };
            let s_idx = hsp.s_idx as usize;
            if ids.get(s_idx).is_some_and(Option::is_none) {
                anyhow::bail!(
                    "subject record {} ({}) has a title whose first word starts with 'Subject_' and whose non-UTF-8 bytes are in no encoding that NCBI BLAST+ recognises; NCBI BLAST+ writes part of a tabular row for such a subject and stops with an exception that names its build's source files (exit 255), which is not supported by LOSAT's BLASTN",
                    s_idx + 1,
                    String::from_utf8_lossy(subjects[s_idx].shown_id())
                );
            }
        }
    }
    Ok(ids
        .into_iter()
        .map(|id| id.unwrap_or_else(|| Arc::from(&b""[..])))
        .collect())
}

fn post_process_hits_and_write(
    mut hit_lists: Vec<Option<BlastnHitList>>,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    query_cov_hsp_perc: f64,
    subject_besthit: bool,
    query_lengths: &[usize],
    outputs: &mut ReportOutputs<'_>,
    verbose: bool,
    query_ids: &[Arc<[u8]>],
    output_formats: &[BlastnOutputFormat],
    query_titles: &[Arc<[u8]>],
    subject_title: &str,
    timing: Option<&BlastnTiming>,
    report: &BlastnReportInputs<'_>,
    unsearched: Vec<bool>,
    query_warnings: &[Vec<u8>],
    mut sam: Vec<Option<SamTarget<'_>>>,
) -> Result<()> {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:183-187
    // ```c
    // typedef struct BlastHSPResults {
    //    Int4 num_queries;
    //    BlastHitList** hitlist_array;
    // } BlastHSPResults;
    // ```
    if verbose {
        let total_hits: usize = hit_lists
            .iter()
            .filter_map(|list| list.as_ref())
            .map(|list| {
                list.hsplist_array
                    .iter()
                    .map(|l| l.hsps.len())
                    .sum::<usize>()
            })
            .sum();
        eprintln!(
            "[INFO] Received {} raw hits total, starting post-processing...",
            total_hits
        );
    }
    let chain_start = std::time::Instant::now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:836-867
    // ```c
    // for (subject_index = 0; subject_index < hit_list->hsplist_count; ++subject_index) {
    //     BlastHSPList * hsp_list = hit_list->hsplist_array[subject_index];
    //     if (hit_options->max_hsps_per_subject) {
    //         Blast_TrimHSPListByMaxHsps(hsp_list, hit_options);
    //     }
    //     if ((hit_options->hsp_filt_opt != NULL) &&
    //         (hit_options->hsp_filt_opt->subject_besthit_opts != NULL)) {
    //         Blast_HSPListSubjectBestHit(program_number,
    //             hit_options->hsp_filt_opt->subject_besthit_opts,
    //             query_info, hsp_list);
    //     }
    // }
    // ```
    // A database search has finished the hit lists of each query batch after its traceback
    // (`run_in_pool`: NCBI reads the definition lines of the batch's final hit lists between
    // the two steps), so they are not finished again.
    if !report.finished_hit_lists {
        for (q_idx, hit_list_opt) in hit_lists.iter_mut().enumerate() {
            let hit_list = match hit_list_opt {
                Some(value) => value,
                None => continue,
            };
            finish_query_hit_list(
                hit_list,
                query_lengths.get(q_idx).copied().unwrap_or_default(),
                hitlist_size,
                max_hsps_per_subject,
                query_cov_hsp_perc,
                subject_besthit,
                report.culling_limit,
            );
            drop_hsps_of_score_zero(hit_list);
        }
    }

    let total_hits: usize = hit_lists
        .iter()
        .filter_map(|list| list.as_ref())
        .map(|list| {
            list.hsplist_array
                .iter()
                .map(|l| l.hsps.len())
                .sum::<usize>()
        })
        .sum();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:877-892
    // ```c
    // static void s_BlastPruneExtraHits(BlastHSPResults* results, Int4 hitlist_size)
    // ```
    if let Some(timing) = timing {
        let post_elapsed = chain_start.elapsed();
        let post_elapsed_ns = post_elapsed.as_nanos() as u64;
        BlastnTiming::record_ns(&timing.traceback_ns, post_elapsed_ns);
        BlastnTiming::record_ns(&timing.traceback_postprocess_ns, post_elapsed_ns);
        timing
            .traceback_calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if verbose {
            eprintln!(
                "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
                post_elapsed.as_secs_f64(),
                total_hits
            );
        }
    } else if verbose {
        eprintln!(
            "[INFO] Post-processing done in {:.2}s, {} hits after filtering, writing output...",
            chain_start.elapsed().as_secs_f64(),
            total_hits
        );
    }
    let write_start = std::time::Instant::now();

    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_hits.c:3077-3106
    // ```c
    // static int s_EvalueCompareHSPLists(const void* v1, const void* v2) { ... }
    // ```
    // NCBI reference: ncbi-blast/c++/src/app/blast/blast_formatter.cpp:429-467
    // ```c
    // CRef<CSearchResultSet> results = m_RmtBlast->GetResultSet();
    // formatter.PrintProlog();
    // ...
    // ITERATE(CSearchResultSet, result, *results) {
    //     ...
    //         formatter.PrintOneResultSet(**result, queries);
    //     ...
    // }
    // ```
    // Each requested format prints the same final result without searching again.
    // `unsearched` marks the queries of the batches without a valid context, which NCBI
    // does not search (outfmt 0 and 7 show them without results).
    // Rendered hits are needed by outfmt 0 and by the caller's hit records; both get the
    // same final hit list, built once.
    // NCBI reference: c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
    // ```c++
    // 	if (seqinfo_src->CanReturnPartialSequence() == true)
    // 	{
    //         	CConstRef<CSeq_loc> subj_loc = seqinfo_src->GetSeqLoc(kOid);
    //         	NON_CONST_ITERATE(vector<CRef<CSeq_align > >, iter, hit_align) {
    //              	   RemapToQueryLoc(*iter, query_loc);
    //              	   if ( !is_ooframe )
    //                    	s_RemapToSubjectLoc(*iter, *subj_loc);
    // ```
    // The reports print record coordinates: the HSPs move by the start of the query's and
    // the subject's interval (`shift_to_records`); their internal offsets stay those of
    // the searched letters, from which the shown rows are read.
    shift_to_records(&mut hit_lists, report.ranges);
    let pairwise_hits =
        if outputs.hits.is_some() || output_formats.contains(&BlastnOutputFormat::Pairwise) {
            // The records as read with `-subject_loc` (never packed), else the searched ones.
            let subject_records = match report.ranges.subjects.records {
                Some(records) => SubjectLetters::new(records, None),
                None => report.subjects,
            };
            // EXPERIMENT (LOSAT_X_RPTSTREAM): the same masks, made once for every batch.
            let subject_masks = report.lcase_masking.then(|| {
                match report
                    .subject_cache
                    .and_then(|cache| cache.masks.as_deref())
                {
                    Some(masks) => std::borrow::Cow::Borrowed(masks),
                    None => std::borrow::Cow::Owned(
                        (0..subject_records.records().len())
                            .map(|index| subject_records.lowercase_masks(index))
                            .collect::<Vec<_>>(),
                    ),
                }
            });
            let mut hits = pairwise_hits(
                &hit_lists,
                report.queries,
                report.subjects,
                &DisplayMasks {
                    query: report.query_masks,
                    query_offsets: &report.ranges.queries.placements,
                    subject: subject_masks.as_deref(),
                },
                report.x_pairwise_pool,
            )?;
            // NCBI reference: c++/src/objtools/align_format/showalign.cpp:2479
            // ```c++
            //     out<<kLengthString<<bsp_handle.GetBioseqLength()<<"\n";
            // ```
            // The length of the subject's record, also where `-subject_loc` searches part
            // of it.
            for hit in &mut hits {
                let s_idx = hit.hit.s_idx as usize;
                hit.subject_length = hit
                    .subject_length
                    .map(|searched| report.ranges.subjects.placements.length(s_idx, searched));
            }
            Some(hits)
        } else {
            None
        };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:1411
    // ```c
    // CBlastFormat::PrintOneResultSet(const blast::CSearchResults& results,
    // ```
    // The caller receives the same final result that every formatter prints.
    if let (Some(hits_sink), Some(hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref()) {
        hits_sink(hits);
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output.
    let pairwise_data = if output_formats.contains(&BlastnOutputFormat::Pairwise) {
        Some(blastn_pairwise_report(
            report,
            unsearched.clone(),
            query_titles,
            subject_title,
        )?)
    } else {
        None
    };
    // EXPERIMENT (LOSAT_X_RPTSTREAM): a subject set's ids made once for every batch; the
    // switch acts only when no subject's id is refused (`x_rpt_stream::SubjectCache`), so
    // they are the ids made here.
    let subject_ids: std::borrow::Cow<'_, [Arc<[u8]>]> =
        if output_formats.iter().any(|format| format.is_tabular()) {
            match (report.database_subject_ids, report.subject_cache) {
                (Some(ids), _) => std::borrow::Cow::Borrowed(ids),
                (None, Some(cache)) => std::borrow::Cow::Borrowed(cache.ids.as_slice()),
                (None, None) => std::borrow::Cow::Owned(tabular_subject_ids(
                    &hit_lists,
                    report.subjects.records(),
                )?),
            }
        } else {
            std::borrow::Cow::Owned(Vec::new())
        };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // ```c++
    //     //get defline
    //     sdl->defline = CDeflineGenerator().GenerateDefline(m_ScopeRef->GetBioseqHandle(*(sdl->id)), sequence::CDeflineGenerator::fLeavePrefixSuffix);
    // ```
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // ```c++
    // 	alnDispParams->title = CDeflineGenerator().GenerateDefline(bsp_handle);
    // ```
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index.
    // EXPERIMENT (LOSAT_X_RPTSTREAM): the same titles, made once for every batch.
    let subject_titles: std::borrow::Cow<'_, [Arc<[u8]>]> =
        if output_formats.contains(&BlastnOutputFormat::Pairwise) {
            match report
                .subject_cache
                .and_then(|cache| cache.titles.as_deref())
            {
                Some(titles) => std::borrow::Cow::Borrowed(titles),
                None => report.subjects.records().iter().map(record_title).collect(),
            }
        } else {
            std::borrow::Cow::Owned(Vec::new())
        };
    let observer = &mut outputs.observer;
    // The warnings are written once, with the first format (`QueryWarnings`).
    let mut warning_lines = Some(query_warnings);
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(output_formats).enumerate()
    {
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a search with several formats, `run_local`).
        if let Some(handle) = report.archive {
            handle.set_muted(output_format != BlastnOutputFormat::Archive);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // The taxonomy database's messages of XML2 and JSON go to the diagnostics at every
        // format.
        let mut format_warnings = (warning_lines.is_some()
            || matches!(
                output_format,
                BlastnOutputFormat::Xml2(_) | BlastnOutputFormat::Archive
            ))
        .then(|| crate::report::query_warnings::QueryWarnings {
            before: warning_lines.take().unwrap_or(&[]),
            sink: &mut *outputs.diagnostics,
        });
        // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:68-93
        // ```c
        // CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...)
        //     : m_FormatType(format_type), ..., m_Outfile(outfile),
        // ```
        let mut writer = format.sink.open()?;
        if output_format == BlastnOutputFormat::Pairwise {
            let hits = pairwise_hits
                .as_deref()
                .expect("pairwise hits are built for outfmt 0");
            let (queries, pairwise_report) = pairwise_data
                .as_ref()
                .expect("the pairwise report data is built for outfmt 0");
            write_blastn_pairwise_report(
                hits,
                &mut writer,
                &PairwiseConfig {
                    program: "blastn".to_string(),
                    show_frame: false,
                    ..PairwiseConfig::default()
                },
                queries,
                &subject_titles,
                pairwise_report,
                probe.as_mut(),
                format_warnings.as_mut(),
                report.x_pairwise_pool,
            )?;
        } else if let BlastnOutputFormat::Seqalign(serial) = output_format {
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
            // The Seq-annot of each query with alignments (`seqalign`).
            let subjects = match (report.database_ordinals, report.database) {
                (Some(ordinals), Some(database)) => SubjectNames::Database {
                    ordinals,
                    title: &database.definition,
                },
                _ => SubjectNames::Records(report.subjects.records()),
            };
            write_blastn_seqalign(
                serial,
                &hit_lists,
                &QueryNames {
                    records: report.queries,
                    searched_lengths: query_lengths,
                    placements: &report.ranges.queries.placements,
                },
                subjects,
                report.ranges.subjects.placements,
                &mut writer,
                format_warnings.as_mut(),
            )?;
        } else if output_format == BlastnOutputFormat::Archive {
            // The BLAST archive (outfmt 11): one document per query batch of the call.
            let handle = report
                .archive
                .expect("a BLAST archive request has its handle");
            let subjects = match (report.database_ordinals, report.database) {
                (Some(ordinals), Some(database)) => SubjectNames::Database {
                    ordinals,
                    title: &database.definition,
                },
                _ => SubjectNames::Records(report.subjects.records()),
            };
            let call = archive::ArchiveCall {
                handle,
                queries: report.queries,
                query_lengths,
                query_placements: &report.ranges.queries.placements,
                subject_records: report.database.is_none().then(|| {
                    report
                        .ranges
                        .subjects
                        .records
                        .unwrap_or_else(|| report.subjects.records())
                }),
                subjects,
                subject_placements: report.ranges.subjects.placements,
                hit_lists: &hit_lists,
                query_masks: report.query_masks,
                dust_masks: report.archive_dust,
                unsearched: &unsearched,
                lcase_masking: report.lcase_masking,
                both_strands: report.both_strands,
                minus_strand: report.minus_strand,
                spans: report.batch_spans,
            };
            archive::write_blastn_archives(
                &call,
                format_warnings.as_mut(),
                &|first| archive_first_context(report, query_lengths, first),
                &mut writer,
            )?;
        } else if let BlastnOutputFormat::Sam(_) = output_format {
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:750-755
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eSAM) {
            //         if (results.HasAlignments()) {
            //             m_SamFormatter->Print(*results.GetSeqAlign());
            //         }
            //         return;
            //     }
            // ```
            // The queries' records are added to the run's document, which NCBI writes when its
            // formatter is destroyed (`SamDocument::finish`): after the last call, or in a call
            // that is the run's only one.
            let target = sam[format_index]
                .as_mut()
                .expect("a SAM document is made for outfmt 17");
            if let Some(warnings) = format_warnings.as_mut() {
                for q_idx in 0..hit_lists.len() {
                    warnings.before_query(q_idx, &mut writer)?;
                }
            }
            let names = QueryNames {
                records: report.queries,
                searched_lengths: query_lengths,
                placements: &report.ranges.queries.placements,
            };
            match (
                report.database_ordinals,
                report.database,
                report.database_reader,
            ) {
                (Some(ordinals), Some(database), Some(reader)) => {
                    // A database's sequences are searched whole; the lengths come from its
                    // index (no letters).
                    let lengths = reader.sam_lengths(ordinals)?;
                    let cache = database::SamLetterCache::default();
                    let letters_of = |rank: usize, from: usize, to: usize| {
                        reader.sam_letters(&cache, ordinals[rank], from, to)
                    };
                    crate::algorithm::blastn::sam::add_blastn_sam(
                        target.document,
                        &hit_lists,
                        &names,
                        report.queries_before,
                        SubjectNames::Database {
                            ordinals,
                            title: &database.definition,
                        },
                        report.ranges.subjects.placements,
                        &crate::algorithm::blastn::sam::SubjectLengths {
                            searched: &|rank| lengths[rank],
                            whole: &|rank| lengths[rank],
                        },
                        &mut crate::algorithm::blastn::sam::BlastnSamLetters {
                            queries: report.queries,
                            query_placements: &report.ranges.queries.placements,
                            subjects: crate::algorithm::blastn::sam::SubjectSource::Database(
                                &letters_of,
                            ),
                        },
                    )?;
                }
                _ => {
                    // The records as read with `-subject_loc` (never packed), else the
                    // searched ones: the letters are asked for in record coordinates.
                    let whole_records = match report.ranges.subjects.records {
                        Some(records) => SubjectLetters::new(records, None),
                        None => report.subjects,
                    };
                    crate::algorithm::blastn::sam::add_blastn_sam(
                        target.document,
                        &hit_lists,
                        &names,
                        report.queries_before,
                        SubjectNames::Records(report.subjects.records()),
                        report.ranges.subjects.placements,
                        &crate::algorithm::blastn::sam::SubjectLengths {
                            searched: &|s_idx| report.subjects.len(s_idx),
                            whole: &|s_idx| {
                                report
                                    .ranges
                                    .subjects
                                    .placements
                                    .length(s_idx, report.subjects.len(s_idx))
                            },
                        },
                        &mut crate::algorithm::blastn::sam::BlastnSamLetters {
                            queries: report.queries,
                            query_placements: &report.ranges.queries.placements,
                            subjects: crate::algorithm::blastn::sam::SubjectSource::Records(
                                whole_records,
                            ),
                        },
                    )?;
                }
            }
            if target.finish {
                target.document.finish(&mut writer)?;
            }
        } else if output_format.is_structured_xml() {
            // The XML2 / JSON documents of outfmt 15 and 16 take the path of outfmt 5: the same
            // per-query records, written by another writer.
            let xml2_kind = output_format.xml2_kind();
            // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1431-1441
            // ```c++
            //     if (m_FormatType == CFormattingArgs::eAsnText
            //       || m_FormatType == CFormattingArgs::eAsnBinary
            //       || m_FormatType == CFormattingArgs::eXml
            //       ...
            //     {
            //         x_PrintStructuredReport(results, queries);
            //         return;
            //     }
            // ```
            // The CBlastOutput document of the queries (`blast_xml`).
            let is_database = report.database.is_some();
            let xml_stats = blastn_xml_stats(report, query_lengths)?;
            let xml_messages: Vec<Vec<String>> = report
                .queries
                .iter()
                .enumerate()
                .map(|(q_idx, query)| {
                    if report.context_karlin[2 * q_idx + report.first_context].is_some() {
                        Vec::new()
                    } else if query.seq().is_empty() {
                        vec![NO_DATA_MESSAGE.to_string()]
                    } else {
                        vec![crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()]
                    }
                })
                .collect();
            let version = format!("BLASTN {NCBI_BLASTN_VERSION}");
            let filter =
                crate::report::blast_xml::filter::filter_string(report.dust, None, None, true);
            let subject_records = report.subjects.records();
            let database_length = |rank: usize| report.subjects.len(rank);
            let record_length = |s_idx: usize| {
                report
                    .ranges
                    .subjects
                    .placements
                    .length(s_idx, report.subjects.len(s_idx))
            };
            // The taxids of the reported sequences and the taxonomy database's names, read for
            // the XML2 and JSON formats only.
            let (xml_taxids, xml_taxonomy, xml_taxdb_notes): (
                Vec<i32>,
                Vec<crate::report::blast_xml::HitTaxonomy>,
                crate::algorithm::common::database::TaxdbNotes,
            ) = match (xml2_kind, report.database_ordinals, report.database_reader) {
                (Some(_), Some(ordinals), Some(database)) if is_database => {
                    database.read_taxids(ordinals, &hit_lists)?
                }
                _ => (Vec::new(), Vec::new(), Default::default()),
            };
            // The masked regions of the queries in record coordinates and the queries of the
            // batches that were not searched, which only XML2 shows (outfmt 5 shows neither).
            let (xml_masks, xml_searched): (
                Vec<Vec<crate::report::blast_xml::QueryMask>>,
                Vec<bool>,
            ) = if xml2_kind.is_some() {
                (
                    (0..report.queries.len())
                        .map(|q_idx| {
                            shown_query_masks(
                                &report.query_masks[q_idx],
                                query_lengths[q_idx],
                                report.ranges.queries.placements.offset(q_idx),
                            )
                            .iter()
                            .map(|mask| crate::report::blast_xml::QueryMask {
                                from: mask.start as i32,
                                to: mask.end as i32 - 1,
                                frame: 0,
                            })
                            .collect()
                        })
                        .collect(),
                    unsearched.iter().map(|&skipped| !skipped).collect(),
                )
            } else {
                (Vec::new(), Vec::new())
            };
            let xml_subjects = match report.database_ordinals {
                Some(ordinals) if is_database => crate::report::blast_xml::XmlSubjects::Database {
                    ordinals,
                    titles: crate::report::blast_xml::XmlTitles::Records(subject_records),
                    // Only XML2 shows the taxid.
                    taxids: &xml_taxids,
                    taxonomy: &xml_taxonomy,
                    length: &database_length,
                },
                _ => crate::report::blast_xml::XmlSubjects::Records {
                    records: subject_records,
                    length: &record_length,
                    protein: false,
                },
            };
            let xml = crate::algorithm::blastn::blast_xml::BlastnXml {
                run: crate::report::blast_xml::XmlRun {
                    program: "blastn",
                    version: &version,
                    reference: if report.megablast {
                        crate::report::blast_xml::REFERENCE_MEGABLAST
                    } else {
                        crate::report::blast_xml::REFERENCE_GAPPED_BLAST
                    },
                    db_name: if is_database {
                        subject_title.as_bytes()
                    } else {
                        b""
                    },
                    matrix: None,
                    expect: report.expect,
                    gap_open: report.gap_open,
                    gap_extend: report.gap_extend,
                    reward_penalty: Some((report.reward, report.penalty)),
                    filter: &filter,
                    db_num: if is_database {
                        report.db_num_seqs as i32
                    } else {
                        0
                    },
                    db_len: if is_database {
                        report.db_len_total as i64
                    } else {
                        0
                    },
                    positives: None,
                    mask_char: b'N',
                },
                queries: report.queries,
                searched_lengths: query_lengths,
                stats: xml_stats,
                messages: xml_messages,
                masks: xml_masks,
                searched: xml_searched,
                subjects: xml_subjects,
                taxdb_notes: xml_taxdb_notes,
                letters: report.subjects,
                epilog: report.epilog,
            };
            match xml2_kind {
                None => crate::algorithm::blastn::blast_xml::write_blastn_xml(
                    &xml,
                    &mut crate::report::blast_xml::BlastXmlWriter::new(report.queries_before),
                    &hit_lists,
                    &mut writer,
                    format_warnings.as_mut(),
                )?,
                Some(kind) => {
                    // The gap costs are written for a gapped search, which LOSAT's BLASTN
                    // always is (it does not take `-ungapped`); `cbs` and the genetic codes
                    // are not those of a nucleotide search (blastxml2_format.cpp:531-546).
                    let params = crate::report::blast_xml2::Xml2Run {
                        target: if is_database {
                            crate::report::blast_xml2::Xml2Target::Database
                        } else {
                            crate::report::blast_xml2::Xml2Target::Subjects(subject_records.len())
                        },
                        gapped: true,
                        cbs: 0,
                        query_gencode: None,
                        db_gencode: None,
                    };
                    let mut document =
                        crate::report::blast_xml2::Xml2Writer::for_run(kind, params, report.out);
                    crate::algorithm::blastn::blast_xml::write_blastn_xml(
                        &xml,
                        &mut document,
                        &hit_lists,
                        &mut writer,
                        format_warnings.as_mut(),
                    )?;
                }
            }
        } else {
            crate::algorithm::blastn::x_tab_par::write_output(
                report.x_tab_pool,
                &hit_lists,
                &mut writer,
                query_ids,
                &subject_ids,
                output_format,
                query_titles,
                subject_title,
                &unsearched,
                report
                    .epilog
                    .then_some(report.queries_before + hit_lists.len()),
                probe.as_mut(),
                format_warnings.as_mut(),
            )?;
        }
        writer.flush()?;
    }
    // Without an output format the warnings are written all at once.
    if let Some(handle) = report.archive {
        handle.set_muted(false);
    }
    if let Some(before) = warning_lines {
        for query_warnings in before {
            outputs.diagnostics.write_all(query_warnings)?;
        }
    }

    // NCBI reference: ncbi-blast/c++/src/algo/blast/format/blast_format.cpp:770-782
    // ```c
    // // tabular formatting just prints each alignment in turn
    // if (m_FormatType == CFormattingArgs::eTabular ||
    //     m_FormatType == CFormattingArgs::eTabularWithComments ||
    //     m_FormatType == CFormattingArgs::eCommaSeparatedValues ||
    //     m_FormatType == CFormattingArgs::eCommaSeparatedValuesWithHeader) {
    //   CBlastTabularInfo tabinfo(m_Outfile, m_CustomOutputFormatSpec, kDelim);
    // ```
    if let Some(timing) = timing {
        let write_elapsed = write_start.elapsed();
        timing.format_ns.fetch_add(
            write_elapsed.as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
        timing
            .format_calls
            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        if verbose {
            eprintln!(
                "[INFO] Output written in {:.2}s",
                write_elapsed.as_secs_f64()
            );
        }
    } else if verbose {
        eprintln!(
            "[INFO] Output written in {:.2}s",
            write_start.elapsed().as_secs_f64()
        );
    }
    Ok(())
}

/// Whether the outfmt 0 prolog of `sink` is written before the search, as NCBI writes it
/// (`write_pairwise_prologs`). A file sink is created when its report is written
/// (`OutputSink::open`), so its prolog stays at the start of the report.
fn prolog_before_search(sink: &OutputSink<'_>) -> bool {
    !matches!(sink, OutputSink::File(_))
}

/// Whether the prolog of `format` comes with its reports (or, for a run that stops before any
/// report, from `write_pairwise_prologs`) and not before the search: the outfmt 0 prolog of a
/// file sink (`prolog_before_search`), and the start of the XML2 and JSON documents of outfmt
/// 16 and 15 for every sink. NCBI writes that start before the search too (`PrintProlog`), but
/// LOSAT's refusals of what NCBI could write (a title outside US-ASCII, DECISIONS 20) come
/// after the search and before the first byte, so the start is written with the first report, or
/// by the path of an NCBI error that stops the run before one.
pub(super) fn prolog_with_report(
    format: &FormatOutput<'_>,
    output_format: BlastnOutputFormat,
) -> bool {
    match output_format {
        BlastnOutputFormat::Pairwise => !prolog_before_search(&format.sink),
        BlastnOutputFormat::Xml2(_) => true,
        _ => false,
    }
}

/// Writes and flushes the outfmt 0 prolog of every outfmt 0 sink whose prolog comes before the
/// search (`before_search`) or, for a run that stops before any report, the prolog of the
/// others (`prolog_with_report`: the file sinks of outfmt 0 and every sink of outfmt 15 and 16).
pub(super) fn write_pairwise_prologs(
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    megablast: bool,
    shown: &ShownSubjects,
    before_search: bool,
) -> Result<()> {
    for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
        let writes_prolog = match output_format {
            BlastnOutputFormat::Pairwise => prolog_before_search(&format.sink) == before_search,
            BlastnOutputFormat::Xml2(_) => !before_search,
            _ => false,
        };
        if writes_prolog {
            write_format_prolog(format, output_format, megablast, shown)?;
        }
    }
    Ok(())
}

/// Writes and flushes the prolog of one sink, for the formats that have one (outfmt 0, 15, 16).
pub(super) fn write_format_prolog(
    format: &mut FormatOutput<'_>,
    output_format: BlastnOutputFormat,
    megablast: bool,
    shown: &ShownSubjects,
) -> Result<()> {
    if output_format == BlastnOutputFormat::Pairwise {
        let mut writer = format.sink.open()?;
        write_blastn_pairwise_prolog(
            &mut writer,
            NCBI_BLASTN_VERSION,
            megablast,
            &shown.name,
            shown.num_seqs,
            shown.total_letters,
            shown.database.as_ref(),
        )?;
        writer.flush()?;
    } else if let Some(kind) = output_format.xml2_kind() {
        let mut writer = format.sink.open()?;
        writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
        writer.flush()?;
    }
    Ok(())
}

/// The results of the query batches searched so far, indexed by query (two contexts per
/// query).
pub(super) struct BatchResults<'a> {
    pub(super) hit_lists: Vec<Option<BlastnHitList>>,
    /// The warnings written before each query's report (`QueryWarnings`).
    pub(super) warnings: &'a [Vec<u8>],
    pub(super) lengths: &'a [usize],
    pub(super) unsearched: Vec<bool>,
    pub(super) query_masks: &'a [Vec<MaskedInterval>],
    pub(super) context_karlin: &'a [Option<ContextKarlin>],
    pub(super) query_eff_searchsp: &'a [i64],
    /// The searched queries of each query batch of the call, in order, for the reports of
    /// every batch written at once (the BLAST archive writes one document per batch); empty
    /// when the call is one batch.
    pub(super) batches: &'a [std::ops::Range<usize>],
}

/// Writes what NCBI has written when a query batch fails, before its error: the reports of
/// the batches before (`before`, the queries searched so far; without the epilog), or, for
/// the first batch, the outfmt 0 prolog of the sinks that write it with the reports (the
/// others wrote it before the search, `write_pairwise_prologs`).
#[allow(clippy::too_many_arguments)]
pub(super) fn write_before_failed_batch(
    args: &BlastnArgs,
    before: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
) -> Result<()> {
    x_qbatch_par::debug_assert_slots_joined();
    if before.is_empty() {
        write_pairwise_prologs(
            outputs,
            output_formats,
            args.task == "megablast",
            subjects.shown,
            false,
        )
    } else {
        write_batch_reports(
            args,
            before,
            subjects,
            ranges,
            outputs,
            output_formats,
            results,
            false,
            None,
        )
    }
}

/// The subjects as the reports read them: what the batches searched, how the reports show
/// it, and for a database the stored titles of the sequences that the final hit lists
/// name (`database::read_titles`; empty for a subject set).
#[derive(Clone, Copy)]
pub(super) struct ReportSubjects<'a> {
    pub(super) source: SubjectSource<'a>,
    pub(super) shown: &'a ShownSubjects,
    pub(super) titles: &'a database::DatabaseTitles,
}

/// Writes the reports of `query_records` (the queries of the batches searched), with
/// NCBI's epilog (`PrintEpilog`) when `epilog`; `x_pool` as `write_reports`'.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_batch_reports(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    x_qbatch_par::debug_assert_slots_joined();
    write_reports(
        args,
        query_records,
        subjects,
        ranges,
        outputs,
        output_formats,
        results,
        epilog,
        x_rpt_stream::Written::default(),
        x_pool,
    )
}

/// `write_batch_reports` after its check, with what earlier calls wrote (`written`: nothing,
/// except for LOSAT_X_RPTSTREAM's calls, `x_rpt_stream`, which also run while LOSAT_X_QBATCHPAR
/// slots search the batches ahead).
///
/// EXPERIMENT (LOSAT_X_PAIRWISEPAR): `x_pool` is the search's pool when no thread of a
/// LOSAT_X_QBATCHPAR slot exists (its threads are idle; the slots have their own pools, and the
/// reports already overlap their searches); `None` keeps the serial writer.
#[allow(clippy::too_many_arguments)]
pub(super) fn write_reports(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    written: x_rpt_stream::Written<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    // Diagnostic (LOSAT_TIMING): the report for the common line (utils/xtiming.rs); the time of
    // every call adds up. With streaming (`x_rpt_stream`) the calls come between the query
    // batches, while other batches are searched: that part is inside `search` too. No NCBI
    // counterpart.
    let t_output = crate::utils::xtiming::clock();
    let result = write_reports_untimed(
        args,
        query_records,
        subjects,
        ranges,
        outputs,
        output_formats,
        results,
        epilog,
        written,
        x_pool,
    );
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, t_output);
    result
}

#[allow(clippy::too_many_arguments)]
fn write_reports_untimed(
    args: &BlastnArgs,
    query_records: &[FastaRecord],
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    mut written: x_rpt_stream::Written<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    let megablast = args.task == "megablast";
    let config = configure_task(args);
    let mut hit_lists = results.hit_lists;
    // The SAM documents (outfmt 17), one per SAM format: the stream's, which only its last call
    // writes (the stream acts for one format), else documents of this call's own, which it
    // writes at its end (the off path's single call, and the calls that write the document of
    // the batches before an NCBI error).
    let mut own_sams: Vec<Option<crate::report::sam::SamDocument>> = Vec::new();
    let mut sam_targets: Vec<Option<SamTarget<'_>>> = Vec::new();
    // No SAM format: no target is made (the vectors stay empty and are never indexed).
    let has_sam = output_formats
        .iter()
        .any(|format| format.sam_spec().is_some());
    match written.sam.take() {
        Some(document) => sam_targets.push(Some(SamTarget {
            document,
            finish: epilog,
        })),
        None if !has_sam => {}
        None => {
            own_sams.extend(output_formats.iter().map(|format| {
                format
                    .sam_spec()
                    .map(|spec| x_rpt_stream::sam_document(args, spec))
            }));
            sam_targets.extend(own_sams.iter_mut().map(|document| {
                document.as_mut().map(|document| SamTarget {
                    document,
                    finish: true,
                })
            }));
        }
    }
    // A database: the reported sequences, named by rank from here on, with their stored
    // titles, and their stored letters where a report shows them (outfmt 0 and the hit
    // records; `database::reported_subjects`).
    let reported = match subjects.source {
        SubjectSource::Prepared(_) => None,
        SubjectSource::Database(database) => {
            let letters = outputs.hits.is_some()
                || output_formats.contains(&BlastnOutputFormat::Pairwise)
                || output_formats
                    .iter()
                    .any(|format| format.is_structured_xml());
            Some(database::reported_subjects(
                database,
                subjects.titles,
                &mut hit_lists,
                letters,
            )?)
        }
    };
    let database_subject_ids = match &reported {
        Some(reported) if output_formats.iter().any(|format| format.is_tabular()) => {
            Some(reported.tabular_subject_ids()?)
        }
        _ => None,
    };
    let subject_letters = match (&reported, subjects.source) {
        (Some(reported), _) => SubjectLetters::new(&reported.records, None),
        (None, SubjectSource::Prepared(prepared)) => prepared.letters(),
        (None, SubjectSource::Database(_)) => SubjectLetters::new(&[], None),
    };
    // The archive's `LCaseMask` holds the DUST intervals of every query, valid or not, whether
    // or not its batch was searched.
    let archive_dust = if output_formats.contains(&BlastnOutputFormat::Archive) {
        let serial = crate::utils::threading::SearchPool::x_serial();
        crate::algorithm::blastn::coordination::apply_dust_masking(
            args,
            query_records,
            x_pool.unwrap_or(&serial),
        )
    } else {
        Vec::new()
    };
    let report_inputs = BlastnReportInputs {
        queries: query_records,
        subjects: subject_letters,
        ranges,
        query_masks: results.query_masks,
        archive_dust: &archive_dust,
        lcase_masking: args.lcase_masking,
        query_eff_searchsp: results.query_eff_searchsp,
        first_context: QueryStrand::from_option(&args.strand).first_context(),
        megablast,
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:2272
        // ```c
        //         if ((m_Program == "megablast" || m_Program == "blastn") && options.GetGapExtensionCost() == 0)
        // ```
        zero_gap_extension_formula: matches!(args.task.as_str(), "megablast" | "blastn"),
        window_size: config.window_size,
        context_karlin: results.context_karlin,
        num_descriptions: args.report_counts.num_descriptions,
        num_alignments: args.report_counts.num_alignments,
        culling_limit: args.culling_limit_value(),
        db_num_seqs: subjects.shown.num_seqs,
        db_len_total: subjects.shown.total_letters,
        database: subjects.shown.database.as_ref(),
        database_subject_ids: database_subject_ids.as_deref(),
        database_ordinals: reported
            .as_ref()
            .map(|reported| reported.ordinals.as_slice()),
        database_reader: match subjects.source {
            SubjectSource::Database(database) => Some(database),
            SubjectSource::Prepared(_) => None,
        },
        finished_hit_lists: reported.is_some(),
        epilog,
        prolog: outputs
            .formats
            .iter()
            .zip(output_formats)
            .any(|(format, &output_format)| {
                output_format == BlastnOutputFormat::Pairwise && !prolog_before_search(&format.sink)
            }),
        queries_before: written.queries,
        subject_cache: written.subjects,
        x_pairwise_pool: x_pool.and_then(crate::report::x_pairwise_par::pool),
        x_tab_pool: x_pool.and_then(crate::algorithm::blastn::x_tab_par::pool),
        out: args.out.as_deref(),
        expect: crate::algorithm::blastn::coordination::determine_evalue(args),
        dust: args.dust.params().map(|(level, window, linker)| {
            (level as i32, window as u32 as i32, linker as u32 as i32)
        }),
        eff_db_length: if args.db_length() > 0 {
            args.db_length()
        } else {
            subjects.shown.total_letters as i64
        },
        archive: args.archive.as_ref(),
        batch_spans: results.batches,
        both_strands: args.strand == "both",
        minus_strand: QueryStrand::from_option(&args.strand) == QueryStrand::Minus,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-496
    // ```c++
    //         if (sid_in->IsLocal()) {
    //             string id_token;
    //             vector<string> title_tokens;
    //             title_tokens =
    //                 NStr::Split(CAlignFormatUtil::GetTitle(bh), " ", title_tokens);
    //             if(title_tokens.empty()){
    //                 id_token = NcbiEmptyString;
    //             } else {
    //                 id_token = title_tokens[0];
    //             }
    //
    //             if (id_token == NcbiEmptyString || parse_local) {
    //                 const CObject_id& obj_id = sid_in->GetLocal();
    //                 if (obj_id.IsStr())
    //                     id_token = obj_id.GetStr();
    // ```
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`).
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records.iter().map(record_title).collect();
    // `User specified sequence set (Input: PATH)` for a subject set, the `-db` value for a
    // database (`ShownSubjects`).
    let subject_title = subjects.shown.name.as_str();
    let hitlist_size = args.report_counts.hitlist_size;
    post_process_hits_and_write(
        hit_lists,
        hitlist_size,
        args.max_hsps_per_subject.unwrap_or(0),
        args.qcov_hsp_perc.unwrap_or(0.0),
        args.subject_besthit,
        results.lengths,
        outputs,
        args.verbose,
        &query_ids,
        output_formats,
        &query_titles,
        subject_title,
        None,
        &report_inputs,
        results.unsearched,
        results.warnings,
        sam_targets,
    )
}

/// The outputs that a query batch's search writes to: the outfmt 0 prolog of a file sink
/// before an error of the first batch (`write_pairwise_prologs`), and nothing otherwise.
/// `run_in_pool` gives them to the batches it searches; a batch searched ahead of the loop
/// (LOSAT_X_QBATCHPAR, never the first batch) has none.
pub(super) type BatchOutputs<'o, 'r> =
    Option<(&'o mut ReportOutputs<'r>, &'o [BlastnOutputFormat])>;

/// Why a search of the first query batch has its outputs (`BatchOutputs`).
pub(super) const FIRST_BATCH_OUTPUTS: &str =
    "run_in_pool searches the first query batch itself, with the outputs";

/// `outputs` again, for one call.
pub(super) fn reborrow_outputs<'o, 'r>(
    outputs: &'o mut BatchOutputs<'_, 'r>,
) -> BatchOutputs<'o, 'r> {
    outputs
        .as_mut()
        .map(|(outputs, formats)| (&mut **outputs, &**formats))
}
