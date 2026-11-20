//! The setup of a BLASTP call before its search (`run_call`, from its requested formats to the
//! shape of its scan): the reports that the call's formats ask for (`requested_reports`), the
//! queries' IDs, titles and lengths (`query_meta`), the query frames (`query_frames_of`), the
//! encoded subjects and their titles (`call_subjects`), the scoring values and X-drops
//! (`call_scoring`), the query batches, the lookup table and the query contexts (`call_lookup`),
//! the database's extent (`db_extent`), the search spaces (`call_search_spaces`), the cutoffs
//! (`prelim_cutoffs`), the statistics of each subject (`subject_parameters`), the hit lists
//! (`prelim_sinks`) and the scan's sizes (`scan_shape`). Each function is the statements of
//! `run_call` at its place, in their order; `run_call` calls them once each, in that order, with
//! its clocks and refusals between them, and its `?` passes their errors on.
//!
//! NCBI sets a search up in one function, the query data first, then the score block and the
//! lookup table:
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:132-136
//! ```c++
//! BlastSetupPreliminarySearchEx(CRef<IQueryFactory> qf,
//!                               CRef<CBlastOptions> options,
//!                               CConstRef<CPssmWithParameters> pssm,
//!                               BlastSeqSrc* seqsrc,
//!                               size_t num_threads)
//! ```
//!
//! This module is a child of `blast_engine` for the private items of that file it reads
//! (`BlastpReport`, `BlastpRanges`, `BlastpReports`, `PrelimOnly`, `BlastpRuntimeLookup`,
//! `BlastpContextLookupBounds`, the format kinds' functions, the query frames' and subjects'
//! encoders, the lookup builders, `x_prelim_block_plan`, `BlastpPrelimCutoffs`,
//! `BlastpSubjectParameters`, the cutoffs', search spaces' and scan's functions) and its sibling
//! modules `culling`, `database`, `psi_round`, `pssm_hooks`, `rps_hooks`, `x_ih_pool`,
//! `x_lookup_ahead`, `x_soft_subject` and `x_subject_meta`.

use std::ops::Range;
use std::sync::Arc;

use anyhow::{bail, Context, Result};

use super::database::BlastpDatabase;
use super::x_lookup_ahead::{XAheadBlock, XAheadLookup};
use super::x_subject_meta::XSubjectMeta;
// `compute_blosum62_ideal_karlin_params` is `blast_engine`'s import (from `stats`), which only this
// module reads.
use super::{
    blastp_aa_lookup, blastp_aa_lookup_on_block_spare, blastp_compressed_lookup,
    blastp_concatenated_query_length, blastp_context_lookup_bounds, blastp_is_archive,
    blastp_is_xml, blastp_offset_array_size, blastp_prelim_cutoffs, blastp_prelim_cutoffs_of,
    blastp_prelim_evalue,
    blastp_preliminary_cutoff_subject_length, blastp_query_anchored_view,
    blastp_query_frames_marked, blastp_seqalign_format, blastp_x_dropoff_init_for_context,
    blastp_xml2_kind, build_blastp_concatenated_query, build_redone_match_heaps,
    compute_blosum62_ideal_karlin_params, compute_diag_table_shape, culling, database,
    prepare_blastp_subjects_preserving_order, psi_round, pssm_hooks, rps_hooks, split_option_space,
    x_blocks_contexts, x_ih_pool, x_lookup_ahead, x_prelim_block_plan, x_soft_subject,
    BlastpContextLookupBounds, BlastpPrelimCutoffs, BlastpRanges, BlastpReport, BlastpReports,
    BlastpRuntimeLookup, BlastpSubjectParameters, BlastpSumStats, PrelimOnly,
};
use crate::algorithm::blastp::args::{
    BlastpCompositionMode, BlastpLookupTableType, ResolvedBlastpArgs,
};
use crate::algorithm::blastp::encoding::EncodedProtein;
use crate::algorithm::blastp::hsp::{BlastCompoHeap, BlastpHitList};
use crate::algorithm::blastp::x_prelim_compact::{x_expanded_lists, XPrelimHsp};
use crate::algorithm::blastp::{tabular, x_query_blocks, x_seed_fused, x_ungapped_sent};
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::x_bpclut_on;
use crate::algorithm::tblastx::lookup::QueryContext;
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::ReportOutputs;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range::RangedBatch;
use crate::config::ScoringMatrix;
use crate::report::record::ReportNeeds;
use crate::report::OutputFormat;
use crate::stats::search_space::SearchSpace;
use crate::stats::setup::SearchSpaceOverflow;
use crate::stats::{
    lookup_protein_gumbel_params, lookup_protein_params, lookup_protein_params_ungapped,
    BlastGumbelBlk, KarlinParams,
};
use crate::utils::threading::SearchPool;

/// What `query_frames_of`, `call_subjects` and `call_lookup` read of `run_call`: each field is
/// `run_call`'s value of the same name, by reference where it is not `Copy` (`subject_records`
/// are those of a FASTA subject set, none for a database; `x_block_inputs` are a query block's
/// subjects, query frames and subject scan ranges, LOSAT_X_BPQBLOCK). `'a` is the lifetime of
/// what the query frames and the subjects borrow; `'s` that of `run_call`'s own values (`args`,
/// `prelim`), whose borrow ends with the last call.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:132-134
/// ```c++
///     CRef<SBlastSetupData> setup_data =
///         BlastSetupPreliminarySearchEx(query_factory, options, pssm, seqsrc,
///                                       num_threads);
/// ```
pub(super) struct CallSetup<'s, 'a> {
    pub(super) args: &'s ResolvedBlastpArgs,
    pub(super) query_records: &'a [FastaRecord],
    pub(super) subject_records: &'a [FastaRecord],
    pub(super) prelim: &'s Option<PrelimOnly<'a>>,
    pub(super) query_chunk: bool,
    #[allow(clippy::type_complexity)]
    pub(super) x_block_inputs: Option<(
        &'a [EncodedProtein],
        &'a [Vec<QueryFrame>],
        Option<&'a [Vec<(i32, i32)>]>,
    )>,
    pub(super) blastp_parallel_pool: &'s SearchPool<'a>,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) blastp_num_threads: usize,
    pub(super) ranges: &'s BlastpRanges<'a>,
    pub(super) reports: &'s BlastpReports<'a>,
}

/// The refusal of the Seq-align formats (outfmt 8, 9, 11, 12) of an ungapped search, whose
/// Dense-diag Seq-align LOSAT's BLASTP does not write (`requested_reports`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1440
/// (`BLASTUngappedHspListToSeqAlign`)
fn ungapped_format_refusal(number: i32) -> String {
    format!("output format {number} with -ungapped (an ungapped search) is not supported by LOSAT's BLASTP")
}

/// The reports that the call's output formats ask for, in the formats' order (`BlastpReport`:
/// each `-outfmt` parsed as the single CLI `-outfmt` is; the query-anchored outfmt 1 to 4 are
/// BLASTP's only), and what they need beyond outfmt 0, 6 and 7 (`ReportNeeds`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2795-2798
/// ```c++
/// CFormattingArgs::ParseFormattingString(const CArgs& args,
///                                        EOutputFormat& fmt_type,
///                                        string& custom_fmt_spec,
///                                        string& custom_delim) const
/// ```
pub(super) fn requested_reports(
    args: &ResolvedBlastpArgs,
    outputs: &ReportOutputs<'_>,
) -> Result<(Vec<BlastpReport>, ReportNeeds)> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:2800-2803
    // ```c
    // if (args[kArgOutputFormat]) {
    //     string fmt_choice =
    //         NStr::TruncateSpaces(args[kArgOutputFormat].AsString());
    // ```
    // Each requested output format is parsed exactly as the single CLI `-outfmt` was.
    let requested_formats = outputs
        .formats
        .iter()
        .map(|format| -> Result<BlastpReport> {
            // The tabular reports with a field list or a delimiter, and outfmt 10 and 20.
            if let Some(spec) = tabular::custom_tabular_spec(format.outfmt) {
                return Ok(BlastpReport::CustomTabular(spec));
            }
            if let Some(serial) = blastp_seqalign_format(format.outfmt) {
                return Ok(BlastpReport::Seqalign(serial));
            }
            if blastp_is_xml(format.outfmt) {
                return Ok(BlastpReport::Xml);
            }
            if blastp_is_archive(format.outfmt) {
                return Ok(BlastpReport::Archive);
            }
            if let Some(kind) = blastp_xml2_kind(format.outfmt) {
                return Ok(BlastpReport::Xml2(kind));
            }
            // The query-anchored outfmt 1 to 4: BLASTP's (`Program::writes`; psiblast and rpsblast
            // refuse them, line RG PS15).
            if let Some(view) = blastp_query_anchored_view(args.app, format.outfmt)? {
                return Ok(BlastpReport::Pairwise(view));
            }
            // Outfmt 0, and outfmt 6 and 7 of the `std` fields and a tab (`TabularSpec::of`).
            let (outfmt, custom_fields) =
                OutputFormat::parse(format.outfmt).map_err(anyhow::Error::msg)?;
            if outfmt == OutputFormat::Pairwise && custom_fields.is_some() {
                bail!("blastp outfmt 0 does not accept custom field lists");
            }
            Ok(match outfmt {
                OutputFormat::Pairwise => {
                    BlastpReport::Pairwise(crate::report::pairwise::PairwiseView::Pairwise)
                }
                _ => BlastpReport::Text(outfmt),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1440
    // (`BLASTUngappedHspListToSeqAlign`: the Seq-align of an ungapped search is a Dense-diag per
    // HSP, which outfmt 8, 9 and 12 write and outfmt 11 archives)
    // The Seq-align formats of an ungapped search are not written yet (the Dense-diag writer of
    // branch `feature/np-o4b-dendiag`): refused by name before the first byte, with BLASTN's
    // wording (`blastn::blast_engine::run::process_options`).
    if args.ungapped {
        for (format, report) in outputs.formats.iter().zip(&requested_formats) {
            if matches!(report, BlastpReport::Seqalign(_) | BlastpReport::Archive) {
                bail!(
                    "{}",
                    ungapped_format_refusal(
                        crate::blastinput::app::parse_formatting_string(format.outfmt)?.number
                    )
                );
            }
        }
    }
    // What the requested formats need beyond outfmt 0, 6 and 7, decided once for the search
    // (`ReportNeeds`): the engine's positives of the composition-adjusted matches
    // (`engine_positives`). The engine's counts of the HSPs that the pyrrolysine recount
    // changes are kept for every format (`EngineCounts`: the Seq-aligns', which outfmt 0 reads
    // back with `-sorthits` and `-sorthsps`).
    let needs = crate::report::record::ReportNeeds {
        // The archive's Seq-aligns are those of outfmt 8.
        seqalign: requested_formats
            .iter()
            .any(|report| matches!(report, BlastpReport::Seqalign(_) | BlastpReport::Archive)),
    };
    Ok((requested_formats, needs))
}

/// The IDs, titles and lengths of the call's queries by query index, as the reports write them
/// (`query_meta`).
///
/// No NCBI counterpart: NCBI's formatter reads each query's ID, title and length from the
/// query's Bioseq where it writes them; LOSAT makes them once for the call.
pub(super) struct BlastpQueryMeta {
    pub(super) query_ids: Vec<Arc<[u8]>>,
    pub(super) query_titles: Vec<Arc<[u8]>>,
    pub(super) query_lengths: Vec<usize>,
}

/// The IDs, titles and lengths of `query_records` (`BlastpQueryMeta`): the titles and lengths that
/// NCBI writes in each query's preamble, and the query IDs of the tabular formats.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1492-1495
/// ```c++
///     CBlastFormatUtil::AcknowledgeBlastQuery(*bioseq, kFormatLineLength,
///                                             m_Outfile, m_BelieveQuery,
///                                             m_IsHTML, kIsTabularOutput,
///                                             results.GetRID());
/// ```
pub(super) fn query_meta(query_records: &[FastaRecord]) -> BlastpQueryMeta {
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:474-504
    // (`s_ReplaceLocalId`)
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`); the reports write the titles' bytes.
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
    let query_lengths: Vec<usize> = query_records
        .iter()
        .map(|record| record.seq().len())
        .collect();
    BlastpQueryMeta {
        query_ids,
        query_titles,
        query_lengths,
    }
}

/// The call's query frames (encoded, with SEG where it applies), in `x_owned_frames` where the
/// call makes them: a query block's from the search of all queries (LOSAT_X_BPQBLOCK), those
/// built ahead of the block (LOSAT_X_BPLKPAHEAD, taken from `x_prebuilt`), psiblast's later
/// rounds' (`psi_round`), or the queries' own.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:149-152
/// ```c++
///     // 1. Initialize the query data (borrow it from the factory)
///     CRef<ILocalQueryData> query_data(qf->MakeLocalQueryData(&*options));
///     retval->m_InternalData->m_Queries = query_data->GetSequenceBlk();
///     retval->m_InternalData->m_QueryInfo = query_data->GetQueryInfo();
/// ```
///
/// NCBI masks the queries (SEG) later, when it makes the score block.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158
/// ```c++
///     Int2 status = BLAST_MainSetUp(opts_memento->m_ProgramType,
/// ```
pub(super) fn query_frames_of<'a>(
    setup: &CallSetup<'_, 'a>,
    x_prebuilt: &mut Option<XAheadLookup>,
    x_owned_frames: &'a mut Vec<Vec<QueryFrame>>,
) -> &'a [Vec<QueryFrame>] {
    let CallSetup {
        args,
        query_records,
        x_block_inputs,
        ..
    } = *setup;
    match (x_block_inputs, x_prebuilt.as_mut()) {
        // EXPERIMENT (LOSAT_X_BPQBLOCK): a block of the preliminary stage takes its frames from the
        // search of all queries (its records, encoded with the same parameters).
        (Some((_, frames, _)), _) => frames,
        (None, Some(prebuilt)) => {
            *x_owned_frames = prebuilt.take_query_frames();
            x_owned_frames
        }
        // psiblast's later rounds: the query as the rounds before left it
        // (`psi_round::RoundQuery`).
        (None, None) if psi_round::searched_other_letters(args) => {
            *x_owned_frames = psi_round::round_query_frames(args);
            x_owned_frames
        }
        // Diagnostic (LOSAT_TIMING): this call's encoding and SEG go to the common line
        // (utils/xtiming.rs); frames built ahead (LOSAT_X_BPLKPAHEAD) are not marked. No NCBI
        // counterpart.
        (None, None) => {
            *x_owned_frames = blastp_query_frames_marked(args, query_records, true);
            x_owned_frames
        }
    }
}

/// The call's encoded subjects by subject index, the block search's subject IDs and titles
/// (`x_meta`, LOSAT_X_BPQBLOCK) and the subjects' titles as the reports write them
/// (`call_subjects`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:337-338
/// ```c++
///     /// The source of subject sequence data
///     CRef<TBlastSeqSrc> m_SeqSrc;
/// ```
pub(super) struct BlastpCallSubjects<'a> {
    pub(super) subjects: &'a [EncodedProtein],
    pub(super) x_meta: Option<&'a XSubjectMeta>,
    pub(super) subject_titles: std::borrow::Cow<'a, [Arc<[u8]>]>,
}

/// The subjects of a FASTA subject set encoded once, in `x_owned_subjects` where the call encodes
/// them (none for a database, whose sequences are read in blocks; the block search's encoded
/// once for all blocks, LOSAT_X_BPQBLOCK), the block search's subject IDs and titles, and the
/// subjects' titles (none for a call of the preliminary stage only).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_query_factory.cpp:101-102
/// ```c++
///     SetupSubjects_OMF(*m_QuerySource, program, &m_SeqBlkVector, &m_MaxLength);
///     m_NumSeqs = static_cast<Uint4>(m_QuerySource->Size());
/// ```
pub(super) fn call_subjects<'a>(
    setup: &CallSetup<'_, 'a>,
    x_owned_subjects: &'a mut Vec<EncodedProtein>,
) -> Result<BlastpCallSubjects<'a>> {
    let CallSetup {
        args,
        subject_records,
        prelim,
        x_block_inputs,
        blastp_parallel_pool,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
        ..
    } = *setup;
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1426
    // ```c
    // /* iterate over all subject sequences */
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
    //        continue;
    //    }
    // ```
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the block search encodes the subjects once
    // (`x_query_blocks::XBlockSubjects`), and so does the search of all queries whose preliminary
    // stage runs in blocks (`PrelimOnly::Block`); the same values.
    let subjects: &[EncodedProtein] = if let Some(encoded) = x_block_inputs
        .map(|(subjects, _, _)| subjects)
        .or_else(|| x_query_blocks::x_block_subjects(blastp_parallel_pool, subject_records))
    {
        encoded
    } else {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        let prepared_subjects = prepare_blastp_subjects_preserving_order(
            subject_records,
            args.num_threads,
            blastp_num_threads,
            blastp_parallel_pool,
        )?;
        #[cfg(any(
            not(feature = "parallel"),
            all(target_arch = "wasm32", not(feature = "wasm-threads"))
        ))]
        let prepared_subjects =
            prepare_blastp_subjects_preserving_order(subject_records, args.num_threads)?;
        *x_owned_subjects = prepared_subjects
            .into_iter()
            .map(|prepared_subject| prepared_subject.encoded)
            .collect();
        x_owned_subjects
    };
    // The subjects' titles by subject index, as the reports write them. EXPERIMENT
    // (LOSAT_X_BPQBLOCK): the block search's titles and tabular IDs, made once for all blocks
    // (`blast_engine/x_subject_meta.rs`); the same values, borrowed.
    let x_meta = x_query_blocks::x_block_subject_meta(blastp_parallel_pool, subject_records);
    let subject_titles: std::borrow::Cow<'_, [Arc<[u8]>]> = match x_meta {
        // A call of the preliminary stage only (`PrelimOnly`) writes no report.
        _ if prelim.is_some() => std::borrow::Cow::Owned(Vec::new()),
        Some(meta) => std::borrow::Cow::Borrowed(meta.titles()),
        None => std::borrow::Cow::Owned(
            subject_records
                .iter()
                .map(|record| Arc::from(record.title()))
                .collect(),
        ),
    };
    Ok(BlastpCallSubjects {
        subjects,
        x_meta,
        subject_titles,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:766-769
// ```c
// kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1833-1839
// ```c
// kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
// hsp->evalue = BLAST_KarlinStoE_simple(hsp->score, kbp[hsp->context],
//                                       hsp_list->hspcnt ? searchsp : 1);
// ```
/// The one owner of the Karlin block of the search's statistics for a context (decision D3 of
/// task O4b): the call's shared gapped block (`shared`: the matrix's gap costs row) in a gapped
/// search, and in an ungapped one the context's own `kbp_std` (`QueryContext::karlin_params`, with
/// the matrix's ungapped alpha and beta, `call_lookup`). The effective lengths, the cutoffs and
/// the XML statistics of the call read it here; the preliminary E-values, the bit scores and the
/// footer of an ungapped search (`prelim_kernel::ungapped_subject`, `traceback::UngappedFinish`,
/// the report's `ungapped_karlin`) read `ctx.karlin_params` themselves, which is the same block.
///
/// The selection is made once for a call by the caller of a per-subject or per-context loop:
/// `block_of` takes it as a constant, so the loops of a gapped search read `shared` as they did
/// before `-ungapped` (`subject_parameters`, `blastp_prelim_cutoffs`).
pub(super) fn context_block<'b>(
    args: &ResolvedBlastpArgs,
    shared: &'b KarlinParams,
    ctx: &'b QueryContext,
) -> &'b KarlinParams {
    if args.ungapped {
        block_of::<true>(shared, ctx)
    } else {
        block_of::<false>(shared, ctx)
    }
}

/// `context_block` for a search whose mode is known at compile time (`UNGAPPED` is
/// `args.ungapped`): the block that a loop of one mode reads, without a test of the mode.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1833-1839 (quoted at
/// `context_block`)
#[inline(always)]
pub(super) fn block_of<'b, const UNGAPPED: bool>(
    shared: &'b KarlinParams,
    ctx: &'b QueryContext,
) -> &'b KarlinParams {
    if UNGAPPED {
        &ctx.karlin_params
    } else {
        shared
    }
}

/// The call's ideal ungapped Lambda (the redo's), its matrix's gapped Karlin-Altschul parameters
/// and its gapped X-drops in raw score, preliminary and final (`call_scoring`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:324-325
/// ```c++
///     /// BLAST score block structure
///     CRef<TBlastScoreBlk> m_ScoreBlk;
/// ```
pub(super) struct BlastpScoring {
    pub(super) ideal_ungapped_lambda: f64,
    pub(super) gapped_params: KarlinParams,
    pub(super) x_drop_gapped: i32,
    pub(super) x_drop_gapped_final: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:517-527
// ```c
//        if (scoring_options->gapped_calculation) {
//        ...
//        } else {
//           ASSERT(sbp->kbp_gap == NULL);
//           /* for ungapped cases we do not have gbp filled */
//           if (sbp->gbp) {
//               sfree(sbp->gbp);
//               sbp->gbp=NULL;
// ```
/// The matrix's gapped Karlin block for the gap costs of `args`, or none for an ungapped search:
/// NCBI makes no `kbp_gap` for it and does not read its gap costs, which therefore need not be in
/// the matrix's table (`lookup_protein_params` stops for costs outside it). The owner for the
/// search (`call_scoring`) and for the reports of a run without a query batch
/// (`write_reports_without_queries`).
pub(super) fn gapped_block(args: &ResolvedBlastpArgs) -> Option<KarlinParams> {
    (!args.ungapped).then(|| lookup_protein_params(&args.scoring))
}

/// The call's `BlastpScoring` of `args`, with the refusal of an X-drop that NCBI's gapped kernels
/// cannot size before the search.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:422-424
/// ```c
/// Int2 BlastExtensionParametersNew(EBlastProgramType program_number,
///         const BlastExtensionOptions* options, BlastScoreBlk* sbp,
///         BlastQueryInfo* query_info, BlastExtensionParameters* *parameters)
/// ```
pub(super) fn call_scoring(args: &ResolvedBlastpArgs) -> Result<BlastpScoring> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2833-2852
    // ```c
    // Blast_ScoreBlkKbpIdealCalc(BlastScoreBlk* sbp)
    // {
    //    Blast_ResFreqStdComp(sbp, stdrfp);
    //    BlastScoreFreqCalc(sbp, sfp, stdrfp, stdrfp);
    //    Blast_KarlinBlkUngappedCalc(sbp->kbp_ideal, sfp);
    // }
    // ```
    // BLOSUM62's block is computed once per process; another matrix's is the one its lookup
    // table's query preparation computes (`ideal_karlin_params_of`). The redo reads its λ only.
    let ideal_ungapped_lambda = if args.scoring.matrix == ScoringMatrix::Blosum62 {
        compute_blosum62_ideal_karlin_params()
            .map_err(anyhow::Error::msg)
            .context("failed to calculate NCBI blastp kbp_ideal parameters")?
            .lambda
    } else {
        crate::algorithm::tblastx::lookup::ideal_karlin_params_of(args.scoring.matrix).lambda
    };
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_hits.c:1833-1839
    // ```c
    // kbp = (gapped_calculation ? sbp->kbp_gap : sbp->kbp);
    // hsp->evalue = BLAST_KarlinStoE_simple(hsp->score, kbp[hsp->context],
    //                                       hsp_list->hspcnt ? searchsp : 1);
    // ```
    // An ungapped search has no gapped Karlin block (`-gapopen` and `-gapextend` are ignored,
    // blast_setup.c:517-527), and so no gapped X-drops: its statistics are those of each
    // context's own block (`context_block`).
    let Some(gapped_params) = gapped_block(args) else {
        return Ok(BlastpScoring {
            ideal_ungapped_lambda,
            gapped_params: KarlinParams::default(),
            x_drop_gapped: 0,
            x_drop_gapped_final: 0,
        });
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:457-463
    // ```c
    // params->gap_x_dropoff = (Int4)
    //     (options->gap_x_dropoff*NCBIMATH_LN2 / min_lambda);
    // params->gap_x_dropoff_final = (Int4)
    //     MAX(options->gap_x_dropoff_final*NCBIMATH_LN2 / min_lambda,
    //         params->gap_x_dropoff);
    // ```
    // The X-drops in bits are `-xdrop_gap` and `-xdrop_gap_final` or BLASTP's defaults
    // (`ResolvedBlastpArgs`), converted as the protein programs convert them (`core::x_drop`); the
    // contexts of a batch share the matrix's gapped Lambda, the smallest one. An X-drop that NCBI's
    // gapped kernels cannot size is refused before the search.
    let x_drop_bits = (args.x_drop_gapped_bits, args.x_drop_final_bits);
    let (x_drop_gapped, x_drop_gapped_final) =
        crate::core::x_drop::gapped_x_dropoffs(x_drop_bits.0, x_drop_bits.1, gapped_params.lambda);
    // The final X-drop reaches the unscaled `ALIGN_EX` only on the ordinary traceback (mode 0,
    // after IDENTITY's reset too; the same choice as `ordinary_traceback` in `run_call`). A PSSM
    // search converts its X-drops with the PSSM's Lambda and checks them there
    // (`pssm_hooks::pssm_prelim`).
    if !args.ungapped && args.psi_search.is_none() {
        crate::core::x_drop::check_gapped_x_dropoffs(
            args.app.name(),
            x_drop_bits,
            (x_drop_gapped, x_drop_gapped_final),
            (args.scoring.gap_open, args.scoring.gap_extend),
            args.comp_based_stats.mode == BlastpCompositionMode::NoCompositionBasedStats
                || args.rps_search.is_some(),
        )?;
    }
    Ok(BlastpScoring {
        ideal_ungapped_lambda,
        gapped_params,
        x_drop_gapped,
        x_drop_gapped_final,
    })
}

/// The call's query batches and its plan of query blocks (LOSAT_X_BPQBLOCK), its lookup table and
/// query contexts, the contexts' length bounds, the concatenated query, the context of each query
/// offset, which queries are valid and which batches are not searched (`call_lookup`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/api/setup_factory.hpp:327-329
/// ```c++
///     /// Lookup table, usually only needed in the preliminary stage of the
///     /// search, but for PHI-BLAST it's also needed in the traceback stage.
///     CRef<TLookupTableWrap> m_LookupTable;
/// ```
pub(super) struct BlastpCallLookup {
    pub(super) query_batches: Vec<RangedBatch>,
    pub(super) x_blocks: Option<Vec<Range<usize>>>,
    pub(super) lookup: BlastpRuntimeLookup,
    pub(super) contexts: Vec<QueryContext>,
    pub(super) context_lookup_bounds: BlastpContextLookupBounds,
    pub(super) concatenated_query: Vec<u8>,
    pub(super) x_context_table: x_seed_fused::XContextTable,
    pub(super) query_valid: Vec<bool>,
    pub(super) query_batch_skipped: Vec<bool>,
}

/// The call's `BlastpCallLookup`: its query batches and block plan, then its lookup table and
/// contexts (rpsblast's contexts for the database's table, psiblast's PSSM table, a block plan's
/// contexts, the table of `-lookup_table_type`, or one built ahead of the block,
/// LOSAT_X_BPLKPAHEAD, taken from `x_prebuilt`), after which the build of the next block's table
/// starts (`x_ahead`), then the bounds, the concatenated query and the context table, and the
/// queries' validity.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:47-53
/// ```c
/// Int2 LookupTableWrapInit(BLAST_SequenceBlk* query,
///         const LookupTableOptions* lookup_options,
///         const QuerySetUpOptions* query_options,
///         BlastSeqLoc* lookup_segments, BlastScoreBlk* sbp,
///         LookupTableWrap** lookup_wrap_ptr, const BlastRPSInfo *rps_info,
///         Blast_Message* *error_msg,
///         BlastSeqSrc* seqsrc)
/// ```
pub(super) fn call_lookup(
    setup: &CallSetup<'_, '_>,
    query_frames: &[Vec<QueryFrame>],
    query_lengths: &[usize],
    x_prebuilt: &mut Option<XAheadLookup>,
    x_ahead: &XAheadBlock<'_>,
) -> Result<BlastpCallLookup> {
    let CallSetup {
        args,
        query_records,
        prelim,
        query_chunk,
        blastp_parallel_pool,
        ranges,
        reports,
        ..
    } = *setup;
    let batch_size = args.query_batch_size()?;
    // The batches count whole records (`seq_range::QueryInput::batches`); the batches
    // here are those that NCBI searches (`read_query_batches`).
    let query_batches = ranges.query_input.batches(batch_size);
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage of the search of all queries in blocks of
    // whole batches, each block a call of its own with its own lookup table (`PrelimOnly::Block`,
    // `x_prelim_block_plan`); `None`: one table for all queries. Planned only by the search of all
    // queries (`BlastpReports::x_plan_blocks`).
    let x_blocks = if prelim.is_none() && reports.x_plan_blocks {
        x_prelim_block_plan(
            args,
            query_records,
            ranges.query_input,
            batch_size,
            blastp_parallel_pool,
        )
    } else {
        None
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2746-2748
    // ```c
    // Boolean check_ideal =
    //    (program == eBlastTypeBlastx || program == eBlastTypeTblastx ||
    //     program == eBlastTypeRpsTblastn);
    // ```
    // Exact blastp (`AaLookupTable`) builds no compressed table: LOSAT_X_BPCLUT does not act.
    if matches!(args.lookup_table_type, BlastpLookupTableType::AaLookupTable) && x_bpclut_on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", false);
    }
    let (lookup, mut contexts) = match args.lookup_table_type {
        // rpsblast: the queries' contexts; the search scans with the database's prebuilt table
        // (`rps_hooks::rps_contexts`).
        _ if args.rps_search.is_some() => rps_hooks::rps_contexts(args, query_frames),
        // psiblast `-in_pssm`: the lookup table of the PSSM (`pssm_hooks::pssm_lookup`).
        _ if args.psi_search.is_some() => pssm_hooks::pssm_lookup(args, query_frames),
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the preliminary stage in blocks (`x_blocks`): each block
        // builds its own table; this call has the contexts of all queries (`x_blocks_contexts`).
        _ if x_blocks.is_some() => x_blocks_contexts(args, query_frames, blastp_parallel_pool),
        // `blastp_aa_lookup`. EXPERIMENT (LOSAT_X_BPLKPAHEAD): or the table and contexts built
        // ahead of this block by the same call on the block search's spare vectors
        // (`blast_engine/x_lookup_ahead.rs`). EXPERIMENT (LOSAT_X_BPLKREUSE): or, where nothing is
        // built ahead, built here on the vectors of the block before
        // (`x_query_blocks::x_bp_aa_block_spare`).
        BlastpLookupTableType::AaLookupTable => match x_prebuilt.take() {
            Some(prebuilt) => x_lookup_ahead::x_use_prebuilt(prebuilt),
            None => match x_query_blocks::x_bp_aa_block_spare(blastp_parallel_pool) {
                Some(spare) => blastp_aa_lookup_on_block_spare(
                    args,
                    query_frames,
                    spare,
                    query_chunk,
                    blastp_parallel_pool,
                ),
                None => blastp_aa_lookup(args, query_frames, query_chunk, blastp_parallel_pool),
            },
        },
        // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_aalookup.c:1349-1352
        // ```c
        // s_CompressedAddNeighboringWords(lookup, new_alphabet->matrix->data,
        //                                 query, locations);
        // s_CompressedLookupFinalize(lookup);
        // ```
        // `blastp_compressed_lookup`. EXPERIMENT (LOSAT_X_BPLKPAHEAD): or the table and contexts built
        // ahead of this block by the same calls (`blast_engine/x_lookup_ahead.rs`).
        BlastpLookupTableType::CompressedAaLookupTable => match x_prebuilt.take() {
            Some(prebuilt) => x_lookup_ahead::x_use_prebuilt(prebuilt),
            None => blastp_compressed_lookup(
                args,
                query_frames,
                &query_batches,
                query_lengths,
                blastp_parallel_pool,
                if query_chunk {
                    x_query_blocks::XLookupFor::QueryChunk
                } else {
                    x_query_blocks::XLookupFor::Search
                },
                if query_chunk {
                    0
                } else {
                    x_query_blocks::x_block_appends(blastp_parallel_pool, query_records)
                },
            )?,
        },
    };
    // EXPERIMENT (LOSAT_X_BPLKPAHEAD): this block's table is ready; the build of the next block's
    // starts now, on the pool, beside this block's search (`blast_engine/x_lookup_ahead.rs`).
    x_ahead.lookup_ready(&lookup);
    let context_lookup_bounds = blastp_context_lookup_bounds(&contexts);
    // EXPERIMENT (LOSAT_X_BPQBLOCK): with the preliminary stage in blocks (`x_blocks`) this call
    // scans no subject (each block's call makes its own query and table), so it makes neither
    // the concatenated query nor the context table (both empty).
    let x_scans = x_blocks.is_none();
    let concatenated_query = if x_scans {
        build_blastp_concatenated_query(&contexts)
    } else {
        Vec::new()
    };
    // The context of every query offset, for the candidates of the fused seed loop (the binary
    // search gives the same index).
    // No NCBI counterpart: a table of BSearchContextInfo's results; it does not change any value
    // NCBI computes.
    let x_context_table = {
        let frame_bases: Vec<i32> = if x_scans {
            contexts.iter().map(|ctx| ctx.frame_base).collect()
        } else {
            Vec::new()
        };
        x_seed_fused::XContextTable::new(&frame_bases, concatenated_query.len())
    };
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/local_blast.cpp:177-180
    // ```c
    //     int status = m_PrelimSearch->CheckInternalData();
    //     if (status != 0)
    //     {
    //          // Search was not run, but we send back an empty CSearchResultSet.
    // ```
    // A query whose Karlin-Altschul parameters cannot be computed is invalid
    // (`is_valid`); a query batch (`GetQueryBatchSize`) of invalid queries only is not
    // searched, and its reports show -1 parameters.
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_setup_cxx.cpp:632-640
    // ```c
    //         } catch (const CException& e) {
    //             ...
    //             CRef<CSearchMessage> m
    //                 (new CSearchMessage(eBlastSevWarning, index, e.GetMsg()));
    //             messages[index].push_back(m);
    //             s_InvalidateQueryContexts(qinfo, index);
    //         }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_stat.c:2771-2772
    // ```c
    //       if ( !contexts[context].is_valid )
    //           continue;
    // ```
    // A query without letters (`Sequence contains no data`, blast_setup.hpp:190-198) is
    // invalidated when the queries are set up, before its Karlin-Altschul parameters would
    // be computed: it is never searched and has no message of its own about them.
    for context in &mut contexts {
        if context.aa_len == 0 {
            context.is_valid = false;
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:3128-3130
    // ```c
    //    else if (num_values > 0) {
    //      (*alpha) = alpha_arr[0];
    //      (*beta) = beta_arr[0];
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:814-819
    // ```c
    //          BLAST_GetAlphaBeta(sbp->name, &alpha, &beta,
    //                             scoring_options->gapped_calculation,
    //                             scoring_options->gap_open,
    //                             scoring_options->gap_extend,
    //                             sbp->kbp_std[index]);
    // ```
    // An ungapped search's effective lengths use each context's own block (`context_block`) with
    // the matrix's ungapped alpha and beta, the first row of its table (the contexts' blocks hold
    // their lambda, K and H only).
    if args.ungapped {
        let ungapped = lookup_protein_params_ungapped(args.scoring.matrix);
        for context in &mut contexts {
            context.karlin_params.alpha = ungapped.alpha;
            context.karlin_params.beta = ungapped.beta;
        }
    }
    let mut query_valid = vec![true; query_records.len()];
    for context in &contexts {
        if !context.is_valid {
            query_valid[context.q_idx as usize] = false;
        }
    }
    let mut query_batch_skipped = vec![false; query_records.len()];
    for batch in &query_batches {
        let range = batch.searched.clone();
        let skipped = query_valid[range.clone()].iter().all(|valid| !valid);
        for flag in &mut query_batch_skipped[range] {
            *flag = skipped;
        }
    }
    Ok(BlastpCallLookup {
        query_batches,
        x_blocks,
        lookup,
        contexts,
        context_lookup_bounds,
        concatenated_query,
        x_context_table,
        query_valid,
        query_batch_skipped,
    })
}

/// The database's extent (`db_extent`): its total length, whether the traceback sorts the
/// results, the gapped Gumbel block (none with `-old_fsc`), the number of subjects, the length of
/// the effective lengths and the shortest subject length of the cutoffs (`min_subject_length`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:904-906
/// ```c
///    Uint4 min_subject_length;
///    Int8 total_length = -1;
///    Int4 num_seqs = -1;
/// ```
pub(super) struct BlastpDbExtent {
    /// NCBI's `Int8` (`BlastpDatabase::total_letters`, or the sum of the subject set's records).
    pub(super) total_db_len: i64,
    pub(super) sorts_results: bool,
    pub(super) gapped_gumbel: Option<BlastGumbelBlk>,
    pub(super) num_subjects: i32,
    pub(super) eff_db_length: i64,
    pub(super) cutoff_subject_length: i32,
}

/// The call's `BlastpDbExtent` (`database` none: the FASTA subject set `subjects`), with the
/// refusal of statistics that count no sequence (`-old_fsc`); `run_call` checks the average
/// subject length after it (its error writes NCBI's prolog first).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:888-889
/// ```c
/// BLAST_GapAlignSetUp(EBlastProgramType program_number,
///     const BlastSeqSrc* seq_src,
/// ```
pub(super) fn db_extent(
    args: &ResolvedBlastpArgs,
    database: Option<&BlastpDatabase>,
    subjects: &[EncodedProtein],
) -> Result<BlastpDbExtent> {
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:766-769
    // ```c
    // kbp_ptr = (scoring_options->gapped_calculation ? sbp->kbp_gap_std : sbp->kbp);
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:182-187
    // ```c++
    // static Int8
    // s_SeqDbGetTotLen(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetTotalLength();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-181
    // ```c++
    //     if(dbscan_mode)
    //     {
    // 	ITERATE(vector<BLAST_SequenceBlk*>, iter, m_ivSeqBlkVec)
    // 	{
    //             m_iTotalLength += (Int8) (*iter)->length;
    // 	}
    //     }
    // ```
    // The total length of a database is the one of its index files, not a sum over
    // sequences read; NCBI's `Int8` on every target (an alias file's `LENGTH` is not bounded by
    // the letters LOSAT holds). A subject set's lengths are in memory (at most `isize::MAX`).
    let total_db_len = match database {
        None => subjects
            .iter()
            .map(|subject| subject.aa_len as i64)
            .sum::<i64>(),
        Some(database) => database.total_letters(),
    };
    // NCBI reference: /mnt/c/Users/genom/GitHub/ncbi-blast/c++/src/algo/blast/core/blast_setup.c:914-921
    // ```c
    // /* Set the database length for new FSC */
    // if (sbp->gbp) {
    //     sbp->gbp->db_length =
    //         (Blast_SubjectIsTranslated(program_number)) ? dbl/3 : dbl;
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:913-923
    // ```c
    //       /* Set the database length for new FSC */
    //       if (sbp->gbp) {
    //           Int8 dbl = total_length;
    //           /* Override the database length or set one (e.g., blast2seq). */
    //           if (eff_len_options->db_length) {
    //               dbl = eff_len_options->db_length;
    //           }
    //           sbp->gbp->db_length =
    //               (Blast_SubjectIsTranslated(program_number))?
    //               dbl/3 : dbl;
    //       }
    // ```
    // `-dbsize` (any value but 0) replaces the total length in the Gumbel block. A database's
    // length is the setup's total before its fallback to the length of OID 0
    // (`DbStatistics::gumbel_db_length`): 0 for a `LENGTH` of 0, negative for one of 2^63 or
    // more, where every Spouge e-value is negative (written 0; line RG D1a-3).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1768-1772 (quoted at
    // `DbStatistics::sorts_results`): whether the traceback sorts the results at its end. A
    // FASTA subject set that finds hits has letters.
    let sorts_results = database.is_none_or(|database| database.statistics().sorts_results());
    let gumbel_db_length = if args.dbsize != 0 {
        args.dbsize
    } else {
        database.map_or(total_db_len, |database| {
            database.statistics().gumbel_db_length
        })
    };
    let gapped_gumbel = if args.has_gumbel() {
        Some(
            lookup_protein_gumbel_params(&args.scoring, gumbel_db_length)
                .context("missing NCBI gumbel parameters for supported blastp scoring")?,
        )
    } else {
        None
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:163-168
    // ```c++
    // static Int4
    // s_SeqDbGetNumSeqs(void* seqdb_handle, void*)
    // {
    //     CSeqDB & seqdb = **(TSeqDBData *) seqdb_handle;
    //     return seqdb.GetNumSeqs();
    // }
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:925-938 (`BLAST_GapAlignSetUp`)
    // ```c
    //       if (total_length > 0) {
    //           num_seqs = BlastSeqSrcGetNumSeqsStats(seq_src);
    //           if (num_seqs <= 0)
    //               num_seqs = BlastSeqSrcGetNumSeqs(seq_src);
    //       } else {
    //           /* Not a database search; each subject sequence is considered
    //              individually */
    //           ...
    //           num_seqs = 1;
    //       }
    // ```
    // A subject set of no letters at all (every record empty) is not a database search to
    // NCBI: it counts one sequence, whatever the number of records (visible only with
    // `-dbsize`, the effective search space of the XML statistics and the outfmt 0 footer).
    // NCBI's `Int4`: negative for an alias file's negative `NSEQ` or one of 2^31 or more.
    let num_subjects: i32 = match database {
        None if total_db_len == 0 => 1,
        None => i32::try_from(subjects.len().max(1)).unwrap_or(i32::MAX),
        Some(database) => rps_hooks::dbsize_volume_num_seqs(args, database)
            .map_or_else(|| database.statistics_num_sequences(), |count| count as i32),
    };
    // Without the Gumbel block (`-searchsp`, `OLD_FSC`, an ungapped search: `has_gumbel`) NCBI's
    // shortest subject length is the statistics' average (blast_setup.c:975, quoted at the check
    // of the average in `run_call`), which divides by their number of sequences.
    if let Some(database) = database.filter(|_| !args.has_gumbel()) {
        database.check_statistics_count(args.app.name())?;
        // RPS-BLAST's first volume of none of the profiles (`rps_hooks::dbsize_volume_num_seqs`)
        // gives the same division by 0: the same refusal (DECISIONS 14 of line RG), not a panic.
        if num_subjects == 0 {
            return Err(crate::algorithm::common::database::refusal(
                crate::algorithm::common::database::STATISTICS_COUNT_ZERO,
                args.app.name(),
            ));
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:716-721
    // ```c
    //    /* use overriding value from effective lengths options or the real value
    //       from effective lengths parameters. */
    //    if (eff_len_options->db_length > 0)
    //       db_length = eff_len_options->db_length;
    //    else
    //       db_length = eff_len_params->real_db_length;
    // ```
    // `-dbsize` above 0 replaces the total length behind the effective lengths.
    let eff_db_length: i64 = if args.dbsize > 0 {
        args.dbsize
    } else {
        total_db_len
    };
    // NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:128-136
    // ```c++
    //             m_NumSeqs     = x_GetNumSeqs();
    //             m_TotalLength = x_GetTotalLength();
    //             m_MaxLength   = x_GetMaxLength();
    //             m_MinLength   = x_GetMinLength();
    //
    //             // Do not bother scanning the db... it would be slow
    //             // FIXME: future implementation should probably have the
    //             // shortest length encoded in index file...
    //             if ( m_MinLength <= 0)  m_MinLength = BLAST_SEQSRC_MINLENGTH;
    // ```
    // A database answers `BLAST_SEQSRC_MINLENGTH` (10) for its shortest sequence (an
    // index file stores none, and the alias key `MIN_SEQ_LENGTH` is rejected by the
    // reader), whatever its sequences are, unless a list made NCBI scan its totals (then
    // the scanned shortest length, line DI-8); a FASTA subject set answers its real minimum,
    // raised to 10.
    // NCBI's `min_subject_length` (blast_setup.c:969-976, quoted at its check in `run_call`),
    // which the hit-saving parameters read: the shortest length with the Gumbel block, else
    // (`-searchsp`, `OLD_FSC`) the statistics' average, which only psiblast's `-sum_stats true`
    // reads (its cutoff, `sum_stats_hit_cutoff`; the cutoffs without the block read the search
    // space).
    let cutoff_subject_length = match (database, gapped_gumbel.is_some()) {
        // `setup_subject_length` without the block: the average as NCBI's `Int4` (0 where it is
        // 0 or the count is 0, which `run_call` refuses for a search after this function).
        (_, false) => {
            crate::stats::setup::setup_subject_length(total_db_len, num_subjects, None).unwrap_or(0)
        }
        (None, true) => blastp_preliminary_cutoff_subject_length(subjects),
        (Some(database), true) => database.cutoff_subject_length(),
    };
    Ok(BlastpDbExtent {
        total_db_len,
        sorts_results,
        gapped_gumbel,
        num_subjects,
        eff_db_length,
        cutoff_subject_length,
    })
}

/// What `call_search_spaces`, `prelim_cutoffs` and `subject_parameters` read of `run_call`: each
/// field is `run_call`'s value of the same name; the arguments, the contexts and the two Karlin
/// blocks by reference, as `run_call`'s closures read them.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:894-895
/// ```c
///     BlastQueryInfo* query_info,
///     BlastScoreBlk* sbp,
/// ```
#[derive(Clone, Copy)]
pub(super) struct CallStatistics<'s> {
    pub(super) args: &'s ResolvedBlastpArgs,
    pub(super) database: Option<&'s BlastpDatabase>,
    pub(super) ranges: &'s BlastpRanges<'s>,
    pub(super) contexts: &'s [QueryContext],
    pub(super) gapped_params: &'s KarlinParams,
    pub(super) gapped_gumbel: &'s Option<BlastGumbelBlk>,
    pub(super) total_db_len: i64,
    pub(super) num_subjects: i32,
    pub(super) eff_db_length: i64,
    pub(super) cutoff_subject_length: i32,
}

/// The call's search spaces (`call_search_spaces`): the setup's of each context and their exact
/// values, whether the preliminary stage computes the statistics of each subject, the option
/// space of a context (which a query chunk's search reads again), the preliminary stage's option
/// spaces where it reads them, and its own spaces where they differ from the setup's.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:699-700
/// ```c
/// Int2 BLAST_CalcEffLengths (EBlastProgramType program_number,
///    const BlastScoringOptions* scoring_options,
/// ```
pub(super) struct BlastpSearchSpaces<F> {
    pub(super) search_spaces: Vec<SearchSpace>,
    pub(super) exact_search_spaces: Vec<Option<i64>>,
    pub(super) per_subject: bool,
    pub(super) option_space_of: F,
    pub(super) preliminary_option_spaces: Option<Vec<i64>>,
    pub(super) preliminary_spaces: Option<(Vec<SearchSpace>, Vec<Option<i64>>)>,
}

/// The call's `BlastpSearchSpaces`: the setup's spaces, then the database's statistics, the
/// options' spaces and the preliminary stage's own spaces.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:946-947
/// ```c
///    if ((status = BLAST_CalcEffLengths(program_number, scoring_options,
///                      *eff_len_params, sbp, query_info, NULL)) != 0)
/// ```
pub(super) fn call_search_spaces<'s>(
    statistics: &CallStatistics<'s>,
) -> Result<BlastpSearchSpaces<impl Fn(&QueryContext) -> Result<i64, SearchSpaceOverflow> + 's>> {
    let CallStatistics {
        args,
        database,
        ranges,
        contexts,
        gapped_params,
        total_db_len,
        num_subjects,
        eff_db_length,
        ..
    } = *statistics;
    // The effective search space of a database or a FASTA subject set is computed in NCBI's
    // `Int8` (`blastp_database_search_space`, the one owner `stats::setup::calc_eff_lengths`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:942-944
    // ```c
    //          } else {
    //              BLAST_Cutoffs(&new_cutoff, &evalue, kbp, searchsp, FALSE, 0);
    //          }
    // ```
    // Without the Gumbel block (`!has_gumbel`: `-searchsp`, `OLD_FSC` or an ungapped search) the hit cutoffs read
    // each context's search space, which the exact computation hands them
    // (`exact_search_spaces`; `blastp_cutoff_score_max`).
    let options_searchsp = args.searchsp.unwrap_or(0);
    // NCBI reference: c++/src/algo/blast/core/blast_setup.c:729-732
    // ```c
    //    if (db_length == 0 &&
    //        !BlastEffectiveLengthsOptions_IsSearchSpaceSet(eff_len_options)) {
    //       return 0;
    //    }
    // ```
    // When every subject is empty, NCBI leaves each query's effective search space and
    // length adjustment at 0 (the outfmt 0 footers print 0).
    let mut exact_search_spaces: Vec<Option<i64>> = Vec::new();
    let search_spaces: Vec<SearchSpace> = contexts
        .iter()
        .enumerate()
        .map(|(context_index, ctx)| -> Result<SearchSpace> {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:676-697
            // ```c
            //     if (eff_len_options->num_searchspaces == 0) {
            //         retval = 0;
            //     } else if (eff_len_options->num_searchspaces == 1) {
            //         ...
            //         retval = eff_len_options->searchsp_eff[0];
            //     } else if (eff_len_options->num_searchspaces > 1) {
            //         ASSERT(context_index < eff_len_options->num_searchspaces);
            //         retval = eff_len_options->searchsp_eff[context_index];
            // ```
            // A query chunk's search takes its batch's search spaces as options
            // (`search_split_query_batches`), by the context's index in the chunk, unless the
            // user's `-searchsp` is set (`SplitQuery_SetEffectiveSearchSpace` leaves it).
            let given_searchsp = match ranges.split_search_spaces {
                Some(batch) => split_option_space(batch, context_index, options_searchsp),
                None => options_searchsp,
            };
            if given_searchsp != 0 {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:826-846
                // ```c
                //          if (effective_search_space == 0) {
                //          ...
                //       query_info->contexts[index].eff_searchsp = effective_search_space;
                // ```
                exact_search_spaces.push(Some(given_searchsp));
                return Ok(SearchSpace {
                    effective_space: given_searchsp as f64,
                    length_adjustment: if ctx.aa_len > 0 {
                        i64::from(
                            crate::stats::setup::calc_eff_lengths(
                                ctx.aa_len as i64,
                                eff_db_length,
                                num_subjects,
                                context_block(args, gapped_params, ctx),
                                options_searchsp,
                                args.app.name(),
                            )?
                            .length_adjustment,
                        )
                    } else {
                        0
                    },
                });
            }
            if eff_db_length == 0 {
                exact_search_spaces.push(None);
                return Ok(SearchSpace {
                    effective_space: 0.0,
                    length_adjustment: 0,
                });
            }
            let (space, exact) = database::blastp_database_search_space(
                ctx.aa_len,
                eff_db_length,
                num_subjects,
                context_block(args, gapped_params, ctx),
                args.app.name(),
            )?;
            exact_search_spaces.push(Some(exact));
            Ok(space)
        })
        .collect::<Result<Vec<_>>>()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-220 (quoted at
    // `stats::setup::preliminary_option_space`): the preliminary stage's setup reads the options'
    // search spaces of the statistics before their fallback to OID 0, which differ from the
    // setup's own only for a total length that is not positive, and which the statistics of each
    // subject read too (`DbStatistics::per_subject`, blast_engine.c:1434-1443).
    let statistics = database.map(|database| database.statistics());
    // RPS-BLAST's engine has no loop over subjects (blast_engine.c:1384-1393).
    let per_subject = args.rps_search.is_none() && statistics.is_some_and(|s| s.per_subject());
    let (option_total, option_count) = statistics
        .map_or((total_db_len, num_subjects), |statistics| {
            (statistics.gumbel_db_length, statistics.request_count)
        });
    // The option space of a context of this search's queries (`preliminary_option_space`): 0 for
    // a context that is not valid (its Karlin-Altschul block failed on the query, or it has no
    // letters) unless `-searchsp` gives it.
    let option_space_of = move |ctx: &QueryContext| {
        crate::stats::setup::preliminary_option_space(
            if ctx.is_valid { ctx.aa_len as i64 } else { 0 },
            option_total,
            option_count,
            args.dbsize,
            options_searchsp,
            false,
            context_block(args, gapped_params, ctx),
            args.app.name(),
        )
    };
    // A query chunk's search reads the whole batch's option spaces by its own context index
    // (`BlastpRanges::split_search_spaces`, made by `search_split_query_batches`), as its setup's
    // spaces above do.
    let preliminary_option_spaces: Option<Vec<i64>> =
        if option_total <= 0 && statistics.is_some() || per_subject {
            Some(match ranges.split_search_spaces {
                Some(batch) => (0..contexts.len())
                    .map(|ctx_idx| split_option_space(batch, ctx_idx, options_searchsp))
                    .collect(),
                None => contexts
                    .iter()
                    .map(option_space_of)
                    .collect::<Result<_, _>>()?,
            })
        } else {
            None
        };
    // The preliminary stage's own spaces where they differ from the setup's (a total length that
    // is not positive): the setup's computation with the option spaces.
    let preliminary_spaces: Option<(Vec<SearchSpace>, Vec<Option<i64>>)> =
        match &preliminary_option_spaces {
            // `BLAST_CalcEffLengths` returns before the contexts for a length of 0 without an
            // option space (blast_setup.c:729-732): the setup's spaces stay.
            Some(option_spaces)
                if option_total <= 0
                    && statistics.is_some()
                    && (eff_db_length != 0 || option_spaces.iter().any(|&space| space != 0)) =>
            {
                let mut spaces = search_spaces.clone();
                let mut exact = exact_search_spaces.clone();
                for (ctx_idx, ctx) in contexts.iter().enumerate() {
                    if !ctx.is_valid || ctx.aa_len == 0 {
                        continue;
                    }
                    let lengths = crate::stats::setup::calc_eff_lengths(
                        ctx.aa_len as i64,
                        eff_db_length,
                        num_subjects,
                        context_block(args, gapped_params, ctx),
                        option_spaces[ctx_idx],
                        args.app.name(),
                    )?;
                    spaces[ctx_idx] = SearchSpace {
                        effective_space: lengths.eff_searchsp as f64,
                        length_adjustment: i64::from(lengths.length_adjustment),
                    };
                    exact[ctx_idx] = Some(lengths.eff_searchsp);
                }
                Some((spaces, exact))
            }
            _ => None,
        };
    Ok(BlastpSearchSpaces {
        search_spaces,
        exact_search_spaces,
        per_subject,
        option_space_of,
        preliminary_option_spaces,
        preliminary_spaces,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:203-209
// ```c++
// s_SeqDbGetAvgLength(void* seqdb_handle, void* ignoreme)
// {
//    Int8 total_length = s_SeqDbGetTotLen(seqdb_handle, ignoreme);
//    Int4 num_seqs = MAX(1, s_SeqDbGetNumSeqs(seqdb_handle, ignoreme));
//
//    return (Int4) (total_length/num_seqs);
// }
// ```
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1372-1374
// ```c
//     BlastInitialWordParametersNew(program_number, word_options,
//       hit_params, lookup_wrap, sbp, query_info,
//       BlastSeqSrcGetAvgSeqLen(seq_src), &word_params);
// ```
/// `BlastSeqSrcGetAvgSeqLen` of the call's subjects, the subject length of the word cutoffs of an
/// ungapped search (`blastp_prelim_cutoffs`): a database's (`SeqDb::average_length`), a FASTA
/// subject set's total over its records (the multi-sequence source's average).
fn average_subject_length(statistics: &CallStatistics<'_>) -> i32 {
    match statistics.database {
        Some(database) => database.statistics().average,
        None => (statistics.total_db_len / i64::from(statistics.num_subjects.max(1))) as i32,
    }
}

/// The preliminary stage's cutoffs (`prelim_cutoffs`): `BlastpPrelimCutoffs`' fields, the smallest
/// hit-saving cutoff of the valid contexts, the size of the hit lists and the preliminary e-value.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:734-738
/// ```c
/// BlastHitSavingParametersNew(EBlastProgramType program_number,
/// ...
///    Int4 avg_subj_length,
/// ```
pub(super) struct BlastpCallCutoffs {
    pub(super) hit_cutoff_scores: Vec<i32>,
    pub(super) word_cutoff_scores: Vec<i32>,
    pub(super) sum_stats: Option<BlastpSumStats>,
    pub(super) hit_cutoff_score_min: i32,
    pub(super) prelim_hitlist_size: usize,
    pub(super) prelim_evalue: f64,
}

/// The call's `BlastpCallCutoffs` with the preliminary stage's exact search spaces.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:982-983
/// ```c
///    if ((status = BlastHitSavingParametersNew(program_number, hit_options, sbp,
///                                              query_info, min_subject_length,
/// ```
pub(super) fn prelim_cutoffs(
    statistics: &CallStatistics<'_>,
    preliminary_exact_spaces: &[Option<i64>],
) -> Result<BlastpCallCutoffs> {
    let CallStatistics {
        args,
        contexts,
        gapped_params,
        gapped_gumbel,
        cutoff_subject_length,
        ..
    } = *statistics;
    // The setup's cutoffs; psiblast's `-sum_stats true` lowers the cutoff of the HSPs a search
    // keeps (`cutoff_score`) to the sum statistics' cutoff where that is smaller, and the
    // preliminary stage links the HSPs of each subject (`BlastpSumStats`). Only an ungapped
    // search's word cutoffs read the sequence source's average subject length.
    let word_subject_length = if args.ungapped {
        average_subject_length(statistics)
    } else {
        cutoff_subject_length
    };
    let BlastpPrelimCutoffs {
        hit_cutoff_scores,
        word_cutoff_scores,
        sum_stats,
    } = blastp_prelim_cutoffs(
        args,
        contexts,
        cutoff_subject_length,
        word_subject_length,
        args.comp_based_stats.mode,
        preliminary_exact_spaces,
        gapped_params,
        gapped_gumbel.as_ref(),
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:979-990
    // ```c
    //       /* scale up the computed cutoffs, and find the smallest one */
    //       for (context = query_info->first_context;
    //                              context <= query_info->last_context; ++context) {
    //
    //          if (query_info->contexts[context].is_valid) {
    // ...
    //             cutoff_min = MIN(cutoff_min,
    //                              params->cutoffs[context].cutoff_score);
    //          }
    //       }
    //       params->cutoff_score_min = cutoff_min;
    // ```
    // `hit_params->cutoff_score_min`, which the redo reads with sum statistics.
    let hit_cutoff_score_min = contexts
        .iter()
        .zip(&hit_cutoff_scores)
        .filter(|(ctx, _)| ctx.is_valid)
        .map(|(_, &cutoff)| cutoff)
        .min()
        .unwrap_or(i32::MAX);

    // NCBI computes the size in its `int` (`blastn/hsp.rs` `get_prelim_hitlist_size`): a
    // -max_target_seqs from 2^31 - 24 wraps to 2..48; a size that is not positive (NCBI
    // crashes) was rejected with the options (`validate_requested_blastp_support`).
    let prelim_hitlist_size =
        usize::try_from(crate::algorithm::blastn::hsp::get_prelim_hitlist_size(
            args.max_target_seqs,
            args.comp_based_stats.is_enabled(),
            !args.ungapped,
        ))
        .ok()
        .filter(|&size| size > 0)
        .context("BLASTP preliminary hit list size is not positive")?;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:1008-1018
    // ```c
    // BLAST_OneSubjectUpdateParameters(..., subject_length, ...);
    // BlastHitSavingParametersUpdate(program_number, sbp, query_info,
    //                                subject_length, 0, hit_params);
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:926-940
    // ```c
    // params->prelim_evalue = evalue;  /* evalue and prelim_evalue same if no CBS. */
    // ...
    // int cbs_stretch = (compositionBasedStats > 1) ? 5 : 1;
    // params->prelim_evalue = cbs_stretch*evalue;
    // ```
    // `params->prelim_evalue = evalue` without a Gumbel block (blast_parameters.c:931-935).
    let prelim_evalue = if gapped_gumbel.is_some() {
        blastp_prelim_evalue(args.evalue, args.comp_based_stats.mode)
    } else {
        args.evalue
    };
    Ok(BlastpCallCutoffs {
        hit_cutoff_scores,
        word_cutoff_scores,
        sum_stats,
        hit_cutoff_score_min,
        prelim_hitlist_size,
        prelim_evalue,
    })
}

/// The statistics of one subject at its length (`BlastpSubjectParameters`), from the preliminary
/// stage's option spaces and spaces: `run_call`'s `process_own_subject` calls it once per subject
/// with letters (`per_subject`) and notes its error, an overflow of NCBI's `Int8` search space.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1001-1002
/// ```c
/// Int2 BLAST_OneSubjectUpdateParameters(EBlastProgramType program_number,
///                     Uint4 subject_length,
/// ```
pub(super) fn subject_parameters<'s>(
    statistics: &CallStatistics<'s>,
    preliminary_option_spaces: &'s Option<Vec<i64>>,
    preliminary_search_spaces: &'s [SearchSpace],
    preliminary_exact_spaces: &'s [Option<i64>],
) -> impl Fn(usize) -> Result<BlastpSubjectParameters, SearchSpaceOverflow> + 's {
    // `CallStatistics` holds references and numbers only: the closure keeps its own copy.
    let statistics_of_call = *statistics;
    let args = statistics_of_call.args;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1015-1025
    // ```c
    //    /* Update cutoff scores in hit saving parameters */
    //    BlastHitSavingParametersUpdate(program_number, sbp, query_info, subject_length, 0, /* FIXME */
    //                                   hit_params);
    //
    //    if (word_params) {
    //       /* Update cutoff scores in initial word parameters */
    //       BlastInitialWordParametersUpdate(program_number, hit_params, sbp, query_info,
    //            subject_length, word_params);
    //       /* Update the parameters for linking HSPs, if necessary. */
    //       BlastLinkHSPParametersUpdate(word_params, hit_params, scoring_options->gapped_calculation);
    //    }
    // ```
    // The statistics of one subject (`per_subject`): each context's effective lengths
    // (`stats::setup::subject_eff_lengths`) and, through the setup's own `blastp_prelim_cutoffs`,
    // the hit-saving cutoffs of `BlastHitSavingParametersUpdate` with the subject's length and a
    // composition-based statistics of 0 (no stretch of the
    // e-value: `params->prelim_evalue = evalue`, blast_parameters.c:868, 937-938), and the
    // word cutoffs of `BlastInitialWordParametersUpdate` (gapped: the gap trigger, at most the
    // hit-saving cutoff before the sum statistics lower it, `cutoff_score_max`; the X-drops keep
    // their initial values, `x_dropoff_init` is not 0). With psiblast's `-sum_stats true` the
    // subject's sum statistics (`BlastpSumStats` at the subject's length, blast_parameters.c:949-976)
    // lower its hit-saving cutoffs, and the preliminary linking reads the subject's small-gap
    // cutoff (`BlastLinkHSPParametersUpdate`, blast_parameters.c:625-646) and search spaces. The
    // traceback reads the setup's values (blast_traceback.c:1598-1650: it updates its own copy of
    // the query information, which the traceback of a list does not read, and its hit-saving
    // cutoffs, which only the complexity adjustment reads, and links with the traceback's
    // parameters, which no word parameters update).
    // The mode is decided here, once for a subject, outside the loops over the contexts: the
    // instance of a gapped search reads the call's gapped block for every context and the word
    // cutoff's gapped rule, as it did before `-ungapped` (`context_block`, `block_of`).
    move |subject_length: usize| -> Result<
        BlastpSubjectParameters,
        crate::stats::setup::SearchSpaceOverflow,
    > {
        let option_spaces = preliminary_option_spaces
            .as_deref()
            .expect("the statistics of each subject have their preliminary option spaces");
        if args.ungapped {
            one_subject::<true>(
                statistics_of_call,
                option_spaces,
                preliminary_search_spaces,
                preliminary_exact_spaces,
                subject_length,
            )
        } else {
            one_subject::<false>(
                statistics_of_call,
                option_spaces,
                preliminary_search_spaces,
                preliminary_exact_spaces,
                subject_length,
            )
        }
    }
}

/// The statistics of one subject (`subject_parameters`) of a search whose mode is `UNGAPPED`
/// (`args.ungapped`): each valid context's effective lengths with the context's block
/// (`block_of`), then the cutoffs (`blastp_prelim_cutoffs_of`).
#[inline(always)]
fn one_subject<const UNGAPPED: bool>(
    statistics: CallStatistics<'_>,
    option_spaces: &[i64],
    preliminary_search_spaces: &[SearchSpace],
    preliminary_exact_spaces: &[Option<i64>],
    subject_length: usize,
) -> Result<BlastpSubjectParameters, SearchSpaceOverflow> {
    let CallStatistics {
        args,
        contexts,
        gapped_params,
        gapped_gumbel,
        num_subjects,
        ..
    } = statistics;
    let mut search_spaces = preliminary_search_spaces.to_vec();
    let mut exact_spaces = preliminary_exact_spaces.to_vec();
    for (ctx_idx, ctx) in contexts.iter().enumerate() {
        if !ctx.is_valid || ctx.aa_len == 0 {
            continue;
        }
        let lengths = crate::stats::setup::subject_eff_lengths(
            ctx.aa_len as i64,
            subject_length as i64,
            num_subjects,
            args.dbsize,
            option_spaces[ctx_idx],
            false,
            block_of::<UNGAPPED>(gapped_params, ctx),
            args.app.name(),
        )?;
        search_spaces[ctx_idx] = SearchSpace {
            effective_space: lengths.eff_searchsp as f64,
            length_adjustment: i64::from(lengths.length_adjustment),
        };
        exact_spaces[ctx_idx] = Some(lengths.eff_searchsp);
    }
    Ok(BlastpSubjectParameters {
        search_spaces,
        cutoffs: blastp_prelim_cutoffs_of::<UNGAPPED>(
            args,
            contexts,
            subject_length as i32,
            subject_length as i32,
            BlastpCompositionMode::NoCompositionBasedStats,
            &exact_spaces,
            gapped_params,
            gapped_gumbel.as_ref(),
        ),
    })
}

/// What the call's preliminary stage and redo write (`prelim_sinks`): the culling writer, the first
/// overflow of a subject's search space, the count of lists freed by the redo, the redo's heaps,
/// the preliminary hit lists (full and compact) and the diagnostic counters at their start.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:286-287
/// ```c++
/// CSetupFactory::CreateHspStream(const CBlastOptionsMemento* opts_memento,
///                                size_t number_of_queries,
/// ```
pub(super) struct BlastpPrelimSinks {
    pub(super) culling: Option<culling::PrelimCulling>,
    pub(super) subject_space_overflow: crate::stats::setup::SubjectSpaceOverflow,
    pub(super) x_prelim_freed: usize,
    pub(super) redone_matches: Vec<BlastCompoHeap>,
    pub(super) preliminary_hit_lists: Vec<Option<BlastpHitList>>,
    pub(super) x_expanded_at_start: u64,
    pub(super) x_pooled_at_start: u64,
    pub(super) x_compact_hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>>,
}

/// The call's `BlastpPrelimSinks` for its `query_count` queries (the culling writer only where the
/// preliminary stage does not run in blocks, `x_blocks`, LOSAT_X_BPQBLOCK).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:291-293
/// ```c++
///     return BlastHSPStreamNew(opts_memento->m_ProgramType,
///                              opts_memento->m_ExtnOpts, TRUE,
///                              static_cast<Int4>(number_of_queries), writer);
/// ```
pub(super) fn prelim_sinks(
    args: &ResolvedBlastpArgs,
    x_blocks: Option<&[Range<usize>]>,
    prelim_hitlist_size: usize,
    query_count: usize,
) -> BlastpPrelimSinks {
    // `-culling_limit`: NCBI's culling writer takes the place of the collector (`CreateHspWriter`,
    // setup_factory.cpp:318-342). The choice is made once for this call: the whole search, a block
    // of query batches, or one chunk of a split query, each with its own writer.
    // EXPERIMENT (LOSAT_X_BPQBLOCK): with the preliminary stage in blocks (`x_blocks`) each block's
    // call has its own writer, as NCBI's batch.
    let culling = (args.culling_limit > 0 && x_blocks.is_none())
        .then(|| culling::PrelimCulling::new(prelim_hitlist_size, args.culling_limit));
    // The first overflow of a subject's search space (`subject_parameters`), noted by the kernel.
    let subject_space_overflow = crate::stats::setup::SubjectSpaceOverflow::default();
    // `x_prelim_free.rs`.
    let x_prelim_freed = 0usize;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3074-3077
    // ```c
    //     inclusion_ethresh = (psiOptions /* this can be NULL for CBl2Seq */
    //                          ? psiOptions->inclusion_ethresh
    //                          : PSI_INCLUSION_ETHRESH);
    //     ASSERT(inclusion_ethresh != 0.0);
    // ```
    // The heaps hold psiblast's larger hit lists (`traceback_hitlist_size`).
    let redone_matches = build_redone_match_heaps(
        query_count,
        args.traceback_hitlist_size(),
        args.inclusion_ethresh,
    );
    let preliminary_hit_lists: Vec<Option<BlastpHitList>> =
        std::iter::repeat_with(|| None).take(query_count).collect();
    // The preliminary hit lists with the compact element (`x_prelim_compact.rs`), the lists the
    // subject scan keeps until the redo.
    let x_expanded_at_start = x_expanded_lists();
    let x_pooled_at_start = x_ih_pool::x_pooled();
    let x_compact_hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>> =
        std::iter::repeat_with(|| None).take(query_count).collect();
    BlastpPrelimSinks {
        culling,
        subject_space_overflow,
        x_prelim_freed,
        redone_matches,
        preliminary_hit_lists,
        x_expanded_at_start,
        x_pooled_at_start,
        x_compact_hit_lists,
    }
}

/// The shape of the call's scan (`scan_shape`): the X-drop of each context's ungapped extensions,
/// the window, the word size, the query of the inlined two-hit extension (LOSAT_X_BPUNGAPSENT),
/// the offset array's and the diagonal table's sizes, the interval tree's query span, and the
/// subjects' scan ranges of `-lcase_masking`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:110-112
/// ```c
/// Int2 BlastExtendWordNew(Uint4 query_length,
///                         const BlastInitialWordParameters * word_params,
///                         Blast_ExtendWord ** ewp_ptr)
/// ```
pub(super) struct BlastpScanShape<'a> {
    pub(super) x_dropoff_per_context: Vec<i32>,
    pub(super) window: i32,
    pub(super) word_size: i32,
    pub(super) x_ungap_sent: Option<x_ungapped_sent::XSentQuery>,
    pub(super) offset_array_size: i32,
    pub(super) diag_array_size: i32,
    pub(super) diag_mask: i32,
    pub(super) total_query_span: i32,
    pub(super) x_soft_subjects: Option<&'a [Vec<(i32, i32)>]>,
}

/// The call's `BlastpScanShape`, with the refusal of a window that the diagonal table cannot hold;
/// the scan ranges that this call makes are written to `x_owned_soft_subjects` and borrowed.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:159-160
/// ```c
/// BlastInitialWordParametersNew(EBlastProgramType program_number,
///    const BlastInitialWordOptions* word_options,
/// ```
pub(super) fn scan_shape<'a>(
    setup: &CallSetup<'_, 'a>,
    lookup: &BlastpRuntimeLookup,
    contexts: &[QueryContext],
    concatenated_query: &[u8],
    x_blocks: Option<&[Range<usize>]>,
    x_owned_soft_subjects: &'a mut Option<Vec<Vec<(i32, i32)>>>,
) -> Result<BlastpScanShape<'a>> {
    let CallSetup {
        args,
        subject_records,
        x_block_inputs,
        ..
    } = *setup;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_parameters.c:195-221
    // ```c
    // for (context = query_info->first_context;
    //                  context <= query_info->last_context; ++context) {
    //    ...
    //    kbp= sbp->kbp[context];
    //    ...
    //    p->cutoffs[context].x_dropoff_init =
    //        (Int4)(sbp->scale_factor *
    //               ceil(word_options->x_dropoff * NCBIMATH_LN2 / kbp->Lambda));
    // }
    // ```
    //
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_setup.c:318-323,417-420
    // ```c
    // sbp->kbp_std[index] = (Blast_KarlinBlk*)
    //     BlastMemDup(sbp->kbp_gap_std[0], sizeof(Blast_KarlinBlk));
    // sbp->kbp = sbp->kbp_std;
    // ...
    // Blast_KarlinBlkCopy(sbp->kbp_std[context], sbp->kbp_ideal);
    // sbp->kbp = sbp->kbp_std;
    // ```
    //
    // Default blastp uses the score-block ungapped Karlin blocks (`sbp->kbp`
    // / `sbp->kbp_std`), not per-query lookup composition estimates.
    let x_dropoff_per_context: Vec<i32> = contexts
        .iter()
        .map(|ctx| blastp_x_dropoff_init_for_context(ctx, args.x_drop_ungapped_bits))
        .collect();
    let window = args.window_size as i32;
    let word_size = lookup.word_length();
    // EXPERIMENT (LOSAT_X_BPUNGAPSENT): this block made ready for the inlined two-hit extension of
    // the fused loop (`x_ungapped_sent.rs`); None when the switch is off or a check fails. Another
    // matrix than BLOSUM62 scores with its own byte table (`XSentQuery::new_of`, item SC2c).
    // EXPERIMENT (LOSAT_X_BPQBLOCK): none where the preliminary stage runs in blocks (`x_blocks`):
    // each block's call makes its own.
    let x_ungap_sent = if x_blocks.is_some() {
        None
    } else if args.scoring.matrix == ScoringMatrix::Blosum62 {
        x_ungapped_sent::XSentQuery::new(
            concatenated_query,
            &x_dropoff_per_context,
            word_size,
            args.word_size,
        )
    } else {
        x_ungapped_sent::XSentQuery::new_of(
            args.scoring.matrix,
            concatenated_query,
            &x_dropoff_per_context,
            word_size,
            args.word_size,
        )
    };
    if x_ungap_sent.is_none() && x_blocks.is_none() && x_ungapped_sent::on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BPUNGAPSENT", false);
    }
    let offset_array_size = blastp_offset_array_size(lookup.longest_chain());
    let query_length = blastp_concatenated_query_length(contexts);
    crate::blastinput::app::check_diag_table_window(query_length, window, "BLASTP")?;
    let (diag_array_size, diag_mask) = compute_diag_table_shape(query_length, window);
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_traceback.c:369-372
    // ```c
    // tree = Blast_IntervalTreeInit(0, query_blk->length + 1,
    //                               0, (subject_length > 0 ? subject_length :
    //                               subject_blk->length / CODON_LENGTH) + 1);
    // ```
    let total_query_span = query_length + 1;
    // `-lcase_masking`: the scan ranges of each subject of a FASTA set (`-db` subjects have no
    // case; `database` masks are other options). Empty for a subject without lower case.
    // EXPERIMENT (LOSAT_X_BPQBLOCK): made once by the search of all queries for its blocks of the
    // preliminary stage, which take them (`PrelimOnly::Block`).
    let x_soft_subjects: Option<&[Vec<(i32, i32)>]> = match x_block_inputs {
        Some((_, _, soft_subjects)) => soft_subjects,
        None => {
            *x_owned_soft_subjects = args
                .lcase_masking
                .then(|| x_soft_subject::x_subject_scan_ranges(subject_records));
            x_owned_soft_subjects.as_deref()
        }
    };
    Ok(BlastpScanShape {
        x_dropoff_per_context,
        window,
        word_size,
        x_ungap_sent,
        offset_array_size,
        diag_array_size,
        diag_mask,
        total_query_span,
        x_soft_subjects,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1397-1440 (quoted at
    // `ungapped_format_refusal`): the text of the refusal of the Seq-align formats of an ungapped
    // search.
    #[test]
    fn the_seqalign_formats_of_an_ungapped_search_are_refused_by_name() {
        assert_eq!(
            ungapped_format_refusal(8),
            "output format 8 with -ungapped (an ungapped search) is not supported by LOSAT's BLASTP"
        );
    }
}
