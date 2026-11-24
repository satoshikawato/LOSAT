//! Main BLASTN run function
//!
//! This module contains the main `run()` function that coordinates the BLASTN
//! search process.

use crate::api::local_blast::{FormatOutput, FormatProbe, OutputSink, ReportOutputs};
use crate::common::GapEditOp;
use crate::utils::glibc_fma_libm::GlibcLibm;
use crate::wording::{ENGINE_ERROR, QUERY_OPTIONS_ERROR};
use anyhow::Result;
use std::sync::Arc;

use crate::config::NuclScoringSpec;
use crate::core::blast_encoding::{
    encode_iupac_to_blastna, encode_iupac_to_blastna_into, encode_subject_ncbi2na_packed,
    COMPRESSION_RATIO,
};

// Import from parent blastn module (super::super:: because we're in blast_engine/)
use super::super::alignment::{
    blast_get_offsets_for_gapped_alignment,
    blast_get_start_for_gapped_alignment_nucl,
    build_blastna_matrix,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_scratch,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80 (BlastGapAlignStruct)
    extend_gapped_heuristic_with_traceback_with_scratch,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2762-2936 (BLAST_GreedyGappedAlignment)
    // (the score-only call goes through `x_greedy_prelim`, on the one-byte copy of the subject)
    greedy_gapped_alignment_with_traceback,
    stats_from_edit_ops,
    // NCBI's `query` of a context and `adjusted_subject` of a traceback (c++/src/algo/blast/core/blast_traceback.c:385-389,
    // 466-469)
    AdjustedSubject,
    ContextQuery,
    GapAlignScratch,
    GreedyAlignScratch,
};
// No NCBI counterpart: dispatch of the preliminary greedy alignment
// (NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4019-4025, BLAST_GreedyGappedAlignment with compressed_subject).
use super::super::alignment::x_fence;
use super::super::alignment::x_greedy_prelim;
use super::super::archive;
use super::super::args::{BlastnArgs, QueryStrand};
use super::super::blast_extend::DiagStruct;
use super::super::coordination::{
    chunk_query_masks, configure_task, finalize_task_config, prepare_sequence_data, query_masks,
    reverse_mask_intervals, subject_metadata_from_records, LookupTables, MaskError, QueryLetters,
    SequenceData, SubjectMetadata,
};
use super::super::disc_lookup::{choose_disc_scan_subject, disc_word_scan_subject};
use super::super::extension::{
    build_compressed_query, build_nucl_score_table, extend_hit_ungapped_approx_ncbi,
    extend_hit_ungapped_exact_ncbi, type_of_word, SmallNaWord, UngappedData,
};
use super::super::filtering::{
    blast_hsp_test_identity_and_length, purge_hsps_with_common_endpoints,
    purge_hsps_with_common_endpoints_ex, reevaluate_hsp_with_ambiguities_gapped_ex,
    subject_best_hit, ReevalParams,
};
use super::super::hsp::{
    adjust_blastn_offsets, evalue_comp, get_prelim_hitlist_size, sort_hsplist_by_evalue,
    sort_hsps_by_score, trim_by_max_hsps, update_best_evalue, BlastnHitList, BlastnHsp,
    BlastnHspList, BlastnOutputFormat, HitList, HitListEntry,
};
use super::super::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use super::super::lookup::{
    build_unmasked_ranges, encode_reverse_complement_blastna_into, DiscTemplates,
};
use super::super::lookup::{x_dc_csr, x_mb_lut_keep};
use super::super::pairwise::{pairwise_hits, pairwise_windows, DisplayMasks, Rows as PairwiseRows};
use super::super::query_split::{
    masked_before_part, restrict_masks, restrict_masks_with_letter_before, split_query_batch,
    splits, QueryChunk, SplitSizes,
};
use super::super::repeat_filter::MasksOutside;
use super::super::scoring::{
    check_greedy_gap_costs, check_losat_limits, check_scoring_options, ContextKarlin,
};
use super::super::seqalign::write_blastn_seqalign;
use super::super::tracing as blastn_trace;
use super::super::ungapped;
use crate::blastinput::app::{self as app, Program};
use crate::core::blast_hits::reap_by_query_coverage;
use crate::report::record::QueryNames;
// EXPERIMENT (LOSAT_X_BNSTORE): the subjects' letters read from the records or from packed subjects.
use super::super::x_subject_store::{
    x_blastna_into_pool, x_bnstore, x_read_subject_stores, SubjectLetters,
};
use crate::blastinput::fasta_reader::{
    read_queries, read_subjects, FastaInputSource, FastaRecord, QueryEnd, QueryRecords,
    ReaderConfig,
};
use crate::blastinput::input_files::{
    check_utf8_file_name, open_input, write_empty_subject_warnings_by_length,
};
use crate::blastinput::query_batch::{BatchCount, BatchSizeMixer, QueryReading, NO_DATA_MESSAGE};
use crate::blastinput::seq_range::{
    cut_queries, cut_subjects, parse_optional_range, subjects_read, Placements, QueryInput,
    RangeRole, SequenceRange,
};
use crate::report::diag_log::Posts;
use crate::report::pairwise::{
    write_blastn_pairwise_prolog, write_blastn_pairwise_report, BlastnIntro, BlastnPairwiseQuery,
    BlastnPairwiseReport, PairwiseConfig, PairwiseDatabase, PairwiseHit,
};
use crate::report::query_warnings::{query_warning, QueryReads};
use crate::sequence::packed_subject::PackedNucleotide;
use crate::stats::KarlinParams;
use crate::utils::dust::MaskedInterval;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1676
// EXPERIMENT (PAIRPAR, graduated on 2026-10-22: no variable): the one-hit seed stage of one subject
// chunk scanned in strips by the pool threads and folded in scan order by the owner (a child
// module, so that it uses the private types and kernels of this file without changing them).
#[path = "x_pair_seed.rs"]
mod x_pair_seed;

// EXPERIMENT (LOSAT_X_BNFEWSUBJ): a `-subject` search whose longest 5 Mb unit is a large part of
// its letters searches its subjects, and a batch of chunks the same way its chunks, one at a time
// on the whole pool, so that the seed stage scans each chunk in strips (`x_pair_seed`).
#[path = "x_few_subj.rs"]
mod x_few_subj;

// EXPERIMENT (LOSAT_X_BNSTREAM): the subjects of a subject set searched by the pool threads in
// OID order with a bounded window of results that the owner hands over in OID order.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
#[path = "x_bnstream.rs"]
mod x_bnstream;

// LOSAT_X_BNLUTDROP (graduated 2026-10-27): a query batch's lookup tables dropped when its preliminary stage
// ends, before its traceback, as NCBI's traceback stage drops its table.
#[path = "x_lut_drop.rs"]
mod x_lut_drop;

// EXPERIMENT (LOSAT_X_BNQPACK): the query records read before the search held packed in a store,
// each batch's letters made again when the loop takes it and dropped once the batch has encoded
// them, where no report reads them.
#[path = "x_query_drop.rs"]
mod x_query_drop;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:310-318
// EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of each query batch written when the batch loop of
// `run_in_pool` takes the batch, and the epilogs after the loop (a child module, so that it reads
// the private types of this file).
#[path = "x_rpt_stream.rs"]
mod x_rpt_stream;

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
// The search of an NCBI BLAST database (`-db`): its reading and its drivers of the
// per-subject kernel of this file (a child module, so that it uses the kernel's private
// types).
#[path = "database.rs"]
mod database;
// `-use_index true`: NCBI's loader of the MegaBLAST database index, its error lines and its
// warning, the index that the search uses, and LOSAT's refusals of the indexed searches it does
// not run (line RG item D5b).
#[path = "mb_index.rs"]
mod mb_index;
// The search of a query batch with a MegaBLAST database index (`-use_index true`): the index's
// seeds and the indexed word finder (line RG item D5b-2).
#[path = "mb_index_seeds.rs"]
mod mb_index_seeds;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_culling.c (whole file)
// BLASTN `-culling_limit`: the culling writer of the preliminary stage and the culling pipe of
// the traceback (a child module, so that it uses the kernel's private types).
#[path = "culling.rs"]
mod culling;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_besthit.c (whole file)
// BLASTN `-best_hit_overhang` and `-best_hit_score_edge`: the best hit writer of the
// preliminary stage and the best hit pipe of the traceback (a child module, so that it uses the
// kernel's private types).
#[path = "besthit.rs"]
mod besthit;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-300
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603 (the query buffer
// of a batch) and c++/src/algo/blast/core/blast_stat.c:2754-2784 (the ungapped Karlin block of
// each context).
// EXPERIMENT (LOSAT_X_QPREPPAR): both steps of `search_query_batch` as ordered parallel maps on
// the batch's search pool (a child module, so that it reads the private serial function and the
// sentinel of this file).
#[path = "x_query_prep_par.rs"]
mod x_query_prep_par;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3850 and
// c++/src/algo/blast/core/blast_traceback.c:375 (the two loops over one subject's HSPs).
// EXPERIMENT (LOSAT_X_QGROUPPAR): both loops of one subject run per query on the batch's search
// pool (a child module, beside the other modules of `search_query_batch`).
#[path = "x_qgroup_par.rs"]
mod x_qgroup_par;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1052-1174 (`s_BlastNaExtend`),
// :1486-1593 (`s_BlastSmallNaExtend`), :631-775 (`s_BlastnDiagTableExtendInitialHit`) and :799-944
// (`s_BlastnDiagHashExtendInitialHit`).
// Steps of one offset pair of a chunk's scan: the word extension of a two-stage table's hit and
// the off-diagonal search of a single word in two-hit mode (a child module, so that it uses the
// private types of this file).
#[path = "pair_extend.rs"]
mod pair_extend;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684 (the loop of
// `BlastNaWordFinder`) with `s_BlastNaExtendDirect` (:963-1030) and the scanners of
// c++/src/algo/blast/core/blast_nascan.c.
// The scan of a chunk with a table that is not two-stage and the steps of each seed it finds
// (a child module, so that it uses the private types of this file).
#[path = "na_scan.rs"]
mod na_scan;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1684 (the loop of
// `BlastNaWordFinder`) with a megablast scanner of c++/src/algo/blast/core/blast_nascan.c
// (`s_MBChooseScanSubject`, :2602) and `s_BlastNaExtend` (na_ungapped.c:1052-1174).
// The scan of a chunk with a two-stage table, which collects the offset pairs and hands each
// batch to the steps of the pairs (a child module, so that it uses the private types of this
// file).
#[path = "two_stage_scan.rs"]
mod two_stage_scan;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:701-918
// (`s_BlastSearchEngineCore`, whose :870-900 set the e-values of the subject's list and reap it),
// :422-602 (`s_BlastSearchEngineOneContext`: the chunk loop, :478-593, and one chunk's word finder
// and HSP list, :487-563) and :1480-1542 (the finish of an ungapped search's list in
// `BLAST_PreliminarySearchEngine`).
// The preliminary search of one subject: its set-up, the chunk loop with its merge, the search of
// one chunk, then the reap by e-value or the `-ungapped` finish (a child module, so that it uses
// the private types of this file).
#[path = "prelim_subject.rs"]
mod prelim_subject;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:503-578 (the chunk loop of
// `s_BlastSearchEngineOneContext` after the word finder),
// c++/src/algo/blast/core/blast_gapalign.c:3664-4114 (`BLAST_GetGappedScore`; :3850 is the loop
// on the initial hits) and :4719-4775 (`BLAST_GetUngappedHSPList`).
// The preliminary HSP list of one subject chunk from its initial hits: the gapped search's loop,
// its dispatch and the purge and sort after it, or the ungapped search's list (a child module, so
// that it uses the private types of this file).
#[path = "prelim_hsp_list.rs"]
mod prelim_hsp_list;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:635-692 (the purges, the
// re-evaluation, the score sort and the final interval tree at the end of
// `Blast_TracebackFromHSPList`) and :199-253 (`s_HSPListPostTracebackUpdate`).
// The last two steps of the traceback of one subject (a child module, so that it uses the
// private types of this file).
#[path = "traceback_finish.rs"]
mod traceback_finish;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1521-1713 (one subject of
// the loop of `BLAST_ComputeTraceback_MT`: the subject it fetches, :1526-1528,1556-1573, the
// traceback of its lists, :1632-1708, and the search without a traceback, :1691-1696).
// The traceback of one subject: its letters, then the traceback of its kept HSPs or that of an
// ungapped search (a child module, so that it uses the private types of this file).
#[path = "traceback_subject.rs"]
mod traceback_subject;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427 (the loop of
// `BLAST_PreliminarySearchEngine` on the subjects) and c++/src/algo/blast/core/blast_traceback.c:1521
// (the loop of `BLAST_ComputeTraceback_MT` on the subjects' batches).
// The search of one subject by stage, the passes over a batch's subjects (a subject set's in subject
// order, a database's by blocks), the preliminary pass and the traceback pass (a child module, so that
// it uses the private types of this file).
#[path = "subject_passes.rs"]
mod subject_passes;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88 (the thread count of a
// preliminary search).
// The batch's clocks, thread request, task configuration and speculative traceback rules (a child
// module, so that it uses the private types of this file).
#[path = "batch_schedule.rs"]
mod batch_schedule;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:375-605 (the loop of
// `Blast_TracebackFromHSPList`, from the start of one HSP's traceback to its insertion in the tree).
// Where the traceback of one HSP starts, its traceback call, the speculative batch of the
// traceback loop, and the alignment of one HSP in the loop and the HSP it makes (a child module,
// so that it uses the private types of this file).
#[path = "traceback_hsp.rs"]
mod traceback_hsp;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:259-723
// (`Blast_TracebackFromHSPList`; :365-375,611 are the sort assertion, the interval tree and the
// loop on the HSPs).
// The traceback of one subject's kept HSPs and its loop over one list of them (a child module,
// so that it uses the private types of this file).
#[path = "traceback_list.rs"]
mod traceback_list;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
// EXPERIMENT (LOSAT_X_QSTREAM): the records of each query batch read from the query file when the
// batch loop of `run_in_pool` needs them, after a first reading that keeps only each record's
// length (a child module, so that it reads the private types of this file).
#[path = "x_query_stream.rs"]
mod x_query_stream;

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279 (the batch reading that
// LOSAT_X_QSTREAM's first reading precedes).
// That first reading cut into parts at record starts and read by the search's threads (a child
// module, beside `x_query_stream`, which calls it).
#[path = "x_query_scan_par.rs"]
mod x_query_scan_par;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:152-273,484-658 (the query
// block), c++/src/algo/blast/core/blast_nalookup.c:1227-1368 (the lookup table) and
// c++/src/algo/blast/core/blast_parameters.c:279-419,830-994 (the cutoffs).
// The query block, the lookup table and the cutoffs of a query batch (a child module, so that it
// uses the private types of this file).
#[path = "batch_setup.rs"]
mod batch_setup;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:455-532
// (`BlastSetup_ScoreBlkInit`), :887-999 (`BLAST_GapAlignSetUp`) and :699-850
// (`BLAST_CalcEffLengths`).
// The Karlin-Altschul blocks, the checks of the set-up and the effective search spaces of a query
// batch's contexts (a child module, so that it uses the private types of this file).
#[path = "context_statistics.rs"]
mod context_statistics;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607 (the set-up of
// `BLAST_PreliminarySearchEngine`), c++/src/algo/blast/core/blast_parameters.c:158-277,422-490
// (the X-dropoffs) and c++/src/algo/blast/api/prelim_stage.cpp:141-203 (the threads).
// The values of a query batch's subject search and the schedule of its subjects (a child module,
// so that it uses the private types and modules of this file).
#[path = "subject_env.rs"]
mod subject_env;

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:708-755 (the XML, XML2 and SAM
// branches of `CBlastFormat::x_PrintStructuredReport`); the refusals of the query-anchored reports have
// no NCBI counterpart.
// Parts of `post_process_hits_and_write`, the refusals of outfmt 1 to 4, the SAM writer and the
// structured documents (outfmt 5, 13 to 16), and the subjects' names that the Seq-align, archive and
// SAM writers take (a child module, so that it uses the private types of this file).
#[path = "report_write.rs"]
mod report_write;

/// EXPERIMENT (LOSAT_X_QPREPPAR): the summary line `main` prints at exit. No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_query_prep_par_print_summary() {
    x_query_prep_par::print_summary();
}

/// EXPERIMENT (LOSAT_X_QGROUPPAR, LOSAT_X_STATS): the summary lines `main` prints at exit.
/// No NCBI counterpart: it prints counters and changes no value NCBI computes.
pub fn x_qgroup_par_print_summary() {
    x_qgroup_par::print_summary();
}

/// EXPERIMENT (LOSAT_X_BNTBRANGELOG): the summary line `main` prints at exit. No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_tb_range_print_summary() {
    database::x_tb_range_print_summary();
}

/// LOSAT_X_BNLUTDROP (graduated 2026-10-27): the summary line `main` prints at exit (LOSAT_TIMING). No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_lut_drop_print_summary() {
    x_lut_drop::print_summary();
}

/// EXPERIMENT (LOSAT_X_BNQPACK): the summary line `main` prints at exit (LOSAT_TIMING). No NCBI
/// counterpart: it prints counters and changes no value NCBI computes.
pub fn x_query_drop_print_summary() {
    x_query_drop::print_summary();
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/ncbi_math.h:160-161
const NCBIMATH_LN2: f64 = std::f64::consts::LN_2; // the same f64 as NCBI's literal

/// The text of a query's outfmt 0 `Query=` line and outfmt 7 `# Query:` line: the record's
/// title bytes (the local ID is not shown without `-parse_deflines`; the writers apply
/// `TruncateSpaces`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:732-734
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:633-635
#[inline]
fn record_title(record: &FastaRecord) -> Arc<[u8]> {
    Arc::from(record.title.as_slice())
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/gapinfo.h:43-61
fn format_gap_edit_ops_for_trace(ops: &[GapEditOp]) -> String {
    let mut out = String::from("[");
    for (idx, op) in ops.iter().enumerate() {
        if idx > 0 {
            out.push(',');
        }
        match *op {
            GapEditOp::Sub(n) => {
                out.push('S');
                out.push_str(&n.to_string());
            }
            GapEditOp::Del(n) => {
                out.push('D');
                out.push_str(&n.to_string());
            }
            GapEditOp::Ins(n) => {
                out.push('I');
                out.push_str(&n.to_string());
            }
        }
    }
    out.push(']');
    out
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1887-1890
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1923-1928
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:4157-4171
fn calculate_blastn_context_statistics(
    raw_score: i32,
    params: &crate::stats::KarlinParams,
    eff_searchsp: i64,
    round_down_evalue_score: bool,
) -> (f64, f64) {
    let log_k = params.k.glibc_ln();
    let bit_score = (params.lambda * (raw_score as f64) - log_k) / NCBIMATH_LN2;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1864-1870
    let evalue_score = if round_down_evalue_score {
        raw_score & !1
    } else {
        raw_score
    };
    let e_value = if params.lambda < 0.0 || params.k < 0.0 || params.h < 0.0 {
        -1.0
    } else {
        (eff_searchsp as f64) * (-(params.lambda) * (evalue_score as f64) + log_k).glibc_exp()
    };
    (bit_score, e_value)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:662-664
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1996-1998
// The HSP survives unless `evalue > cutoff`, the same comparison as C (a NaN
// e-value survives; LOSAT rejects a NaN `-evalue` before the search).
#[allow(clippy::neg_cmp_op_on_partial_ord)]
fn hsp_survives_evalue_reap(evalue: f64, cutoff: f64) -> bool {
    !(evalue > cutoff)
}

/// The reap of a subject's preliminary HSP list after its e-values are computed (and kept in the
/// HSPs), chosen once for the batch: by the e-value threshold, or, under matrix-only scoring, by
/// the raw score.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:891-899 (the same at
/// 1527-1535 for an ungapped search's finish)
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) enum PrelimReap {
    /// `s_Blast_HSPListReapByPrelimEvalue`: an HSP whose e-value is above the threshold
    /// (`hit_params->prelim_evalue`, the e-value option) is dropped (`hsp_survives_evalue_reap`).
    Evalue(f64),
    /// `Blast_HSPListReapByRawScore`: an HSP whose raw score is below the hit-saving option's
    /// cutoff score (`-min_raw_gapped_score`, 0 when not given) is dropped.
    RawScore(i32),
}

impl PrelimReap {
    /// The reap of a search: by raw score under matrix-only scoring, else by e-value.
    ///
    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2087-2097
    /// (`Blast_HSPListReapByRawScore`)
    /// `min_raw_gapped_score` is `TaskConfig::min_raw_gapped_score` (the option above 0): a cutoff
    /// of 0 or less keeps every HSP, whose scores are positive.
    pub(super) fn new(
        matrix_only: bool,
        min_raw_gapped_score: Option<i32>,
        evalue_threshold: f64,
    ) -> Self {
        if matrix_only {
            PrelimReap::RawScore(min_raw_gapped_score.unwrap_or(0))
        } else {
            PrelimReap::Evalue(evalue_threshold)
        }
    }

    /// Each HSP's e-value (`evalue_of`) kept in it, and the HSPs that the reap drops removed, in
    /// place and in order. The choice is made outside the loop over the HSPs.
    #[inline]
    fn reap(self, hits: &mut Vec<PrelimHit>, evalue_of: impl Fn(&PrelimHit) -> f64) {
        match self {
            PrelimReap::Evalue(threshold) => hits.retain_mut(|hit| {
                hit.prelim_evalue = evalue_of(hit);
                hsp_survives_evalue_reap(hit.prelim_evalue, threshold)
            }),
            PrelimReap::RawScore(cutoff) => hits.retain_mut(|hit| {
                hit.prelim_evalue = evalue_of(hit);
                hit.prelim_score >= cutoff
            }),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4155-4160
const MAX_SUBJECT_OFFSET: i32 = 90000;
const MAX_TOTAL_GAPS: i32 = 3000;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:54
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:192
// The constants of NCBI's chunks and the chunks of hard ranges have one owner,
// `core::subject_split`, which the refusals of `-db_hard_mask` read too.
use crate::core::subject_split::{HardSplitState, DBSEQ_CHUNK_OVERLAP, MAX_DBSEQ_LEN};
// NCBI reference (598d8ae6): c++/include/algo/blast/core/lookup_wrap.h:119
const OFFSET_ARRAY_SIZE: usize = 4096;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:120-121
const NUCL_SENTINEL: u8 = 0x0F;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1534-1539
const OVERLAP_DIAG_CLOSE: i32 = 10;

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:50-54
const DIAGHASH_NUM_BUCKETS: usize = 512;
const DIAGHASH_CHAIN_LENGTH: usize = 256;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:122-143
struct SubjectSplitState {
    full_right: i32,
    hard_ranges: [(i32, i32); 1],
    soft_ranges: Vec<(i32, i32)>,
    hm_index: usize,
    sm_index: usize,
    offset: i32,
    next: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-276
#[derive(Clone)]
enum SeqRanges {
    Inline { range: (i32, i32) },
    SoftRangesFull,
    SoftRangesSlice { start_idx: usize, count: usize },
}

impl SeqRanges {
    /// The chunk's ranges: the one range, the subject's soft ranges, or a cut of them
    /// (`cut_soft_ranges`). The first two are what every chunk of an unmasked subject reads, so
    /// they are taken where the chunk is searched; the cut stays a call.
    #[inline(always)]
    fn as_slice<'a>(
        &'a self,
        soft_ranges: &'a [(i32, i32)],
        scratch: &'a mut Vec<(i32, i32)>,
        chunk_offset: i32,
        chunk_length: i32,
    ) -> &'a [(i32, i32)] {
        match self {
            SeqRanges::Inline { range } => std::slice::from_ref(range),
            SeqRanges::SoftRangesFull => soft_ranges,
            SeqRanges::SoftRangesSlice { start_idx, count } => Self::cut_soft_ranges(
                soft_ranges,
                scratch,
                *start_idx,
                *count,
                chunk_offset,
                chunk_length,
            ),
        }
    }

    /// The soft ranges `start_idx..start_idx + count` of the subject in the chunk's coordinates,
    /// clipped to the chunk, in `scratch`.
    #[inline(never)]
    fn cut_soft_ranges<'a>(
        soft_ranges: &[(i32, i32)],
        scratch: &'a mut Vec<(i32, i32)>,
        start_idx: usize,
        count: usize,
        chunk_offset: i32,
        chunk_length: i32,
    ) -> &'a [(i32, i32)] {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:184-198
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:298-308
        let start = start_idx;
        scratch.clear();
        if scratch.capacity() < count {
            scratch.reserve(count.saturating_sub(scratch.capacity()));
        }
        let end_idx = start.saturating_add(count);
        for &(left, right) in soft_ranges[start..end_idx].iter() {
            scratch.push((left - chunk_offset, right - chunk_offset));
        }
        if let Some(first) = scratch.first_mut() {
            if first.0 < 0 {
                first.0 = 0;
            }
        }
        if let Some(last) = scratch.last_mut() {
            if last.1 > chunk_length {
                last.1 = chunk_length;
            }
        }
        scratch.as_slice()
    }
}

#[derive(Clone)]
struct SubjectChunk {
    offset: usize,
    length: usize,
    seq_ranges: SeqRanges,
    masked: bool,
    overlap: usize,
}

enum SubjectChunkStatus {
    Done,
    NoRange,
    Ok(SubjectChunk),
}

impl SubjectSplitState {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:146-184
    fn new(subject_len: usize, soft_ranges: Vec<(i32, i32)>) -> Self {
        let full_right = subject_len as i32;
        Self {
            full_right,
            // One range, as NCBI's `&(backup->full_range)`: an array, not a vector per subject.
            hard_ranges: [(0, full_right)],
            soft_ranges,
            hm_index: 0,
            sm_index: 0,
            offset: 0,
            next: 0,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:221-310 (the arithmetic of
    // `core::subject_split::HardSplitState::next_chunk` on `&full_range`, kept here for the default
    // path's code: one function for both cost the chunk loop 14 instructions per subject, review 2
    // S1 of M4a-2b)
    fn next_chunk(&mut self, subject_masked: bool, chunk_overlap: usize) -> SubjectChunkStatus {
        if self.next >= self.full_right {
            return SubjectChunkStatus::Done;
        }

        let residual = (self.next as usize) % COMPRESSION_RATIO;
        self.offset = self.next - residual as i32;
        let offset = self.offset;

        let hard_right = self
            .hard_ranges
            .get(self.hm_index)
            .map(|r| r.1)
            .unwrap_or(self.full_right);

        let length = if offset + (MAX_DBSEQ_LEN as i32) < hard_right {
            self.next = offset + MAX_DBSEQ_LEN as i32 - chunk_overlap as i32;
            MAX_DBSEQ_LEN
        } else {
            let len = (hard_right - offset).max(0) as usize;
            self.hm_index = self.hm_index.saturating_add(1);
            self.next = if self.hm_index < self.hard_ranges.len() {
                self.hard_ranges[self.hm_index].0
            } else {
                self.full_right
            };
            len
        };

        let no_chunking = offset == 0 && residual == 0 && self.next == self.full_right;
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-276
        let seq_ranges = if no_chunking {
            if subject_masked {
                SeqRanges::SoftRangesFull
            } else {
                SeqRanges::Inline {
                    range: (0, length as i32),
                }
            }
        } else if !subject_masked {
            SeqRanges::Inline {
                range: (residual as i32, length as i32),
            }
        } else {
            debug_assert_eq!(residual, 0);
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:279-291
            let start = offset;
            let end = offset + length as i32;
            let mut i = self.sm_index;
            while i < self.soft_ranges.len() && self.soft_ranges[i].1 < start {
                i += 1;
            }
            let start_idx = i;
            while i < self.soft_ranges.len() && self.soft_ranges[i].0 < end {
                i += 1;
            }
            let count = i - start_idx;
            self.sm_index = i.saturating_sub(1);
            if count == 0 {
                return SubjectChunkStatus::NoRange;
            }
            SeqRanges::SoftRangesSlice { start_idx, count }
        };

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:579-584
        let range_start = if self.hm_index < self.hard_ranges.len() {
            self.hard_ranges[self.hm_index].0
        } else {
            self.full_right
        };
        let overlap = if offset == range_start {
            0
        } else {
            chunk_overlap
        };

        SubjectChunkStatus::Ok(SubjectChunk {
            offset: offset as usize,
            length,
            seq_ranges,
            masked: subject_masked,
            overlap,
        })
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:264-276 (cited at
// `SubjectSplitState::next_chunk`): the "no chunking" branch takes the soft ranges, the full
// range of a hard-masked subject (blast_engine.c:162-164); else `[residual, length)`.
/// `SubjectSplitState::next_chunk` of a subject with hard ranges (`-db_hard_mask`, series NP,
/// M4a-2b decision H1): the next chunk that the worker's `db_hard_split` cuts (NCBI's arithmetic,
/// `core::subject_split::HardSplitState::next_chunk`, which the refusals of `HardMaskShape` read
/// too), scanned from its residual (`SeqRanges::Inline`), of `right - offset` letters: 0 for a
/// hard range of no letters (a mask at 0, masks that touch), whose `Blast_ExtendWordExit` still
/// moves the diagonal offset on by the window (`collect_prelim_hits_for_chunk`). NCBI's "no
/// chunking" branch scans the full range instead, which is the chunk's `[0, length)` but for a
/// chunk from the start that a mask running to the end follows: LOSAT scans that chunk alone
/// (series NP, M4a-2c: the approved NCBI defect of `blastinput::db_mask::TAIL_MASK_SEARCHED_WORD`;
/// below that word `HardMaskShape::TailMask` refuses the sequence before the search). After the
/// last hard range NCBI reads the overlap one `SSeqRange` past the list (another sequence's mask
/// end, or unset memory; SURVEY §7.3 of M4a-2): LOSAT takes `dbseq_chunk_overlap` there, as the
/// unmasked search does for every sequence over `MAX_DBSEQ_LEN` letters (`(&full_range)[1]`); it
/// decides a merge across a mask only in the shape `HardMaskShape::ShortLastMask` refuses
/// (decision H5), while the final chunk of a last range over `MAX_DBSEQ_LEN` letters merges with
/// the chunk before it in the same range (decision B-d). Cold and out of line: the other
/// subjects' chunk loops never call it.
#[cold]
#[inline(never)]
fn next_hard_chunk(hard: &mut HardSplitState, chunk_overlap: usize) -> SubjectChunkStatus {
    match hard.next_chunk(chunk_overlap) {
        None => SubjectChunkStatus::Done,
        Some(chunk) => SubjectChunkStatus::Ok(SubjectChunk {
            offset: chunk.offset,
            length: chunk.letters(),
            seq_ranges: SeqRanges::Inline {
                range: (chunk.residual as i32, chunk.length),
            },
            masked: true,
            overlap: chunk.overlap,
        }),
    }
}
/// The next chunk of a subject: of its hard ranges (`HARD`: a database sequence of a search with
/// `-db_hard_mask`, whose ranges the worker's `db_hard_split` holds), else of its
/// `SubjectSplitState`. `HARD` is a constant of the chunk loop's instance
/// (`prelim_subject::collect_subject_prelim`), so the loops of the other subjects test nothing.
#[inline(always)]
fn next_subject_chunk<const HARD: bool>(
    split: &mut SubjectSplitState,
    hard: &mut HardSplitState,
    subject_masked: bool,
    chunk_overlap: usize,
) -> SubjectChunkStatus {
    if HARD {
        next_hard_chunk(hard, chunk_overlap)
    } else {
        split.next_chunk(subject_masked, chunk_overlap)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4163-4191
#[inline]
fn adjust_subject_range(
    subject_offset: &mut i32,
    subject_length: &mut i32,
    query_offset: i32,
    query_length: i32,
) -> i32 {
    let mut start_shift: i32 = 0;
    if *subject_length < MAX_SUBJECT_OFFSET {
        return start_shift;
    }

    let s_offset = *subject_offset;
    let max_extension_left = query_offset + MAX_TOTAL_GAPS;
    let max_extension_right = query_length - query_offset + MAX_TOTAL_GAPS;

    if s_offset <= max_extension_left {
        start_shift = 0;
    } else {
        start_shift = s_offset - max_extension_left;
        *subject_offset = max_extension_left;
    }

    let adjusted_end = std::cmp::min(*subject_length, s_offset + max_extension_right);
    *subject_length = adjusted_end - start_shift;
    start_shift
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-182
#[inline]
fn advance_diag_table_offset(
    diag_table_offset: &mut i32,
    diag_window: i32,
    subject_len: usize,
    use_array_indexing: bool,
    hit_level_array: &mut [DiagStruct],
    hit_len_array: &mut [u8],
    diag_hash: &mut DiagHashTable,
) {
    if use_array_indexing {
        if *diag_table_offset >= i32::MAX / 4 {
            *diag_table_offset = diag_window;
            let init_diag = DiagStruct {
                last_hit: -diag_window,
                flag: 0,
            };
            hit_level_array.fill(init_diag);
            if !hit_len_array.is_empty() {
                hit_len_array.fill(0);
            }
        } else {
            *diag_table_offset = diag_table_offset.saturating_add(subject_len as i32 + diag_window);
        }
    } else if diag_hash.offset >= i32::MAX / 4 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:174-179
        diag_hash.occupancy = 1;
        diag_hash.offset = diag_hash.window;
        diag_hash.backbone.fill(0);
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): every chain is empty now. No NCBI
        // counterpart (empty while the switch is off).
        diag_hash.max_level.fill(i32::MIN);
        // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): no bucket has an insert now. No NCBI
        // counterpart (empty while the switch is off).
        diag_hash.x_last.fill(x_diag_lv::LvLast::EMPTY);
        #[cfg(debug_assertions)]
        {
            diag_hash.x_last_fresh = true;
        }
    } else {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:180-181
        diag_hash.offset = diag_hash
            .offset
            .saturating_add(subject_len as i32 + diag_hash.window);
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:62-71
#[derive(Clone, Copy, Default)]
struct DiagHashCell {
    diag: i32,
    level: i32,
    hit_saved: bool,
    hit_len: i32,
    next: u32,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:97-106
struct DiagHashTable {
    num_buckets: usize,
    occupancy: u32,
    capacity: u32,
    backbone: Vec<u32>,
    chain: Vec<DiagHashCell>,
    offset: i32,
    window: i32,
    // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): per bucket, an upper bound of the levels of
    // its cells (the largest level `insert` wrote since the table was made or cleared;
    // `i32::MIN` for none). Empty while the switch is off. No NCBI counterpart.
    max_level: Vec<i32>,
    // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): per bucket, the latest insert and an upper
    // bound of the levels of the bucket's other cells. Empty unless the switch and
    // LOSAT_X_DIAGLV are on and the table takes two hits. No NCBI counterpart. It is the memo of
    // the cells only from a dc-megablast strips stage's `x_lv2_reset` to the end of that stage
    // (and on a new or cleared table): the inserts of the reference fold (`insert_plain`) do not
    // keep it, so only that stage's folds may read it (`retrieve_below::<true>`) or keep it
    // (`insert_memo`, `insert_inl::<true>`).
    x_last: Vec<x_diag_lv::LvLast>,
    // Debug builds only (LOSAT_X_DIAGLV2): whether `x_last` is the memo of the cells (set by
    // `new`, the clear and `x_lv2_reset`, ended by an insert that does not keep a memo the table
    // has); the memo's reads and upkeep assert it.
    #[cfg(debug_assertions)]
    x_last_fresh: bool,
}

impl DiagHashTable {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:122-134
    fn new(window: i32) -> Self {
        let capacity = DIAGHASH_CHAIN_LENGTH as u32;
        // EXPERIMENT (LOSAT_X_DIAGLV): the level bounds when the switch is on.
        let max_level = x_diag_lv::new_levels();
        Self {
            num_buckets: DIAGHASH_NUM_BUCKETS,
            occupancy: 1,
            capacity,
            backbone: vec![0; DIAGHASH_NUM_BUCKETS],
            chain: vec![DiagHashCell::default(); capacity as usize],
            offset: window,
            window,
            max_level,
            // EXPERIMENT (LOSAT_X_DIAGLV2): the memo when the switch is on.
            x_last: x_diag_lv::new_last(window),
            #[cfg(debug_assertions)]
            x_last_fresh: true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:361-381
    fn retrieve(&self, diag: i32) -> Option<(i32, i32, bool)> {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &self.chain[index as usize];
            if cell.diag == diag {
                return Some((cell.level, cell.hit_len, cell.hit_saved));
            }
            index = cell.next;
        }
        None
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:342-449
    // The insert of a dc-megablast strips stage with the memo of LOSAT_X_DIAGLV2 (`x_last`,
    // valid from the stage's `x_lv2_reset`): the statements of `insert_plain` and the memo's
    // upkeep (`insert_inl::<true>`).
    fn insert_memo(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        self.insert_inl::<true>(diag, level, hit_len, hit_saved, s_off, window_size)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
    /// (`s_BlastDiagHashInsert`, cited above `insert_memo`).
    /// `insert_memo` without the memo of LOSAT_X_DIAGLV2 (the code before that switch): the reference
    /// folds (`pair_extend::extend_offset_pairs` for a two-stage table, `na_scan::scan_na_lookup`
    /// for the other tables), which never read the memo (it lives in the strips stages of
    /// dc-megablast only, rebuilt at each stage's start: `x_diag_lv.rs`), the one-hit strips folds
    /// (`x_pair_seed.rs` `fold`, `x_pair_seed_na.rs` `fold_na`), and the strips fold of
    /// dc-megablast in its instance with LOSAT_X_FOLDLEAN off on a table without the memo.
    #[inline(never)]
    pub(super) fn insert_plain(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        self.insert_inl::<false>(diag, level, hit_len, hit_saved, s_off, window_size)
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:396-449
    /// (`s_BlastDiagHashInsert`, cited above `insert_memo`).
    /// The body of `insert_memo` (`MEMO`) and `insert_plain`, inlined into its caller; the strips fold
    /// of dc-megablast inlines it with LOSAT_X_FOLDLEAN (`x_pair_seed_disc.rs` `fold_two_hit`).
    /// The same statements. Without `MEMO` the memo (LOSAT_X_DIAGLV2) is not kept.
    #[inline(always)]
    pub(super) fn insert_inl<const MEMO: bool>(
        &mut self,
        diag: i32,
        level: i32,
        hit_len: i32,
        hit_saved: bool,
        s_off: i32,
        window_size: i32,
    ) -> bool {
        let bucket =
            ((diag as u32).wrapping_mul(0x9E370001) % DIAGHASH_NUM_BUCKETS as u32) as usize;
        // (Debug builds: the memo's upkeep needs the memo of the cells; an insert that does not
        // keep a memo the table has ends it until the next stage's `x_lv2_reset`.)
        #[cfg(debug_assertions)]
        {
            if MEMO {
                assert!(
                    self.x_last.is_empty() || self.x_last_fresh,
                    "LOSAT_X_DIAGLV2: an insert keeps a stale memo (no x_lv2_reset since an insert without it)"
                );
            } else if !self.x_last.is_empty() {
                self.x_last_fresh = false;
            }
        }
        // EXPERIMENT (LOSAT_X_DIAGLV, `x_diag_lv.rs`): this insert writes `level` into one cell
        // of this bucket (the match, a take-over or the append below). No NCBI counterpart;
        // with the switch off `max_level` is empty and this is one length test.
        if let Some(max_level) = self.max_level.get_mut(bucket) {
            if *max_level < level {
                *max_level = level;
            }
            // EXPERIMENT (LOSAT_X_DIAGLV2, `x_diag_lv.rs`): this insert is the bucket's latest.
            // No NCBI counterpart; with the switch off `x_last` is empty and this is one length
            // test (none without `MEMO`: the reference fold's `insert_plain`).
            if MEMO {
                if let Some(last) = self.x_last.get_mut(bucket) {
                    last.note(diag, level, hit_len, hit_saved);
                }
            }
        }
        let mut index = self.backbone[bucket];
        while index != 0 {
            let cell = &mut self.chain[index as usize];
            if cell.diag == diag {
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            if s_off - cell.level > window_size {
                cell.diag = diag;
                cell.level = level;
                cell.hit_len = hit_len;
                cell.hit_saved = hit_saved;
                return true;
            }
            index = cell.next;
        }

        if self.occupancy == self.capacity {
            self.capacity = self.capacity.saturating_mul(2);
            self.chain
                .resize(self.capacity as usize, DiagHashCell::default());
        }

        let cell_index = self.occupancy;
        let cell = &mut self.chain[cell_index as usize];
        cell.diag = diag;
        cell.level = level;
        cell.hit_len = hit_len;
        cell.hit_saved = hit_saved;
        cell.next = self.backbone[bucket];
        self.backbone[bucket] = cell_index;
        self.occupancy = self.occupancy.saturating_add(1);
        true
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:833-839
// EXPERIMENT (LOSAT_X_DIAGLV): a fold's retrieve skips the walk of
// a bucket whose cells all lie below the level that could change a decision of the seed (a
// child module, so that it uses the private table of this file).
#[path = "x_diag_lv.rs"]
mod x_diag_lv;

#[inline(always)]
fn diag_hash_insert_window(window_size: usize, scan_range: usize, word_length: usize) -> i32 {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:825-858, 939-941
    // The clamp at line 858 is inside the two-hit off-diagonal block. In BLASTN
    // one-hit mode (window_size == 0) that block is skipped, so the insert sees
    // the original negative Delta.
    let window = window_size as i32;
    let delta = (scan_range as i32).min(window - word_length as i32);
    window + delta + 1
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:135-150
#[derive(Clone, Copy)]
struct OffsetPair {
    q_off: usize,
    s_off: usize,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:630-637
    // Each seed must carry the current subject range boundary used by the
    // original NCBI extension call.
    s_range: usize,
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_util.h:51-55
#[inline(always)]
fn packed_base_at(packed: &[u8], pos: usize) -> u8 {
    let byte = packed[pos / COMPRESSION_RATIO];
    let shift = 2 * (3 - (pos % COMPRESSION_RATIO));
    (byte >> shift) & 0x03
}

/// No NCBI counterpart (bounds checks only): the `N` bytes of `packed` from `i` read as one
/// big-endian number (`s[0] << 24 | s[1] << 16 | s[2] << 8 | s[3]` of the scans for `N` = 4,
/// `s[0] << 16 | s[1] << 8 | s[2]` for 3), or `None` when they pass its end: one test of the
/// length where indexing each byte tested it `N` times.
#[inline(always)]
fn x_be_bytes<const N: usize>(packed: &[u8], i: usize) -> Option<u32> {
    let b = packed.get(i..)?.first_chunk::<N>()?;
    Some(b.iter().fold(0u32, |v, &x| (v << 8) | x as u32))
}

/// `x_be_bytes` at an index that the scans keep as `i32` (`None` below 0).
#[inline(always)]
fn x_be_bytes_from<const N: usize>(packed: &[u8], i: i32) -> Option<u32> {
    x_be_bytes::<N>(packed, usize::try_from(i).ok()?)
}

/// No NCBI counterpart (bounds checks only): one round of a scan's loop body, the positions it
/// takes until its reads are aligned as at its start again (`steps` positions, `bytes` bytes on),
/// and the end of its reads (they end before byte `reads_end` of its first byte).
#[derive(Clone, Copy)]
struct XRound {
    steps: usize,
    bytes: usize,
    reads_end: usize,
}

/// No NCBI counterpart (bounds checks only): the number of whole rounds (`XRound`) a scan at
/// position `pos` (`step` letters apart, the last position at most `end`) with its round's
/// first byte at `byte` runs with every read inside `packed_len` bytes. Those rounds read with
/// `x_be4` and test neither the range nor the length per position, as NCBI's loop bodies (whose
/// subject has a sentinel byte); the scan's own loop, with its tests, takes the positions before
/// the first round and after the last, so every position gets the same word in the same order.
#[inline(always)]
fn x_whole_rounds(
    round: XRound,
    step: usize,
    pos: usize,
    end: usize,
    byte: usize,
    packed_len: usize,
) -> usize {
    let by_pos = end
        .checked_sub(pos + (round.steps - 1) * step)
        .map_or(0, |d| d / (round.steps * step) + 1);
    let by_len = packed_len
        .checked_sub(byte + round.reads_end)
        .map_or(0, |d| d / round.bytes + 1);
    by_pos.min(by_len)
}

/// No NCBI counterpart (bounds checks only): `s[i] << 24 | s[i + 1] << 16 | s[i + 2] << 8 |
/// s[i + 3]` of the scans, with no test of the length (`s[i] << 16 | s[i + 1] << 8 | s[i + 2]` is
/// this `>> 8`).
///
/// # Safety
/// `i + 4 <= packed.len()` (`x_whole_rounds`).
#[inline(always)]
unsafe fn x_be4(packed: &[u8], i: usize) -> u32 {
    debug_assert!(i + 4 <= packed.len());
    // SAFETY: the caller's: the 4 bytes from `i` are in `packed`.
    u32::from_be_bytes(unsafe { packed.as_ptr().add(i).cast::<[u8; 4]>().read_unaligned() })
}

/// No NCBI counterpart (bounds checks only): `packed[i]` with no test of the length.
///
/// # Safety
/// `i < packed.len()` (`x_whole_rounds`).
#[inline(always)]
unsafe fn x_byte(packed: &[u8], i: usize) -> u32 {
    debug_assert!(i < packed.len());
    // SAFETY: the caller's.
    u32::from(unsafe { *packed.get_unchecked(i) })
}

#[inline]
fn packed_kmer_at(packed: &[u8], start: usize, k: usize) -> u64 {
    let mut kmer = 0u64;
    for i in 0..k {
        let code = packed_base_at(packed, start + i) as u64;
        kmer = (kmer << 2) | code;
    }
    kmer
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:459-487
#[inline(always)]
fn packed_kmer_at_seed_mask(packed: &[u8], start: usize, lut_word_length: usize) -> u64 {
    if lut_word_length == 0 {
        return 0;
    }

    let shift = 2i32 * (16i32 - (start % COMPRESSION_RATIO) as i32 - lut_word_length as i32);
    if !(0..=24).contains(&shift) {
        return packed_kmer_at(packed, start, lut_word_length);
    }

    let byte_idx = start / COMPRESSION_RATIO;
    let shift_u = shift as u32;

    match shift {
        8 | 10 | 12 | 14 => {
            if byte_idx + 3 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                let s2 = packed[byte_idx + 2] as u32;
                (((s0 << 24) | (s1 << 16) | (s2 << 8)) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        16 | 18 | 20 | 22 => {
            if byte_idx + 2 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                (((s0 << 24) | (s1 << 16)) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        24 => {
            if byte_idx < packed.len() {
                packed[byte_idx] as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
        _ => {
            if byte_idx + 4 <= packed.len() {
                let s0 = packed[byte_idx] as u32;
                let s1 = packed[byte_idx + 1] as u32;
                let s2 = packed[byte_idx + 2] as u32;
                let s3 = packed[byte_idx + 3] as u32;
                (((s0 << 24) | (s1 << 16) | (s2 << 8) | s3) >> shift_u) as u64
            } else {
                packed_kmer_at(packed, start, lut_word_length)
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1482-1623
fn scan_subject_kmers_range_mb_any<F>(
    packed: &[u8],
    subject_len: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    if subject_len < lut_word_length || !(9..=12).contains(&lut_word_length) {
        return;
    }

    let packed_len = packed.len();
    if packed_len == 0 {
        return;
    }

    let mask = (1u64 << (2 * lut_word_length)) - 1;

    if scan_step.is_multiple_of(COMPRESSION_RATIO) && !subject_masked && packed_len >= 3 {
        let shift = 2 * (12 - lut_word_length);
        let step_bytes = scan_step / COMPRESSION_RATIO;
        let mut byte_idx = start / COMPRESSION_RATIO;
        let max_byte_idx = packed_len.saturating_sub(3);
        let end_byte_idx = (end / COMPRESSION_RATIO).min(max_byte_idx);
        // The rounds (`x_whole_rounds`): one word per round, its read aligned as the first.
        let round = XRound {
            steps: 1,
            bytes: step_bytes,
            reads_end: 4,
        };
        let first = byte_idx * COMPRESSION_RATIO;
        let end_pos = end_byte_idx * COMPRESSION_RATIO;
        // (the three bytes are the four-byte read `>> 8`: one shift)
        let shift4 = shift + 8;
        if byte_idx <= end_byte_idx {
            for _ in 0..x_whole_rounds(round, scan_step, first, end_pos, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's read ends before `byte_idx + 4`.
                let idx = unsafe { x_be4(packed, byte_idx) };
                if !on_kmer(
                    byte_idx * COMPRESSION_RATIO,
                    ((idx >> shift4) as u64) & mask,
                ) {
                    return;
                }
                byte_idx += step_bytes;
            }
        }
        while byte_idx <= end_byte_idx {
            let Some(idx) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            let kmer = ((idx >> shift) as u64) & mask;
            if !on_kmer(byte_idx * COMPRESSION_RATIO, kmer) {
                return;
            }
            byte_idx = byte_idx.saturating_add(step_bytes);
        }
        let mut pos = byte_idx * COMPRESSION_RATIO;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, lut_word_length);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    if lut_word_length > 9 {
        if packed_len >= 4 {
            let max_byte_idx = packed_len.saturating_sub(4);
            let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
            let fast_end = end.min(max_fast_pos);
            let mut pos = start;
            // The words up to `fast_end`, counted once: each word's read at its byte, at most
            // `max_byte_idx` there, with no test of the length (as `x_whole_rounds`).
            if pos <= fast_end {
                for _ in 0..(fast_end - pos) / scan_step + 1 {
                    let shift = 2 * (16 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                    // SAFETY: `pos <= fast_end`, so `pos / 4 <= max_byte_idx = packed_len - 4`.
                    let idx = unsafe { x_be4(packed, pos / COMPRESSION_RATIO) };
                    if !on_kmer(pos, ((idx >> shift) as u64) & mask) {
                        return;
                    }
                    pos += scan_step;
                }
            }
            while pos <= end {
                let kmer = packed_kmer_at(packed, pos, lut_word_length);
                if !on_kmer(pos, kmer) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
            }
            return;
        }
    } else if packed_len >= 3 {
        let max_byte_idx = packed_len.saturating_sub(3);
        let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
        let fast_end = end.min(max_fast_pos);
        let mut pos = start;
        // The words up to `fast_end`, counted once: each word's three bytes at its byte, at most
        // `max_byte_idx` there, with no test of the length (as `x_whole_rounds`).
        if pos <= fast_end {
            for _ in 0..(fast_end - pos) / scan_step + 1 {
                let byte_idx = pos / COMPRESSION_RATIO;
                let shift = 2 * (12 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                // SAFETY: `pos <= fast_end`, so `byte_idx + 2 <= max_byte_idx + 2 = packed_len - 1`.
                let idx = unsafe {
                    (x_byte(packed, byte_idx) << 16)
                        | (x_byte(packed, byte_idx + 1) << 8)
                        | x_byte(packed, byte_idx + 2)
                };
                if !on_kmer(pos, ((idx >> shift) as u64) & mask) {
                    return;
                }
                pos += scan_step;
            }
        }
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, lut_word_length);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mut pos = start;
    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, lut_word_length);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1639-1701
#[inline(always)]
fn scan_subject_kmers_range_mb_9_1<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 9;
    debug_assert_eq!(scan_step, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;
    let mut init_index: u32 = 0;
    let mut init_valid = false;
    let mut rounds_done = false;

    while pos <= end {
        if !rounds_done && !init_valid && state == 0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, the four words of one byte.
            rounds_done = true;
            let round = XRound {
                steps: 4,
                bytes: 1,
                reads_end: 4,
            };
            for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + 4`.
                let w = unsafe { x_be4(packed, byte_idx) } >> 8;
                for (k, shift) in [6u32, 4, 2, 0].into_iter().enumerate() {
                    if !on_kmer(pos + k, ((w >> shift) as u64) & mask) {
                        return;
                    }
                }
                pos += 4;
                byte_idx += 1;
            }
            continue;
        }
        if !init_valid {
            let Some(v) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            init_index = v;
            init_valid = true;
        }

        match state {
            0 => {
                let index = ((init_index >> 6) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let index = ((init_index >> 4) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = ((init_index >> 2) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                let index = (init_index as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 0;
                byte_idx = byte_idx.saturating_add(1);
                init_valid = false;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1715-1755
#[inline(always)]
fn scan_subject_kmers_range_mb_9_2<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 9;
    debug_assert_eq!(scan_step, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_shift2 = pos % COMPRESSION_RATIO == 2;
    let mut init_index: u32 = 0;
    let mut init_valid = false;
    let mut rounds_done = false;

    while pos <= end {
        if !rounds_done && !init_valid && !use_shift2 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, the two words of one byte.
            rounds_done = true;
            let round = XRound {
                steps: 2,
                bytes: 1,
                reads_end: 4,
            };
            for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + 4`.
                let w = unsafe { x_be4(packed, byte_idx) } >> 8;
                if !on_kmer(pos, ((w >> 6) as u64) & mask) {
                    return;
                }
                if !on_kmer(pos + scan_step, ((w >> 2) as u64) & mask) {
                    return;
                }
                pos += 2 * scan_step;
                byte_idx += 1;
            }
            continue;
        }
        if !init_valid {
            let Some(v) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            init_index = v;
            init_valid = true;
        }

        if !use_shift2 {
            let index = ((init_index >> 6) as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            use_shift2 = true;
        } else {
            let index = ((init_index >> 2) as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(1);
            init_valid = false;
            use_shift2 = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1768-1832
#[inline(always)]
fn scan_subject_kmers_range_mb_10_1<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;
    let mut init_index: u32 = 0;
    let mut init_valid = false;
    let mut rounds_done = false;

    while pos <= end {
        if !rounds_done && !init_valid && state == 0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, the four words of one byte (the
            // fourth with the next byte shifted in).
            rounds_done = true;
            let round = XRound {
                steps: 4,
                bytes: 1,
                reads_end: 4,
            };
            for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + 4`.
                let w = unsafe { x_be4(packed, byte_idx) };
                for (k, shift) in [12u32, 10, 8, 6].into_iter().enumerate() {
                    if !on_kmer(pos + k, ((w >> shift) as u64) & mask) {
                        return;
                    }
                }
                pos += 4;
                byte_idx += 1;
            }
            continue;
        }
        if !init_valid {
            let Some(v) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            init_index = v;
            init_valid = true;
        }

        match state {
            0 => {
                let index = ((init_index >> 4) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let index = ((init_index >> 2) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = (init_index as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                let Some(s3) = x_be_bytes::<1>(packed, byte_idx + 3) else {
                    break;
                };
                init_index = (init_index << 8) | s3;
                let index = ((init_index >> 6) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 0;
                byte_idx = byte_idx.saturating_add(1);
                init_valid = false;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1845-1884
#[inline(always)]
fn scan_subject_kmers_range_mb_10_2<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 3 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_shift0 = pos % COMPRESSION_RATIO == 2;
    let mut init_index: u32 = 0;
    let mut init_valid = false;
    let mut rounds_done = false;

    while pos <= end {
        if !rounds_done && !init_valid && !use_shift0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, the two words of one byte.
            rounds_done = true;
            let round = XRound {
                steps: 2,
                bytes: 1,
                reads_end: 4,
            };
            for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + 4`.
                let w = unsafe { x_be4(packed, byte_idx) } >> 8;
                if !on_kmer(pos, ((w >> 4) as u64) & mask) {
                    return;
                }
                if !on_kmer(pos + scan_step, (w as u64) & mask) {
                    return;
                }
                pos += 2 * scan_step;
                byte_idx += 1;
            }
            continue;
        }
        if !init_valid {
            let Some(v) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            init_index = v;
            init_valid = true;
        }

        if !use_shift0 {
            let index = ((init_index >> 4) as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            use_shift0 = true;
        } else {
            let index = (init_index as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(1);
            init_valid = false;
            use_shift0 = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1898-1964
#[inline(always)]
fn scan_subject_kmers_range_mb_10_3<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 10;
    debug_assert_eq!(scan_step, 3);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 5 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let packed_len_i32 = packed_len as i32;
    let mut pos = start;
    let mut byte_idx = (pos / COMPRESSION_RATIO) as i32;
    let mut state = 0u8;
    let mut init_index: u32 = 0;
    let mut use_fast = true;

    match pos % COMPRESSION_RATIO {
        1 => {
            if byte_idx + 1 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 8)
                    | (packed[(byte_idx + 1) as usize] as u32);
                byte_idx -= 2;
                state = 3;
            }
        }
        2 => {
            if byte_idx + 2 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 16)
                    | ((packed[(byte_idx + 1) as usize] as u32) << 8)
                    | (packed[(byte_idx + 2) as usize] as u32);
                byte_idx -= 1;
                state = 2;
            }
        }
        3 => {
            if byte_idx + 2 >= packed_len_i32 {
                use_fast = false;
            } else {
                init_index = ((packed[byte_idx as usize] as u32) << 16)
                    | ((packed[(byte_idx + 1) as usize] as u32) << 8)
                    | (packed[(byte_idx + 2) as usize] as u32);
                state = 1;
            }
        }
        _ => {
            state = 0;
        }
    }

    if !use_fast {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mut rounds_done = false;
    while pos <= end {
        if !rounds_done && state == 0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, four words over three bytes, the
            // accumulator shifting in the fourth and the fifth byte.
            rounds_done = true;
            if let Ok(mut b) = usize::try_from(byte_idx) {
                let round = XRound {
                    steps: 4,
                    bytes: 3,
                    reads_end: 5,
                };
                for _ in 0..x_whole_rounds(round, scan_step, pos, end, b, packed_len) {
                    // SAFETY: `x_whole_rounds`: this round's reads end before `b + 5`.
                    let (w, s4) = unsafe { (x_be4(packed, b), x_byte(packed, b + 4)) };
                    let w5 = (w << 8) | s4;
                    for (k, word) in [w >> 12, w >> 6, w, w5 >> 2].into_iter().enumerate() {
                        if !on_kmer(pos + k * scan_step, (word as u64) & mask) {
                            return;
                        }
                    }
                    init_index = w5;
                    pos += 4 * scan_step;
                    b += 3;
                }
                byte_idx = b as i32;
                continue;
            }
        }
        match state {
            0 => {
                let Some(v) = x_be_bytes_from::<3>(packed, byte_idx) else {
                    break;
                };
                init_index = v;
                let index = ((init_index >> 4) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 1;
            }
            1 => {
                let Some(b) = x_be_bytes_from::<1>(packed, byte_idx + 3) else {
                    break;
                };
                init_index = (init_index << 8) | b;
                let index = ((init_index >> 6) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 2;
            }
            2 => {
                let index = (init_index as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                state = 3;
            }
            _ => {
                let Some(b) = x_be_bytes_from::<1>(packed, byte_idx + 4) else {
                    break;
                };
                init_index = (init_index << 8) | b;
                let index = ((init_index >> 2) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx += 3;
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1978-2041
#[inline(always)]
fn scan_subject_kmers_range_mb_11_1mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 1);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut state = pos % COMPRESSION_RATIO;
    let mut rounds_done = false;

    while pos <= end {
        if !rounds_done && state == 0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, one word at each of the four
            // places of a byte.
            rounds_done = true;
            let sb = scan_step_byte;
            let round = XRound {
                steps: 4,
                bytes: 4 * sb + 1,
                reads_end: 3 * sb + 4,
            };
            for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + 3 sb + 4`.
                let words = unsafe {
                    [
                        x_be4(packed, byte_idx) >> 10,
                        x_be4(packed, byte_idx + sb) >> 8,
                        x_be4(packed, byte_idx + 2 * sb) >> 6,
                        x_be4(packed, byte_idx + 3 * sb) >> 4,
                    ]
                };
                for (k, word) in words.into_iter().enumerate() {
                    if !on_kmer(pos + k * scan_step, (word as u64) & mask) {
                        return;
                    }
                }
                pos += 4 * scan_step;
                byte_idx += 4 * sb + 1;
            }
            continue;
        }
        match state {
            0 => {
                let Some(idx) = x_be_bytes::<3>(packed, byte_idx) else {
                    break;
                };
                let index = ((idx >> 2) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 1;
            }
            1 => {
                let Some(idx) = x_be_bytes::<3>(packed, byte_idx) else {
                    break;
                };
                let index = (idx as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 2;
            }
            2 => {
                let Some(idx) = x_be_bytes::<4>(packed, byte_idx) else {
                    break;
                };
                let index = ((idx >> 6) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte);
                state = 3;
            }
            _ => {
                let Some(idx) = x_be_bytes::<4>(packed, byte_idx) else {
                    break;
                };
                let index = ((idx >> 4) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx = byte_idx.saturating_add(scan_step_byte + 1);
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2055-2107
#[inline(always)]
fn scan_subject_kmers_range_mb_11_2mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 2);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 4 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let (top_shift, bottom_shift) = if start.is_multiple_of(2) {
        (2u32, 6u32)
    } else {
        (0u32, 4u32)
    };
    let mut pos = start;
    let mut byte_idx = pos / COMPRESSION_RATIO;
    let mut use_bottom = pos % COMPRESSION_RATIO >= 2;
    let mut rounds_done = false;

    while pos <= end {
        if !use_bottom {
            if !rounds_done {
                // The rounds (`x_whole_rounds`): NCBI's loop body, a top read, then a bottom
                // read.
                rounds_done = true;
                let sb = scan_step_byte;
                let round = XRound {
                    steps: 2,
                    bytes: 2 * sb + 1,
                    reads_end: sb + 4,
                };
                // (the top read's three bytes are the four-byte read `>> 8`: one shift)
                let top_shift4 = top_shift + 8;
                for _ in 0..x_whole_rounds(round, scan_step, pos, end, byte_idx, packed_len) {
                    // SAFETY: `x_whole_rounds`: this round's reads end before `byte_idx + sb + 4`.
                    let (top, bottom) =
                        unsafe { (x_be4(packed, byte_idx), x_be4(packed, byte_idx + sb)) };
                    if !on_kmer(pos, ((top >> top_shift4) as u64) & mask) {
                        return;
                    }
                    if !on_kmer(pos + scan_step, ((bottom >> bottom_shift) as u64) & mask) {
                        return;
                    }
                    pos += 2 * scan_step;
                    byte_idx += 2 * sb + 1;
                }
                continue;
            }
            let Some(idx) = x_be_bytes::<3>(packed, byte_idx) else {
                break;
            };
            let index = ((idx >> top_shift) as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(scan_step_byte);
            use_bottom = true;
        } else {
            let Some(idx) = x_be_bytes::<4>(packed, byte_idx) else {
                break;
            };
            let index = ((idx >> bottom_shift) as u64) & mask;
            if !on_kmer(pos, index) {
                return;
            }
            pos = pos.saturating_add(scan_step);
            byte_idx = byte_idx.saturating_add(scan_step_byte + 1);
            use_bottom = false;
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2120-2188
#[inline(always)]
fn scan_subject_kmers_range_mb_11_3mod4<F>(
    packed: &[u8],
    subject_len: usize,
    scan_step: usize,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    const LUT_WORD_LENGTH: usize = 11;
    debug_assert_eq!(scan_step % COMPRESSION_RATIO, 3);
    if subject_len < LUT_WORD_LENGTH || start > end || start + LUT_WORD_LENGTH > subject_len {
        return;
    }

    let packed_len = packed.len();
    if packed_len < 5 {
        let mut pos = start;
        while pos <= end {
            let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
            if !on_kmer(pos, kmer) {
                return;
            }
            pos = pos.saturating_add(scan_step);
        }
        return;
    }

    let mask = (1u64 << (2 * LUT_WORD_LENGTH)) - 1;
    let scan_step_byte = scan_step / COMPRESSION_RATIO;
    let mut pos = start;
    let mut byte_idx = (pos / COMPRESSION_RATIO) as i32;
    let mut state: u8 = match pos % COMPRESSION_RATIO {
        1 => {
            byte_idx -= 2;
            3
        }
        2 => {
            byte_idx -= 1;
            2
        }
        3 => 1,
        _ => 0,
    };

    let mut rounds_done = false;
    while pos <= end {
        if !rounds_done && state == 0 {
            // The rounds (`x_whole_rounds`): NCBI's loop body, one word at each of the four
            // places of a byte.
            rounds_done = true;
            if let Ok(mut b) = usize::try_from(byte_idx) {
                let sb = scan_step_byte;
                let round = XRound {
                    steps: 4,
                    bytes: 4 * sb + 3,
                    reads_end: 3 * sb + 6,
                };
                for _ in 0..x_whole_rounds(round, scan_step, pos, end, b, packed_len) {
                    // SAFETY: `x_whole_rounds`: this round's reads end before `b + 3 sb + 6`.
                    let words = unsafe {
                        [
                            x_be4(packed, b) >> 10,
                            x_be4(packed, b + sb) >> 4,
                            x_be4(packed, b + 2 * sb + 1) >> 6,
                            x_be4(packed, b + 3 * sb + 2) >> 8,
                        ]
                    };
                    for (k, word) in words.into_iter().enumerate() {
                        if !on_kmer(pos + k * scan_step, (word as u64) & mask) {
                            return;
                        }
                    }
                    pos += 4 * scan_step;
                    b += 4 * sb + 3;
                }
                byte_idx = b as i32;
                continue;
            }
        }
        match state {
            0 => {
                let Some(v) = x_be_bytes_from::<3>(packed, byte_idx) else {
                    break;
                };
                let idx = v;
                let index = ((idx >> 2) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 1;
            }
            1 => {
                let Some(idx) = x_be_bytes_from::<4>(packed, byte_idx) else {
                    break;
                };
                let index = ((idx >> 4) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 2;
            }
            2 => {
                let Some(idx) = x_be_bytes_from::<4>(packed, byte_idx + 1) else {
                    break;
                };
                let index = ((idx >> 6) as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx += scan_step_byte as i32;
                state = 3;
            }
            _ => {
                let Some(idx) = x_be_bytes_from::<3>(packed, byte_idx + 2) else {
                    break;
                };
                let index = (idx as u64) & mask;
                if !on_kmer(pos, index) {
                    return;
                }
                pos = pos.saturating_add(scan_step);
                byte_idx += (scan_step_byte + 3) as i32;
                state = 0;
            }
        }
    }

    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, LUT_WORD_LENGTH);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1635-1676
#[derive(Clone, Copy)]
enum MbScanSubjectKind {
    Any,
    Scan9_1,
    Scan9_2,
    Scan10_1,
    Scan10_2,
    Scan10_3,
    Scan11_1Mod4,
    Scan11_2Mod4,
    Scan11_3Mod4,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2602-2677
#[inline(always)]
fn choose_mb_scan_subject_kind(lut_word_length: usize, scan_step: usize) -> MbScanSubjectKind {
    match lut_word_length {
        9 => {
            if scan_step == 1 {
                MbScanSubjectKind::Scan9_1
            } else if scan_step == 2 {
                MbScanSubjectKind::Scan9_2
            } else {
                MbScanSubjectKind::Any
            }
        }
        10 => {
            if scan_step == 1 {
                MbScanSubjectKind::Scan10_1
            } else if scan_step == 2 {
                MbScanSubjectKind::Scan10_2
            } else if scan_step == 3 {
                MbScanSubjectKind::Scan10_3
            } else {
                MbScanSubjectKind::Any
            }
        }
        11 => match scan_step % COMPRESSION_RATIO {
            0 => MbScanSubjectKind::Any,
            1 => MbScanSubjectKind::Scan11_1Mod4,
            2 => MbScanSubjectKind::Scan11_2Mod4,
            _ => MbScanSubjectKind::Scan11_3Mod4,
        },
        12 => MbScanSubjectKind::Any,
        _ => MbScanSubjectKind::Any,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2602-2677
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1651-1667
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:2994-3006
#[inline(always)]
fn select_mb_scan_kind(
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
) -> Option<MbScanSubjectKind> {
    let base = if (9..=12).contains(&lut_word_length) {
        Some(choose_mb_scan_subject_kind(lut_word_length, scan_step))
    } else {
        None
    };
    if subject_masked {
        base.map(|_| MbScanSubjectKind::Any)
    } else {
        base
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:56-66
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:88-91
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:119-123
#[inline(always)]
fn mask_lookup_index(index: u64, lut_word_length: usize) -> u64 {
    let mask = (1u64 << (2 * lut_word_length)) - 1;
    index & mask
}

fn scan_subject_kmers_range<F>(
    packed: &[u8],
    subject_len: usize,
    lut_word_length: usize,
    scan_step: usize,
    subject_masked: bool,
    start: usize,
    end: usize,
    on_kmer: &mut F,
) where
    F: FnMut(usize, u64) -> bool,
{
    if subject_len < lut_word_length || lut_word_length == 0 {
        return;
    }
    if start > end || start + lut_word_length > subject_len {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:165-176
    debug_assert!(scan_step > 0);

    if lut_word_length <= 8 {
        let mask = (1u64 << (2 * lut_word_length)) - 1;
        let packed_len = packed.len();
        if packed_len == 0 {
            return;
        }
        if lut_word_length > 5 {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:172-209
            if scan_step.is_multiple_of(COMPRESSION_RATIO)
                && !subject_masked
                && start.is_multiple_of(COMPRESSION_RATIO)
                && packed_len >= 2
            {
                let shift = 2 * (8 - lut_word_length);
                let step_bytes = scan_step / COMPRESSION_RATIO;
                let mut byte_idx = start / COMPRESSION_RATIO;
                let max_byte_idx = packed_len.saturating_sub(2);
                let end_byte_idx = (end / COMPRESSION_RATIO).min(max_byte_idx);
                while byte_idx <= end_byte_idx {
                    let idx = ((packed[byte_idx] as u16) << 8) | packed[byte_idx + 1] as u16;
                    let kmer = ((idx >> shift) as u64) & mask;
                    if !on_kmer(byte_idx * COMPRESSION_RATIO, kmer) {
                        return;
                    }
                    byte_idx += step_bytes;
                }
                let mut pos = byte_idx * COMPRESSION_RATIO;
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    if !on_kmer(pos, kmer) {
                        return;
                    }
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }

            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:210-240
            if packed_len >= 3 {
                let max_byte_idx = packed_len.saturating_sub(3);
                let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
                let fast_end = end.min(max_fast_pos);
                let mut pos = start;
                while pos <= fast_end {
                    let byte_idx = pos / COMPRESSION_RATIO;
                    let shift = 2 * (12 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                    let idx = ((packed[byte_idx] as u32) << 16)
                        | ((packed[byte_idx + 1] as u32) << 8)
                        | (packed[byte_idx + 2] as u32);
                    let kmer = ((idx >> shift) as u64) & mask;
                    if !on_kmer(pos, kmer) {
                        return;
                    }
                    pos = pos.saturating_add(scan_step);
                }
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    if !on_kmer(pos, kmer) {
                        return;
                    }
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }
        } else {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:245-260
            if packed_len >= 2 {
                let max_byte_idx = packed_len.saturating_sub(2);
                let max_fast_pos = max_byte_idx * COMPRESSION_RATIO + (COMPRESSION_RATIO - 1);
                let fast_end = end.min(max_fast_pos);
                let mut pos = start;
                while pos <= fast_end {
                    let byte_idx = pos / COMPRESSION_RATIO;
                    let shift = 2 * (8 - ((pos % COMPRESSION_RATIO) + lut_word_length));
                    let idx = ((packed[byte_idx] as u16) << 8) | packed[byte_idx + 1] as u16;
                    let kmer = ((idx >> shift) as u64) & mask;
                    if !on_kmer(pos, kmer) {
                        return;
                    }
                    pos = pos.saturating_add(scan_step);
                }
                while pos <= end {
                    let kmer = packed_kmer_at(packed, pos, lut_word_length);
                    if !on_kmer(pos, kmer) {
                        return;
                    }
                    pos = pos.saturating_add(scan_step);
                }
                return;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nascan.c:1598-1605
    // Fallback: advance by scan_step and recompute the k-mer at each position.
    let mut pos = start;
    while pos <= end {
        let kmer = packed_kmer_at(packed, pos, lut_word_length);
        if !on_kmer(pos, kmer) {
            return;
        }
        pos = pos.saturating_add(scan_step);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/masksubj.inl:43-58
#[inline(always)]
fn determine_scanning_offsets(
    seq_ranges: &[(i32, i32)],
    word_length: i32,
    lut_word_length: i32,
    range: &mut [i32; 3],
) -> bool {
    while range[1] > range[2] {
        range[0] += 1;
        if range[0] >= seq_ranges.len() as i32 {
            return false;
        }
        let (left, right) = seq_ranges[range[0] as usize];
        range[1] = left + word_length - lut_word_length;
        range[2] = right - lut_word_length;
    }
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1674
/// The scan range a subject's scan starts from (`seq_ranges` not empty): the whole subject when
/// it is unmasked, else its first range; `determine_scanning_offsets` moves on to the next.
#[inline(always)]
fn first_scan_range(
    seq_ranges: &[(i32, i32)],
    subject_len: usize,
    word_length: usize,
    lut_word_length: usize,
    subject_masked: bool,
) -> [i32; 3] {
    if subject_masked {
        let (left, right) = seq_ranges[0];
        [
            0,
            left + word_length as i32 - lut_word_length as i32,
            right - lut_word_length as i32,
        ]
    } else {
        [0, 0, subject_len as i32 - lut_word_length as i32]
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1673-1676
#[inline(always)]
fn scan_subject_kmers_with_offsets<G>(
    seq_ranges: &[(i32, i32)],
    word_length: usize,
    lut_word_length: usize,
    mut scan_range: [i32; 3],
    mut scan_fn: G,
) where
    G: FnMut(usize, usize),
{
    while determine_scanning_offsets(
        seq_ranges,
        word_length as i32,
        lut_word_length as i32,
        &mut scan_range,
    ) {
        let start = scan_range[1];
        let end = scan_range[2];
        if start >= 0 && end >= 0 {
            scan_fn(start as usize, end as usize);
        }
        scan_range[1] = scan_range[2] + 1;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1674
fn scan_subject_kmers_with_ranges<F>(
    packed: &[u8],
    subject_len: usize,
    word_length: usize,
    lut_word_length: usize,
    scan_step: usize,
    seq_ranges: &[(i32, i32)],
    subject_masked: bool,
    mb_scan_kind: Option<MbScanSubjectKind>,
    mut on_kmer: F,
) where
    F: FnMut(usize, usize, u64),
{
    if seq_ranges.is_empty() {
        return;
    }

    let scan_range = first_scan_range(
        seq_ranges,
        subject_len,
        word_length,
        lut_word_length,
        subject_masked,
    );

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1635-1676
    match mb_scan_kind {
        Some(MbScanSubjectKind::Any) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_any(
                    packed,
                    subject_len,
                    lut_word_length,
                    scan_step,
                    subject_masked,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan9_1) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_9_1(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan9_2) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_9_2(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_1) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_1(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_2) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_2(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan10_3) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_10_3(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_1Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_1mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_2Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_2mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        Some(MbScanSubjectKind::Scan11_3Mod4) => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range_mb_11_3mod4(
                    packed,
                    subject_len,
                    scan_step,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
        None => scan_subject_kmers_with_offsets(
            seq_ranges,
            word_length,
            lut_word_length,
            scan_range,
            |start, end| {
                let s_range = end.saturating_add(lut_word_length);
                scan_subject_kmers_range(
                    packed,
                    subject_len,
                    lut_word_length,
                    scan_step,
                    subject_masked,
                    start,
                    end,
                    &mut |kmer_start, kmer| {
                        on_kmer(kmer_start, s_range, kmer);
                        true
                    },
                );
            },
        ),
    }
}

/// Scans a subject for the discontiguous words of a discontiguous megablast table, with
/// the template length as both the word length and the lookup word length; calls
/// `on_word(s_off, s_range, index, index2)` for every word in NCBI's order.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1624-1634
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1684
fn scan_subject_disc_words_with_ranges<F>(
    packed: &[u8],
    subject_len: usize,
    disc: &DiscTemplates,
    seq_ranges: &[(i32, i32)],
    subject_masked: bool,
    mut on_word: F,
) where
    F: FnMut(usize, usize, u64, Option<u64>),
{
    if seq_ranges.is_empty() {
        return;
    }
    let template_length = disc.template_length;
    let scan_range = first_scan_range(
        seq_ranges,
        subject_len,
        template_length,
        template_length,
        subject_masked,
    );
    let kind = choose_disc_scan_subject(disc.template_type, disc.two_templates);
    scan_subject_kmers_with_offsets(
        seq_ranges,
        template_length,
        template_length,
        scan_range,
        |start, end| {
            let s_range = end.saturating_add(template_length);
            disc_word_scan_subject(
                kind,
                packed,
                disc.template_type,
                disc.second_template_type,
                template_length as i32,
                [start as i32, end as i32],
                &mut |s_off, index, index2| {
                    on_word(s_off as usize, s_range, index as u64, index2.map(u64::from))
                },
            );
        },
    );
}

/// Structure to hold ungapped hit data for batch processing
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c - init_hsp_array is sorted by score descending
/// before gapped extension
#[derive(Clone, Copy)]
struct UngappedHit {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3908-3913
    context_idx: u32,
    query_idx: u32,
    query_frame: i32,
    query_context_offset: i32,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:135-150
    // Initial word hit offsets (seed) stored in init_hsp->offsets.qs_offsets.
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-147
    // The offsets are 32-bit, as NCBI's: those in the query are within one context (whose
    // start, `query_context_offset`, is an `i32`), those in the subject within one chunk of at
    // most `MAX_DBSEQ_LEN` letters. A subject's hits and the scratch of their sort by score
    // (`X_SORT`) then hold 44 bytes per hit each, not 72 (pyani's fragments of one genome x
    // another, about a million hits per chunk: 69 MiB each).
    seed_q_off: u32,
    seed_s_off: u32,
    // Ungapped extension results (0-based coordinates).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3901-3904
    // `qs` is NCBI's context-relative `Int4` `q_start` (`UngappedHit::saved`), which is negative
    // for a hit that starts left of its context, and `qe` its `q_end`; `qe - qs` is NCBI's
    // `length`. The subject offsets are never negative (the left walks stop at the subject's
    // first base).
    qs: i32,    // query start
    qe: i32,    // query end (exclusive)
    ss: u32,    // subject start in search_seq coordinates
    se: u32,    // subject end in search_seq coordinates
    score: i32, // ungapped score
}

impl UngappedHit {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:752-756 (the same at
    // 923-927)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2384-2392
    /// The initial hit that NCBI saves for the word hit at `q_off` (an offset in the concatenated
    /// query) and `s_off` of context `context_idx` (`ctx`), with the ungapped extension
    /// `ungapped` (offsets in the concatenated query), and whose offsets
    /// `s_AdjustInitialHSPOffsets` makes context-relative. NCBI's extensions run on the
    /// concatenated query (na_ungapped.c:164-206 and 284-310), so a saved hit can start left of
    /// its context: where the exact walk passes two sentinels (half of `INT4_MIN` each,
    /// blast_stat.c:1124-1127) and its `Int4` sum wraps above 0 (a raw X-drop of 2^30 or more,
    /// such as the `INT_MIN` of `-xdrop_ungap 1e10`, does not stop the walk at the first one),
    /// or where the approximate walk, which reads the sentinel as bases, keeps a hit below the
    /// cutoff (`off_found`). Its `q_start` is then a negative `Int4` (the `ASSERT` above is not
    /// compiled into NCBI's release), and `qs` here is that value.
    #[inline(always)]
    fn saved(
        context_idx: u32,
        ctx: &QueryContext,
        q_off: usize,
        s_off: usize,
        ungapped: &UngappedData,
    ) -> Self {
        let qs = ungapped.q_start as i32 - ctx.query_offset;
        UngappedHit {
            context_idx,
            query_idx: ctx.query_idx,
            query_frame: ctx.frame,
            query_context_offset: ctx.query_offset,
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_def.h:142-145
            // The word hit lies in its context, so its offset is not negative.
            seed_q_off: (q_off - ctx.query_offset as usize) as u32,
            seed_s_off: s_off as u32,
            qs,
            qe: qs + ungapped.length as i32,
            ss: ungapped.s_start as u32,
            se: (ungapped.s_start + ungapped.length) as u32,
            score: ungapped.score,
        }
    }

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-147
    /// NCBI's `ungapped_data->length`.
    #[inline(always)]
    fn length(&self) -> i32 {
        self.qe - self.qs
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4017
    /// The seed of the greedy alignment, in NCBI's `Int4` (its `Uint4` fields hold these
    /// values, which `BLAST_GreedyGappedAlignment` takes as `Int4`): the middle of the hit.
    /// The query offset is negative where more of the hit lies left of its context than in it.
    #[inline(always)]
    fn greedy_seed(&self) -> (i32, i32) {
        let half = self.length() / 2;
        (self.qs + half, self.ss as i32 + half)
    }
}

// No NCBI counterpart: only where LOSAT keeps its buffers (no value changes).
/// EXPERIMENT (LOSAT_X_BNMAPPED, a switch of the strict set, `utils/xswitch.rs`: on unless it is
/// `0`): the large per-subject buffers of BLASTN (the
/// blocks and sequences read from a database and their BLASTNA letters, `database.rs`; the
/// ungapped hits of a subject as the seed stages push them, `x_push_ungapped`, and the scratch of
/// their sort by score, outside the per-query groups of `LOSAT_X_QGROUPPAR` and the merge of the
/// folded strips, `x_pair_seed_part.rs`) are kept in buffers that glibc maps
/// from 1 MiB on, so that they are returned when they are dropped instead of staying in the
/// arena of the worker thread that freed them; the database's buffers and the sort's scratch go
/// to shared pools (`utils::x_mapped::X_BYTES`, `X_SORT`) when they are done, so that the next
/// thread that needs one takes it with its pages resident. A worker of the preliminary stage
/// reads its blocks below 1 MiB into one buffer that it keeps from one block to the next, only its
/// bytes past those of the last block zeroed (`database.rs` `x_read_block`; a block of 1 MiB or
/// more, one long sequence, takes a buffer of the pool). The hits are not pooled: a thread
/// keeps its scratch, and with it its hit buffer, from one subject to the next (as NCBI keeps its
/// `init_hitlist`), so idle hit buffers in a pool only added to what was in use (Sakai x EDL933
/// `-subject` blastn at 4 threads, max RSS 298-327 MiB with a pool of 64 MiB, 245-264 without,
/// 238-261 with the switch off).
fn x_bn_mapped() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNMAPPED"))
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_MBHITPF, in the strict set: on unless it is `0`) the
/// first extension of a hit reads the query at a place the hit's cell chose, a cache miss for a
/// query of megabases; the loop over the offset pairs prefetches the line of the hit `X_HIT_AHEAD`
/// places on. It changes no value.
fn x_mb_hit_prefetch() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| {
        let on = crate::utils::xswitch::on("LOSAT_X_MBHITPF");
        crate::utils::xreach::x_reach!("LOSAT_X_MBHITPF", on);
        on
    })
}
const X_HIT_AHEAD: usize = 16;

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNPVONCE, in the strict set: on unless it is `0`) the
/// scan of a table laid out as packed CSR tests each word once, in the head's presence bits, when
/// those fit the second-level cache (4^11 cells or fewer); a larger head keeps the small presence
/// vector in front. The words it keeps are the words with hits either way (the presence vector
/// only skips words without hits). It changes no value.
fn x_bn_pv_once(head: &x_dc_csr::XPacked<'_>) -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let on = *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNPVONCE"));
    let once = on && head.x_bits_bytes() <= X_PV_ONCE_MAX_BITS;
    crate::utils::xreach::x_reach!("LOSAT_X_BNPVONCE", once);
    once
}
const X_PV_ONCE_MAX_BITS: usize = 1 << 20;

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNEXTFIRST, in the strict set: on unless it is `0`)
/// `pair_extend::extend_offset_pairs` extends a megablast table's lookup hit to its word before it
/// looks up the hit's context, so a hit whose word is not found costs only its extension, as in
/// NCBI's `s_BlastNaExtend`, where the table's lookup width is below the word size. The words and
/// every value are the same.
pub(super) fn x_bn_ext_first() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNEXTFIRST"))
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNRANKONCE, in the strict set: on unless it is `0`)
/// the scans of a contiguous table laid out as CSR take each word's run from the element of the
/// head that the ring hands over with the word (`x_dc_csr::XRing`): a word of a packed head is
/// ranked once, where it enters the ring, not again where it leaves. The words, their runs and
/// every value are the same. Recorded (`LOSAT_X_REACH`) where a scan chooses the table's form.
pub(super) fn x_bn_rank_once(csr: bool) -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    let on = *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNRANKONCE"));
    if on {
        crate::utils::xreach::x_reach!("LOSAT_X_BNRANKONCE", csr);
    }
    on && csr
}
/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNTBRANGE) the traceback of a database search decodes
/// a sequence of at least `x_tb_range_min` bases in the ranges of its HSPs only, as NCBI's partial
/// fetch does (`database::traceback_listed`).
pub(super) fn x_bn_tb_range() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNTBRANGE"))
}

/// No NCBI counterpart: the shortest sequence that LOSAT_X_BNTBRANGE decodes in ranges, 1 Mi
/// bases (below it the whole BLASTNA is small, and an allocation below glibc's mmap threshold
/// may come from an arena that writes it); `LOSAT_X_BNTBRANGEMIN=<bases>` sets it (0 = every
/// sequence, for the tests: the fixtures' sequences are short).
pub(super) fn x_tb_range_min() -> usize {
    static MIN: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *MIN.get_or_init(|| {
        std::env::var("LOSAT_X_BNTBRANGEMIN")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .unwrap_or(crate::utils::x_mapped::X_MAPPED_FROM)
    })
}

/// No NCBI counterpart: records for `LOSAT_X_REACH` whether a buffer of `bytes` takes
/// LOSAT_X_BNMAPPED's own path (mapped, and pooled where it has a pool) or stays as `Vec` keeps
/// it (below 1 MiB, or not glibc).
pub(super) fn x_bn_mapped_reach(bytes: usize) {
    crate::utils::xreach::x_reach!(
        "LOSAT_X_BNMAPPED",
        crate::utils::x_mapped::x_is_mapped(bytes)
    );
}

/// No NCBI counterpart: EXPERIMENT (LOSAT_X_BNMAPPED) the pool of the scratch of the hits' sort
/// by score (`x_mapped::x_sort_by_mapped`), whose idle buffers may hold `x_sort_idle_bytes`.
/// Without it each sort maps its scratch again and touches every page of it (Sakai x EDL933
/// `-subject` blastn, where the hits of one chunk are sorted at a time: x1.009-1.027 of the time
/// with the switch off). The idle buffers are freed when a query batch's preliminary stage ends
/// (`x_clear_idle` next to `x_lut_drop::after_prelim`); two searches of one process at once share
/// the pool, so one's clear may free the other's idle buffer (time only).
static X_SORT: crate::utils::x_mapped::XPool<UngappedHit> =
    crate::utils::x_mapped::XPool::new(x_sort_idle_bytes);

/// No NCBI counterpart: the built-in bound of `x_sort_idle_bytes`, in MiB: one scratch of the
/// largest hit lists measured (Sakai x EDL933 `-subject` blastn sorts 10 lists of 1.0-11.4 MiB;
/// with 8 each of the 5 above 9 MiB is mapped again, x1.009 of the time).
const X_SORT_IDLE_MIB: usize = 16;

/// No NCBI counterpart: the idle bound of `X_SORT` (`LOSAT_X_SORTIDLE` MiB, a tuning variable,
/// 0 to 65,536; built-in `X_SORT_IDLE_MIB`).
fn x_sort_idle_bytes() -> usize {
    static IDLE: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *IDLE.get_or_init(|| crate::utils::x_mapped::x_idle_bytes("LOSAT_X_SORTIDLE", X_SORT_IDLE_MIB))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-355
/// `hits.push(hit)`, growing as `Vec` does; with LOSAT_X_BNMAPPED with the capacity of
/// `x_mapped::x_reserve_mapped`.
#[inline(always)]
fn x_push_ungapped(hits: &mut Vec<UngappedHit>, hit: UngappedHit) {
    if hits.len() == hits.capacity() {
        x_grow_ungapped(hits);
    }
    hits.push(hit);
}

#[cold]
#[inline(never)]
fn x_grow_ungapped(hits: &mut Vec<UngappedHit>) {
    if x_bn_mapped() {
        crate::utils::x_mapped::x_reserve_mapped(hits, 1);
    } else {
        hits.reserve(1);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:271-300
// NCBI's `ungapped_data->q_start` is an offset in the concatenated query (the
// scan and the ungapped extension run on the whole query block, and the sort
// comes before `s_AdjustInitialHSPOffsets`), so hits of different contexts
// never tie on it; `qs` here is NCBI's context-relative value (negative for a
// hit that starts left of its context), and the context offset is added back:
// the `Int4` of the concatenated query, at least 0.
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:206-207
fn score_compare_ungapped_hits(a: &UngappedHit, b: &UngappedHit) -> std::cmp::Ordering {
    let a_len = a.length();
    let b_len = b.length();
    let a_q_start = a.query_context_offset + a.qs;
    let b_q_start = b.query_context_offset + b.qs;

    b.score
        .cmp(&a.score)
        .then_with(|| a.ss.cmp(&b.ss))
        .then_with(|| b_len.cmp(&a_len))
        .then_with(|| a_q_start.cmp(&b_q_start))
        .then_with(|| b_len.cmp(&a_len))
}

/// Preliminary gapped HSP data for greedy traceback.
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4031
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4058-4077
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:96-101
#[derive(Clone, Debug)]
struct PrelimHit {
    context_idx: u32,
    query_idx: u32,
    query_frame: i32,
    query_context_offset: i32,
    /// NCBI's `Int4` `query.offset` and `query.end` of the HSP in its context. They lie outside
    /// the context (`prelim_qs < 0`, or `prelim_qe` past the context's length) for an HSP from a
    /// hit that passes an end of its context (`UngappedHit::saved`): an `-ungapped` HSP
    /// (`prelim_hsp_list::ungapped_hsp_list`) or a greedy alignment from a seed outside the
    /// context (`prelim_hsp_list::x_prelim_extend`).
    prelim_qs: i32,
    prelim_qe: i32,
    prelim_ss: usize,
    prelim_se: usize,
    prelim_score: i32,
    /// NCBI's `Int4` `query.gapped_start` (outside the context where the greedy alignment's
    /// estimate of its start lies outside it) and `subject.gapped_start`.
    seed_qs: i32,
    seed_ss: usize,
    /// The e-value of the preliminary stage (set before its reap); a merge of two HSPs keeps
    /// the first one's (`s_BlastMergeTwoHSPs` does not change `evalue`).
    prelim_evalue: f64,
}

// No NCBI counterpart: the sizes of the hit records that each subject's lists hold (NCBI's
// `BlastInitHSP` + `BlastUngappedData`, `BlastHSP`; `PrelimHit`'s on 64-bit targets, where its
// `usize` fields are 8 bytes).
const _: () = assert!(std::mem::size_of::<UngappedHit>() == 44);
#[cfg(target_pointer_width = "64")]
const _: () = assert!(std::mem::size_of::<PrelimHit>() == 64);

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
fn score_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    b.prelim_score
        .cmp(&a.prelim_score)
        .then_with(|| a.prelim_ss.cmp(&b.prelim_ss))
        .then_with(|| b.prelim_se.cmp(&a.prelim_se))
        .then_with(|| a.prelim_qs.cmp(&b.prelim_qs))
        .then_with(|| b.prelim_qe.cmp(&a.prelim_qe))
}

type PrelimHitCompare = fn(&PrelimHit, &PrelimHit) -> std::cmp::Ordering;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2286-2320
fn query_offset_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    a.context_idx
        .cmp(&b.context_idx)
        .then_with(|| a.prelim_qs.cmp(&b.prelim_qs))
        .then_with(|| a.prelim_ss.cmp(&b.prelim_ss))
        .then_with(|| b.prelim_score.cmp(&a.prelim_score))
        .then_with(|| b.prelim_qe.cmp(&a.prelim_qe))
        .then_with(|| b.prelim_se.cmp(&a.prelim_se))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2351-2386
fn query_end_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    a.context_idx
        .cmp(&b.context_idx)
        .then_with(|| a.prelim_qe.cmp(&b.prelim_qe))
        .then_with(|| a.prelim_se.cmp(&b.prelim_se))
        .then_with(|| b.prelim_score.cmp(&a.prelim_score))
        .then_with(|| b.prelim_qs.cmp(&a.prelim_qs))
        .then_with(|| b.prelim_ss.cmp(&a.prelim_ss))
}

fn qsort_prelim_hits_by(hsps: &mut [PrelimHit], compare: PrelimHitCompare) {
    if hsps.len() <= 1 {
        return;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1379-1381
    // Keep the NCBI comparator fields and call timing while using the same
    // stable Rust sort on every target.
    hsps.sort_by(compare);
}

fn sort_prelim_hits_by_score(hsps: &mut [PrelimHit]) {
    if hsps.len() <= 1 {
        return;
    }

    let sorted = hsps
        .windows(2)
        .all(|pair| score_compare_prelim_hits(&pair[0], &pair[1]) != std::cmp::Ordering::Greater);
    if sorted {
        return;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:555
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1382
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:358-365
    // Preliminary HSPs drive traceback interval-tree insertion order. The
    // target-neutral stable sort preserves incoming order for comparator ties.
    qsort_prelim_hits_by(hsps, score_compare_prelim_hits);
}

fn purge_prelim_hits_with_common_endpoints(mut hsps: Vec<PrelimHit>) -> Vec<PrelimHit> {
    if hsps.len() <= 1 {
        return hsps;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:540-557
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2535
    // Preliminary BLASTN HSPs run with `purge` enabled, so common-start/common-end
    // duplicates are deleted, not trimmed. The full traceback stage still uses
    // the later trim (`purge` disabled) then delete (`purge` enabled) passes over final HSPs with edit scripts.
    qsort_prelim_hits_by(&mut hsps, query_offset_compare_prelim_hits);
    // In place: `dedup_by` drops an HSP that matches the last one kept, as the copy into a
    // second vector did.
    hsps.dedup_by(|hsp, prev| {
        prev.context_idx == hsp.context_idx
            && prev.prelim_qs == hsp.prelim_qs
            && prev.prelim_ss == hsp.prelim_ss
    });
    let mut kept = hsps;

    if kept.len() <= 1 {
        return kept;
    }

    qsort_prelim_hits_by(&mut kept, query_end_compare_prelim_hits);
    kept.dedup_by(|hsp, prev| {
        prev.context_idx == hsp.context_idx
            && prev.prelim_qe == hsp.prelim_qe
            && prev.prelim_se == hsp.prelim_se
    });

    kept
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits_priv.h:66-72
#[inline]
fn contained_in_hsp<T: PartialOrd>(a: T, b: T, c: T, d: T, e: T, f: T) -> bool {
    a <= c && b >= c && d <= f && e >= f
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1487-1531
fn merge_two_prelim_hits(hsp1: &mut PrelimHit, hsp2: &PrelimHit, allow_gap: bool) -> bool {
    // NCBI's `Int4` query offsets (`PrelimHit::prelim_qs`, which may lie outside the context).
    let (q1, q2) = (hsp1.prelim_qs as isize, hsp2.prelim_qs as isize);
    if !allow_gap && (hsp1.prelim_ss as isize - hsp2.prelim_ss as isize - q1 + q2) != 0 {
        return false;
    }

    // BLASTN subject frame is always +1, so no frame mismatch handling needed here.

    let (qe1, qe2) = (hsp1.prelim_qe as isize, hsp2.prelim_qe as isize);
    if contained_in_hsp(
        q1,
        qe1,
        q2,
        hsp1.prelim_ss as isize,
        hsp1.prelim_se as isize,
        hsp2.prelim_ss as isize,
    ) || contained_in_hsp(
        q1,
        qe1,
        qe2,
        hsp1.prelim_ss as isize,
        hsp1.prelim_se as isize,
        hsp2.prelim_se as isize,
    ) {
        let len1 = (qe1 - q1).max(0) as f64;
        let len2 = (qe2 - q2).max(0) as f64;
        let score_density = (hsp1.prelim_score as f64 + hsp2.prelim_score as f64) / (len1 + len2);

        hsp1.prelim_qs = hsp1.prelim_qs.min(hsp2.prelim_qs);
        hsp1.prelim_qe = hsp1.prelim_qe.max(hsp2.prelim_qe);
        hsp1.prelim_ss = hsp1.prelim_ss.min(hsp2.prelim_ss);
        hsp1.prelim_se = hsp1.prelim_se.max(hsp2.prelim_se);

        if hsp2.prelim_score > hsp1.prelim_score {
            hsp1.seed_qs = hsp2.seed_qs;
            hsp1.seed_ss = hsp2.seed_ss;
            hsp1.prelim_score = hsp2.prelim_score;
        }

        let new_len = (hsp1.prelim_qe - hsp1.prelim_qs).max(0) as f64;
        let score_from_density = (score_density * new_len) as i32;
        if score_from_density > hsp1.prelim_score {
            hsp1.prelim_score = score_from_density;
        }
        return true;
    }

    false
}

/// NCBI's preliminary subject-best-hit filter (`-subject_besthit`) over the combined HSP
/// list of a subject, applied after each chunk is merged, before traceback: the list holds the
/// HSPs of every query, each with the length of its own context.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:586-591
fn prelim_subject_best_hit(combined: &mut Vec<PrelimHit>, query_lengths: &[usize]) {
    crate::core::blast_hits::subject_besthit(combined, true, |hit| {
        crate::core::blast_hits::BestHitExtent {
            context: i64::from(hit.context_idx),
            query_frame: hit.query_frame,
            offset: i64::from(hit.prelim_qs),
            end: i64::from(hit.prelim_qe),
            query_length: query_lengths[hit.query_idx as usize] as i64,
        }
    });
}

/// Where the second of two preliminary HSP lists continues the first: the start of a subject
/// chunk, or the offsets of a query chunk's contexts (plus and minus strand; `Blast_HSPListsMerge`
/// takes a negative `contexts_per_query` for the subject).
#[derive(Clone, Copy)]
enum HspListSplit {
    Subject { offset: usize },
    Query { offsets: [i32; 2] },
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2857-2862
// `hsp_num_max` is `INT4_MAX` for BLASTN (`BlastHspNumMax`), so every HSP is kept.
fn merge_prelim_hit_lists(
    combined: &mut Vec<PrelimHit>,
    mut incoming: Vec<PrelimHit>,
    split: HspListSplit,
    chunk_overlap_size: usize,
    allow_gap: bool,
) -> Vec<PrelimHit> {
    if incoming.is_empty() {
        return incoming;
    }
    if combined.is_empty() {
        combined.append(&mut incoming);
        return incoming;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2890-2945
    // The swaps put the HSPs of the overlap strip first, in order, and reorder the others;
    // the stable sort by score below keeps that order for equal HSPs.
    let overlap = chunk_overlap_size as i64;
    let in_strip = |hsp: &PrelimHit, combined_side: bool| match split {
        HspListSplit::Subject { offset } => {
            if combined_side {
                hsp.prelim_se > offset
            } else {
                hsp.prelim_ss < offset.saturating_add(chunk_overlap_size)
            }
        }
        HspListSplit::Query { offsets } => {
            let split_offset = offsets[(hsp.context_idx % 2) as usize] as i64;
            if split_offset < 0 {
                return false;
            }
            let ends_after = hsp.prelim_qe as i64 > split_offset;
            let starts_before = (hsp.prelim_qs as i64) < split_offset + overlap;
            if (hsp.query_frame >= 0) == combined_side {
                ends_after
            } else {
                starts_before
            }
        }
    };
    let mut hspcnt1 = 0usize;
    for index1 in 0..combined.len() {
        if in_strip(&combined[index1], true) {
            combined.swap(hspcnt1, index1);
            hspcnt1 += 1;
        }
    }
    let mut hspcnt2 = 0usize;
    for index2 in 0..incoming.len() {
        if in_strip(&incoming[index2], false) {
            incoming.swap(hspcnt2, index2);
            hspcnt2 += 1;
        }
    }

    if hspcnt1 > 0 && hspcnt2 > 0 {
        let mut deleted = vec![false; incoming.len()];
        for index1 in 0..hspcnt1 {
            let hsp1_context = combined[index1].context_idx;
            for index2 in 0..hspcnt2 {
                if deleted[index2] || incoming[index2].context_idx != hsp1_context {
                    continue;
                }
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2982-2990
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1465-1477
                // `hsp1` may already have been expanded by an earlier merge in this
                // inner loop, so recompute its diagonals every iteration.
                // NCBI's `Int4` query offsets (`PrelimHit::prelim_qs`, `prelim_qe`).
                let end_diag = |hsp: &PrelimHit| hsp.prelim_qe as i64 - hsp.prelim_se as i64;
                let start_diag = |hsp: &PrelimHit| hsp.prelim_qs as i64 - hsp.prelim_ss as i64;
                let hsp1_first = matches!(split, HspListSplit::Subject { .. })
                    || combined[index1].query_frame >= 0;
                let (end, start) = if hsp1_first {
                    (end_diag(&combined[index1]), start_diag(&incoming[index2]))
                } else {
                    (end_diag(&incoming[index2]), start_diag(&combined[index1]))
                };
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2991-2995
                if (end - start).abs() < OVERLAP_DIAG_CLOSE as i64
                    && merge_two_prelim_hits(&mut combined[index1], &incoming[index2], allow_gap)
                {
                    deleted[index2] = true;
                }
            }
        }

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3000-3001
        let mut delete_idx = 0usize;
        incoming.retain(|_| {
            let keep = !deleted[delete_idx];
            delete_idx += 1;
            keep
        });
    }

    append_prelim_hit_list(combined, incoming)
}

/// NCBI's `s_BlastHSPListsCombineByScore` when every HSP is kept (`hsp_num_max` is
/// `INT4_MAX`): the HSPs of the second list follow those of the first, then the list is sorted
/// by score. `Blast_HSPListAppend` is this alone.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2757-2766
fn append_prelim_hit_list(
    combined: &mut Vec<PrelimHit>,
    mut incoming: Vec<PrelimHit>,
) -> Vec<PrelimHit> {
    if incoming.is_empty() {
        return incoming;
    }
    if combined.is_empty() {
        combined.append(&mut incoming);
        return incoming;
    }
    combined.append(&mut incoming);
    sort_prelim_hits_by_score(combined);
    incoming
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1414-1435
fn evalue_compare_prelim_hits(a: &PrelimHit, b: &PrelimHit) -> std::cmp::Ordering {
    evalue_comp(a.prelim_evalue, b.prelim_evalue).then_with(|| score_compare_prelim_hits(a, b))
}

/// A subject's preliminary HSPs of one query: NCBI's `BlastHSPList` in the HSP stream of the
/// preliminary stage, kept in the query's hit list of `prelim_hitlist_size` subjects.
struct PrelimHspList {
    oid: u32,
    hsps: Vec<PrelimHit>,
    best_evalue: f64,
}

impl HitListEntry for PrelimHspList {
    fn oid(&self) -> u32 {
        self.oid
    }

    fn hsp_count(&self) -> usize {
        self.hsps.len()
    }

    fn best_evalue(&self) -> f64 {
        self.best_evalue
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1740-1748
    fn update_best_evalue(&mut self) {
        let mut best = i32::MAX as f64;
        for hsp in &self.hsps {
            if hsp.prelim_evalue < best {
                best = hsp.prelim_evalue;
            }
        }
        self.best_evalue = best;
    }

    fn first_score(&self) -> Option<i32> {
        self.hsps.first().map(|hsp| hsp.prelim_score)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1442-1455
    fn sort_by_evalue(&mut self) {
        let sorted = self.hsps.windows(2).all(|pair| {
            evalue_compare_prelim_hits(&pair[0], &pair[1]) != std::cmp::Ordering::Greater
        });
        if !sorted {
            qsort_prelim_hits_by(&mut self.hsps, evalue_compare_prelim_hits);
        }
    }
}

/// The hit lists of the preliminary stage, one per query.
type PrelimHitLists = Vec<Option<HitList<PrelimHspList>>>;

/// NCBI's HSP collector of the preliminary stage: each subject's HSPs (in subject order) are
/// split into one list per query, in their order, and saved in the query's hit list, which
/// holds at most `prelim_hitlist_size` subjects (`Blast_HitListUpdate` drops the list with
/// the worst preliminary e-value).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:116-161
#[cfg(test)]
fn collect_prelim_hit_lists(
    subject_hits: Vec<Vec<PrelimHit>>,
    num_queries: usize,
    prelim_hitlist_size: usize,
) -> PrelimHitLists {
    let mut collector = PrelimHitCollector::new(num_queries, prelim_hitlist_size);
    for (oid, hits) in subject_hits.into_iter().enumerate() {
        collector.add(oid, hits);
    }
    collector.hit_lists
}

/// EXPERIMENT (LOSAT_X_BNSTREAM, a strict-set switch: on unless it is `0`)
/// : with two or more threads, the subjects of a
/// query batch are claimed by the threads in OID order and their results go to the hit lists in
/// OID order while the threads search the next subjects (`x_bnstream::search_in_order`); a
/// thread claims a subject only while it is fewer than `X_BNSTREAM_SUBJECTS_PER_THREAD`
/// subjects per thread past the last one handed over, so that the HSPs of that window are held
/// at a time instead of those of every subject. The hit lists are the same.
// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes.
fn x_bnstream() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNSTREAM"))
}

/// No NCBI counterpart: the subjects per thread of the window of `LOSAT_X_BNSTREAM`.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
const X_BNSTREAM_SUBJECTS_PER_THREAD: usize = 64;

/// EXPERIMENT (LOSAT_X_BNHITGROUP, a strict-set switch: on unless it is `0`): a subject's HSPs
/// go to lists kept by the collector and only the queries of those HSPs are visited, instead of a
/// new list per query of the batch, each visited: that cost a 24-byte cell and a visit per query
/// for every subject with an HSP, on the thread that takes the subjects in OID order (reads 1e6 x
/// 16S, 33,000 queries per batch x 28,036 subjects: 53 % of the cycles at 8 threads). The hit
/// lists are the same.
// No NCBI counterpart: reads the switch once; it does not change any value NCBI computes. NCBI
// visits every query per subject (hspfilter_collector.c, above); AGENTS.md rule 9 allows less work.
fn x_bnhitgroup() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BNHITGROUP"))
}

/// The preliminary hit lists of a batch filled one subject at a time, in OID order
/// (`collect_prelim_hit_lists`), so that only the HSP lists the hit lists keep stay in memory
/// instead of those of every subject.
struct PrelimHitCollector {
    hit_lists: PrelimHitLists,
    prelim_hitlist_size: usize,
    /// LOSAT_X_BNHITGROUP: one empty list per query from the first subject with an HSP, and the
    /// queries of the subject being added.
    x_slots: Vec<Vec<PrelimHit>>,
    x_touched: Vec<u32>,
    /// LOSAT_X_BNHITGROUP: the HSPs of each query in the subject being added (0 between calls).
    x_counts: Vec<usize>,
}

impl PrelimHitCollector {
    fn new(num_queries: usize, prelim_hitlist_size: usize) -> Self {
        let mut hit_lists: PrelimHitLists = Vec::with_capacity(num_queries);
        hit_lists.resize_with(num_queries, || None);
        Self {
            hit_lists,
            prelim_hitlist_size,
            x_slots: Vec::new(),
            x_touched: Vec::new(),
            x_counts: Vec::new(),
        }
    }

    /// The HSPs of subject `oid`; the subjects must come in OID order.
    fn add(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        if hits.is_empty() {
            return;
        }
        if x_bnhitgroup() {
            self.x_add_grouped(oid, hits);
        } else {
            self.add_per_query(oid, hits);
        }
    }

    /// `LOSAT_X_BNHITGROUP=0`: a list per query of the batch, as NCBI's `hsp_list_array`.
    fn add_per_query(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        let mut per_query: Vec<Vec<PrelimHit>> = vec![Vec::new(); self.hit_lists.len()];
        for hit in hits {
            per_query[hit.query_idx as usize].push(hit);
        }
        for (query, hsps) in per_query.into_iter().enumerate() {
            if !hsps.is_empty() {
                self.update(query, oid, hsps);
            }
        }
    }

    /// EXPERIMENT (LOSAT_X_BNHITGROUP): the same lists, made in the collector's slots: only the
    /// queries of the subject's HSPs are visited. A query's hit list does not read another's, so
    /// that the queries may come in the order of their first HSPs; the HSPs of a query keep their
    /// order. A slot gives its list away, so that the lists grow as those of `add_per_query`.
    fn x_add_grouped(&mut self, oid: usize, hits: Vec<PrelimHit>) {
        if self.x_slots.is_empty() {
            self.x_slots.resize_with(self.hit_lists.len(), Vec::new);
            self.x_counts.resize(self.hit_lists.len(), 0);
        }
        let mut touched = std::mem::take(&mut self.x_touched);
        // No NCBI counterpart: the queries' lists are made at their sizes, so that no HSP is
        // copied as a list grows (a slot gives its list away, so each subject starts it empty).
        for hit in &hits {
            let count = &mut self.x_counts[hit.query_idx as usize];
            if *count == 0 {
                touched.push(hit.query_idx);
            }
            *count += 1;
        }
        for &query in &touched {
            let count = std::mem::take(&mut self.x_counts[query as usize]);
            self.x_slots[query as usize].reserve_exact(count);
        }
        for hit in hits {
            self.x_slots[hit.query_idx as usize].push(hit);
        }
        for &query in &touched {
            let hsps = std::mem::take(&mut self.x_slots[query as usize]);
            self.update(query as usize, oid, hsps);
        }
        touched.clear();
        self.x_touched = touched;
    }

    fn update(&mut self, query: usize, oid: usize, hsps: Vec<PrelimHit>) {
        let prelim_hitlist_size = self.prelim_hitlist_size;
        self.hit_lists[query]
            .get_or_insert_with(|| HitList::new(prelim_hitlist_size))
            .update(PrelimHspList {
                oid: oid as u32,
                hsps,
                best_evalue: 0.0,
            });
    }
}

/// The preliminary HSPs that the traceback reads, for each subject: the kept lists of the
/// queries, a query after another. NCBI's `BlastHSPStreamClose` sorts the kept lists by
/// subject, and the traceback reads each subject's lists (`BlastHSPStreamBatchRead`); the
/// traceback of one query's list does not depend on the lists of the other queries.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:196-203
fn prelim_hits_by_subject(hit_lists: PrelimHitLists, num_subjects: usize) -> Vec<Vec<PrelimHit>> {
    // No NCBI counterpart: each subject's list is made at its size, so that it is not copied as it
    // grows (with one subject, every kept HSP of the batch was copied again at each doubling).
    let mut sizes = vec![0usize; num_subjects];
    for list in hit_lists
        .iter()
        .flatten()
        .flat_map(|hit_list| &hit_list.hsplist_array)
    {
        sizes[list.oid as usize] += list.hsps.len();
    }
    let mut by_subject: Vec<Vec<PrelimHit>> = sizes.into_iter().map(Vec::with_capacity).collect();
    for hit_list in hit_lists.into_iter().flatten() {
        for list in hit_list.hsplist_array {
            by_subject[list.oid as usize].extend(list.hsps);
        }
    }
    by_subject
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:196-203
/// `prelim_hits_by_subject` for a database: the ordinal ids of the kept lists, ascending,
/// each with its kept HSPs in the same order (the lists of the queries, a query after
/// another), so that the traceback reads only these sequences
/// (`database::traceback_listed`).
fn listed_prelim_hits(hit_lists: PrelimHitLists) -> Vec<(u32, Vec<PrelimHit>)> {
    let mut by_subject: std::collections::BTreeMap<u32, Vec<PrelimHit>> =
        std::collections::BTreeMap::new();
    for hit_list in hit_lists.into_iter().flatten() {
        for list in hit_list.hsplist_array {
            by_subject.entry(list.oid).or_default().extend(list.hsps);
        }
    }
    by_subject.into_iter().collect()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:84-105
struct SubjectScratch {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
    hits: Vec<BlastnHsp>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4077
    prelim_hits: Vec<PrelimHit>,
    ungapped_hits: Vec<UngappedHit>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
    greedy_align_scratch: GreedyAlignScratch,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
    // Every look-ahead helper has its own scratch, as each NCBI call has the `gap_align` structure
    // that holds its working memory. A helper never uses the caller's scratch. Reuse only.
    // No NCBI counterpart: scratch memory of the look-ahead helpers.
    x_ahead_scratch: Vec<(GapAlignScratch, GreedyAlignScratch)>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:509-512
    // No NCBI counterpart: the scratch of the speculative traceback batches, one pair per thread
    // of the pool, kept with the rest of this scratch so that each subject's traceback reuses it
    // instead of making it again. Reuse only, as `x_ahead_scratch`.
    x_spec_scratch: Vec<(GapAlignScratch, GreedyAlignScratch)>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:184-198
    seq_ranges_scratch: Vec<(i32, i32)>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321,399-405
    // The buffers of a database search with `-db_soft_mask` or `-db_hard_mask` (decision D4 of
    // M4a-2): the masks of the block of the database sequence being searched (read once for all
    // the block's sequences, review 1 F3 of M4a-2a), its soft ranges (the complement of its masks;
    // with `-db_hard_mask` its one full range, NCBI's `&(backup->full_range)`) and its hard ranges
    // (the complement, `-db_hard_mask` only, with the state of its chunks: `HardSplitState`, whose
    // ranges are empty outside a hard-masked sequence's preliminary search; M4a-2b), kept from one
    // sequence to the next as NCBI keeps one `TSequenceRanges` per thread. Empty, and never
    // touched, without the options.
    db_masks: crate::seqdb::MaskBlock,
    db_soft_ranges: Vec<(i32, i32)>,
    db_hard_split: HardSplitState,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:77-95
    hit_level_array: Vec<DiagStruct>,
    hit_len_array: Vec<u8>,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:97-106
    diag_hash: DiagHashTable,
    // EXPERIMENT (LOSAT_X_PAIRPARFOLD, `x_pair_seed_part.rs`): the allocations of the partition
    // tables and buffers of the partitioned fold of the strips paths. A partitioned stage on the
    // chains copies `diag_hash` into the partition tables, so no table content carries over from
    // one stage to the next; with LOSAT_X_DCLATEST the stages on the maps copy no table, and the
    // maps (each partition's and the owner's, `PartState::latest`, `PartScratch::latest`) carry
    // their entries from stage to stage, dead at the next stage's start (`DiagLatest::begin`).
    // Empty while the switch is off. No NCBI counterpart.
    x_pair_fold: x_pair_seed::PartScratch,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:159-176
    diag_table_offset: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    offset_pairs: Vec<OffsetPair>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_itree.c:146-170
    // The interval tree of the subject's preliminary gapped stage and of its traceback, kept with
    // this scratch and given each subject's ranges with `reset_with_bounds` (the same state as a
    // new tree) instead of being made again for every subject. Reuse only.
    x_tree: BlastIntervalTree,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:488-491
impl SubjectScratch {
    fn new(offset_array_size: usize, window_size: usize) -> Self {
        Self {
            hits: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4012-4077
            prelim_hits: Vec::new(),
            ungapped_hits: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
            greedy_align_scratch: GreedyAlignScratch::new(),
            // No NCBI counterpart: the helper scratch list starts empty and is filled when helpers are used;
            // it does not change any value NCBI computes.
            x_ahead_scratch: Vec::new(),
            x_spec_scratch: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:184-198
            seq_ranges_scratch: Vec::new(),
            db_masks: crate::seqdb::MaskBlock::default(),
            db_soft_ranges: Vec::new(),
            db_hard_split: HardSplitState::default(),
            hit_level_array: Vec::new(),
            hit_len_array: Vec::new(),
            diag_hash: DiagHashTable::new(window_size as i32),
            // No NCBI counterpart: empty until a partitioned stage runs (LOSAT_X_PAIRPARFOLD).
            x_pair_fold: x_pair_seed::PartScratch::default(),
            diag_table_offset: window_size as i32,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
            offset_pairs: vec![
                OffsetPair {
                    q_off: 0,
                    s_off: 0,
                    s_range: 0,
                };
                offset_array_size
            ],
            // No NCBI counterpart: every user sets the subject's ranges first (`reset_with_bounds`).
            x_tree: BlastIntervalTree::new(0, 1, 0, 1),
        }
    }
}

#[derive(Clone)]
struct QueryContext {
    query_idx: u32,
    frame: i32,
    query_offset: i32,
    /// The context's letter count. The letters are held once, in BLASTNA, in the query buffer.
    len: usize,
    masks: Vec<MaskedInterval>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2424-2429
/// NCBI's query of context `ctx` (`query_tmp`, `ContextQuery`) in the concatenated query
/// `concat`, NCBI's `query->sequence` without its first sentinel (`encoded_query_concat_blastna`
/// or the same letters without masks).
#[inline(always)]
fn context_query<'a>(concat: &'a [u8], ctx: &QueryContext) -> ContextQuery<'a> {
    ContextQuery {
        buf: concat,
        start: ctx.query_offset.max(0) as usize,
        len: ctx.len as i32,
    }
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:503-512
/// The result of the traceback call of NCBI, of the dynamic-programming kind or the greedy
/// kind, computed before the ordered loop reaches the HSP.
///
/// No NCBI counterpart: a traceback computed ahead of the ordered
/// loop, by the dynamic-programming kernel or by the greedy one.
#[derive(PartialEq)]
enum XSpecTraceback<D, G> {
    Dp(D),
    Greedy(G),
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:356-357
/// The working memory of one gapped alignment call (the `gap_align` of NCBI, and the greedy aligner's
/// auxiliary memory), borrowed from the caller or from a helper.
///
/// No NCBI counterpart: the scratch memory one gapped alignment may use.
struct XAheadScratch<'a> {
    gap: &'a mut GapAlignScratch,
    greedy: &'a mut GreedyAlignScratch,
}

#[path = "x_ctx_o1.rs"]
mod x_ctx_o1;

struct QueryContextIndex {
    offsets: Vec<usize>,
    min_length: usize,
    max_length: usize,
    // Direct mapping from query offset to context index for fast lookup.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243 (BSearchContextInfo)
    direct_map: Vec<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
    // NCBI finds the context of a query offset by this binary search over the `query_offset` of the
    // contexts. `direct_map` is a table of its answer for every offset. The runs below hold the same
    // answer as intervals: for an offset inside the sequence of context `c` it is `c`, else 0.
    // No NCBI counterpart: what `direct_map` holds, as the list of
    // its non-zero runs `(start, end, context)` sorted by `start`, plus the
    // length it would have. Four bytes per query base become a few words.
    x_ranges: Vec<(usize, usize, u32)>,
    x_total_len: usize,
    // LOSAT_X_CTXO1 (graduated 2026-10-27, always on): what `direct_map` holds, in constant time
    // from the layout of the contexts (`x_ctx_o1.rs`). When it is built it answers the offsets the table covers, the runs
    // only those of its crowded blocks, and the table is not built.
    x_o1: Option<x_ctx_o1::CtxO1>,
}

impl QueryContextIndex {
    fn new(contexts: &[QueryContext], one_strand: bool) -> Self {
        let offsets: Vec<usize> = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize)
            .collect();
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:229-243
        // With one strand searched `min_length` is 0, so `bsearch_context_info` searches all
        // contexts (the bounds of the other case assume every context has a length).
        let min_length = if one_strand {
            0
        } else {
            contexts
                .iter()
                .map(|ctx| ctx.len)
                .filter(|&len| len > 0)
                .min()
                .unwrap_or(0)
        };
        let max_length = contexts.iter().map(|ctx| ctx.len).max().unwrap_or(0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243 (BSearchContextInfo)
        let total_len = contexts
            .iter()
            .map(|ctx| ctx.query_offset.max(0) as usize + ctx.len)
            .max()
            .unwrap_or(0);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-227,235-241
        // Builds the runs (start, end, context) from the same context offsets and lengths the table is
        // filled from. If two sequences would overlap, the table's "later context wins" rule would not be
        // a list of disjoint runs, so the runs are dropped and the table is used.
        // No NCBI counterpart: the runs the loop below would fill.
        // Later contexts overwrite earlier ones in `direct_map`, so the list
        // only stands in for it when no two runs overlap; context 0 needs no
        // run because the map starts out as zeros.
        // `x_mode` 1: the runs answer and the table is not built; 0: the table answers (overlapping
        // runs).
        let mut x_mode = 1u8;
        // LOSAT_X_CTXO1 (graduated 2026-10-27, always on): the constant-time form, where the layout
        // allows it (`build` answers `None` otherwise: a rule of the layout and not a switch). It
        // takes the place of the runs (kept only for the offsets of its crowded blocks) and of the
        // table.
        let x_o1 = x_ctx_o1::CtxO1::build(contexts, total_len);
        // the runs answer what the form leaves (offsets in a crowded block of `x_ctx_o1.rs`)
        let x_o1_alone = x_o1.as_ref().is_some_and(|o1| !o1.needs_runs());
        let mut x_ranges: Vec<(usize, usize, u32)> = Vec::new();
        if x_mode != 0 && !x_o1_alone {
            for (idx, ctx) in contexts.iter().enumerate() {
                let start = ctx.query_offset.max(0) as usize;
                let end = start.saturating_add(ctx.len);
                if start < end && end <= total_len && idx != 0 {
                    x_ranges.push((start, end, idx as u32));
                }
            }
            let mut all: Vec<(usize, usize)> = contexts
                .iter()
                .map(|ctx| {
                    let start = ctx.query_offset.max(0) as usize;
                    (start, start.saturating_add(ctx.len))
                })
                .filter(|&(start, end)| start < end && end <= total_len)
                .collect();
            all.sort_unstable();
            if all.windows(2).any(|w| w[0].1 > w[1].0) {
                x_mode = 0;
                x_ranges.clear();
            }
            x_ranges.sort_unstable();
        }
        let mut direct_map = Vec::new();
        if x_mode != 1 {
            direct_map = vec![0u32; total_len];
            for (idx, ctx) in contexts.iter().enumerate() {
                let start = ctx.query_offset.max(0) as usize;
                let end = start.saturating_add(ctx.len);
                if start < end && end <= direct_map.len() {
                    direct_map[start..end].fill(idx as u32);
                }
            }
        }
        Self {
            offsets,
            min_length,
            max_length,
            direct_map,
            x_ranges,
            x_total_len: if x_mode == 0 || x_o1_alone {
                0
            } else {
                total_len
            },
            x_o1,
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:235-242
    /// The context of offset `n` from the runs: the last run that starts at or before `n`, if `n` is
    /// inside it. This is the value the table holds at `n`.
    ///
    /// `direct_map[n]` for `n < x_total_len`, from the runs.
    #[inline]
    fn x_direct(&self, n: usize) -> usize {
        // the last run that starts at or before `n`
        let after = self.x_ranges.partition_point(|&(start, _, _)| start <= n);
        if after == 0 {
            return 0;
        }
        let (_, end, context) = self.x_ranges[after - 1];
        if n < end {
            context as usize
        } else {
            0
        }
    }

    /// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219
    /// The context of query offset `n`. LOSAT_X_CTXO1 (graduated 2026-10-27): the constant-time form answers
    /// here, inlined into the seed loops; any other form is the out-of-line `x_context_for_offset_ref`.
    #[inline(always)]
    fn context_for_offset(&self, n: usize) -> usize {
        if let Some(o1) = &self.x_o1 {
            if let Some(context) = o1.get(n) {
                return context;
            }
        }
        self.x_context_for_offset_ref(n)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243 (BSearchContextInfo)
    // The forms besides the constant-time one (LOSAT_X_CTXO1), out of line as `context_for_offset` was.
    #[inline(never)]
    fn x_context_for_offset_ref(&self, n: usize) -> usize {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-243
        // The code after this block (the table, then `bsearch_context_info`) ports `BSearchContextInfo`.
        // With the runs the same context is found without the table.
        if n < self.x_total_len {
            return self.x_direct(n);
        }
        if n < self.direct_map.len() {
            return self.direct_map[n] as usize;
        }
        bsearch_context_info(n, self)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-242
fn bsearch_context_info(n: usize, index: &QueryContextIndex) -> usize {
    let size = index.offsets.len();
    if size == 0 {
        return 0;
    }
    let mut b = 0usize;
    let mut e = size;
    if index.min_length > 0 && index.max_length > 0 {
        b = (n / (index.max_length + 1)).min(size - 1);
        e = (n / (index.min_length + 1) + 1).min(size);
    }
    while b + 1 < e {
        let m = (b + e) / 2;
        if index.offsets[m] > n {
            e = m;
        } else {
            b = m;
        }
    }
    b
}

/// The query length that picks the diagonal container (hash above 8000 letters, array
/// otherwise): NCBI's `query_offset + query_length` of the LAST context. With `-strand plus`
/// the last context is the empty minus context of the last query, one past the end of the plus
/// context, so this is the searched block's length + 1; when the last context has letters it
/// is the searched block's length (`query_concat_length`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:225-229
fn diag_container_query_length(query_contexts: &[QueryContext]) -> usize {
    query_contexts
        .last()
        .map_or(0, |last| last.query_offset.max(0) as usize + last.len)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_encoding.c:120-121
/// The query buffer with its sentinels (`query_concat_length + 2` bytes), every context
/// written straight into it in BLASTNA: for query `q` with context offsets
/// `context_offsets[2 * q]` (plus) and `context_offsets[2 * q + 1]` (minus), the plus
/// context is the BLASTNA code of the letters and the minus context that of their reverse
/// complement. The searched block (`query_concat_length` bytes) is `[1..n + 1]` of it
/// (LOSAT's memory: no second copy of the block, no IUPAC copy of either strand). The buffer is
/// `with_sentinels`, the one that the batch loop keeps (`x_query_sentinels`).
fn fill_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    with_sentinels: &mut Vec<u8>,
) {
    debug_assert_eq!(2 * query_letters.count(), context_offsets.len());

    x_query_sentinels(with_sentinels, query_concat_length + 2);

    for (q, offsets) in context_offsets.chunks_exact(2).enumerate() {
        query_letters.read(q, |letters| {
            let length = letters.len();
            // A context that is empty or does not end inside the block stays sentinels.
            // (Twin: `place` and the two encoder calls below are copied in x_query_prep_par.rs
            // (LOSAT_X_QPREPPAR) and, with the strand switches, in
            // `fill_strand_query_blastna_buffer`; change them there too.)
            let place = |offset: i32| -> Option<std::ops::Range<usize>> {
                let start = offset.max(0) as usize;
                let end = start.saturating_add(length);
                (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
            };
            if let Some(range) = place(offsets[0]) {
                encode_iupac_to_blastna_into(letters, &mut with_sentinels[range]);
            }
            if let Some(range) = place(offsets[1]) {
                encode_reverse_complement_blastna_into(letters, &mut with_sentinels[range]);
            }
        });
    }
}

/// `fill_query_blastna_buffer` into a buffer of its own (the tests).
#[cfg(test)]
fn build_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
) -> Vec<u8> {
    let mut with_sentinels = Vec::new();
    fill_query_blastna_buffer(
        query_letters,
        context_offsets,
        query_concat_length,
        &mut with_sentinels,
    );
    with_sentinels
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-603
/// `buffer` as `length` sentinels: the query buffer of a batch before its contexts are written.
/// NCBI allocates the buffer of each batch; LOSAT's batch loop keeps one buffer from batch to
/// batch and each batch writes all of it again (LOSAT's memory). A buffer freed by each batch
/// leaves a hole of its size under the blocks allocated after it, and glibc raises its mmap
/// threshold to the size of the freed buffer: the next batch's buffer, a few KiB larger, did not
/// fit the hole once small blocks had split it, was mapped on its own, and the 8.5 MiB hole stayed
/// resident (31,168 cDNAs x 3 Mb, dc-megablast `-db`, 1 thread: +9.5 MiB in half the runs, set off
/// by the length of the output file's name). The kept buffer has the capacity of
/// `x_mapped_capacity`, so from 1 MiB glibc maps it and its free leaves the threshold as it is.
/// A batch writes all of its bytes again while the kept buffer is long enough and no longer than
/// the larger of the batch's own capacity and that of 1 MiB (glibc's 32 MiB ceiling and a page;
/// 1 MiB elsewhere), and otherwise takes a new block: a batch below 1 MiB writes the mapped
/// buffer of a longer batch again, while a batch shorter than one past the ceiling does not keep
/// that batch's buffer. So the kept pages beyond a batch's length are at most 32 MiB, and pages
/// fault once, when a batch first reaches them. (Taking a new block for each batch below 1 MiB
/// mapped the buffer again for each longer batch after it: the cDNAs above, 16 threads: 18 maps of
/// 32 MiB instead of 1, +40,000 minor faults, +5 % wall.) Bytes past the batch's length are not
/// read.
fn x_query_sentinels(buffer: &mut Vec<u8>, length: usize) {
    use crate::utils::x_mapped::{x_mapped_capacity, X_MAPPED_FROM};
    let capacity = x_mapped_capacity(length, 1);
    if buffer.capacity() < length
        || buffer.capacity() > capacity.max(x_mapped_capacity(X_MAPPED_FROM, 1))
    {
        // The old contents are not kept: the old block is freed first, then a new one, not a copy.
        drop(std::mem::take(buffer));
        *buffer = Vec::with_capacity(capacity);
    } else {
        buffer.clear();
    }
    buffer.resize(length, NUCL_SENTINEL);
}

/// `-soft_masking false`: writes the masks of the contexts into the searched BLASTNA copy of the
/// query (the buffer with its sentinels) and returns the copy without them, or `None` without
/// masks. Identities are counted on the returned letters, everything else reads the buffer,
/// where a masked letter is an N like one of the input: its word is not in the lookup table, the
/// extensions and the Karlin blocks read N, and the small table's compressed query reads it as
/// a base (`build_compressed_query`).
///
/// A context's own list of masks is cleared: of the lookup tables that keep `masked_locations`
/// for `s_TypeOfWord`, only the small one does so without `mask_at_hash`, and the seed loops take
/// that table's checks from `small_na_word`. `query_context_masks` (the lookup table's words)
/// stays.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1369-1405
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:407-411,1269-1272
#[cold]
#[inline(never)]
fn mask_query_letters(
    buffer: &mut [u8],
    contexts: &mut [QueryContext],
    context_masks: &[Vec<MaskedInterval>],
) -> Option<Vec<u8>> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:38
    const BLASTNA_N: u8 = 14;
    if context_masks.iter().all(Vec::is_empty) {
        return None;
    }
    let nomask = buffer.to_vec();
    for (ctx, masks) in contexts.iter_mut().zip(context_masks) {
        let first = ctx.query_offset.max(0) as usize + 1;
        for mask in masks {
            buffer[first + mask.start..first + mask.end].fill(BLASTNA_N);
        }
        ctx.masks.clear();
    }
    Some(nomask)
}

/// The query buffer of a search of one strand (`-strand plus` or `minus`): as
/// `fill_query_blastna_buffer`, but a context is written only when its length is above 0,
/// with its own length (the strand that is not searched has a context of length 0 and takes no
/// room, so it must not be written over the next query).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:71-81
fn fill_strand_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    strand: QueryStrand,
    with_sentinels: &mut Vec<u8>,
) {
    debug_assert_eq!(2 * query_letters.count(), context_offsets.len());
    let (plus_on, minus_on) = strand.searched();

    x_query_sentinels(with_sentinels, query_concat_length + 2);

    for (q, offsets) in context_offsets.chunks_exact(2).enumerate() {
        query_letters.read(q, |letters| {
            let place = |offset: i32, length: usize| -> Option<std::ops::Range<usize>> {
                let start = offset.max(0) as usize;
                let end = start.saturating_add(length);
                (start < end && end <= query_concat_length).then(|| start + 1..end + 1)
            };
            if plus_on {
                if let Some(range) = place(offsets[0], letters.len()) {
                    encode_iupac_to_blastna_into(letters, &mut with_sentinels[range]);
                }
            }
            if minus_on {
                if let Some(range) = place(offsets[1], letters.len()) {
                    encode_reverse_complement_blastna_into(letters, &mut with_sentinels[range]);
                }
            }
        });
    }
}

/// `fill_strand_query_blastna_buffer` into a buffer of its own (the tests).
#[cfg(test)]
fn build_strand_query_blastna_buffer<L: QueryLetters + ?Sized>(
    query_letters: &L,
    context_offsets: &[i32],
    query_concat_length: usize,
    strand: QueryStrand,
) -> Vec<u8> {
    let mut with_sentinels = Vec::new();
    fill_strand_query_blastna_buffer(
        query_letters,
        context_offsets,
        query_concat_length,
        strand,
        &mut with_sentinels,
    );
    with_sentinels
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_util.c:190-203
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:402-406
fn compute_lookup_query_stats(
    contexts: &[QueryContext],
    context_masks: &[Vec<MaskedInterval>],
) -> (usize, usize) {
    debug_assert_eq!(contexts.len(), context_masks.len());
    let mut approx_entries = 0i64;
    let mut max_q_off = 0usize;

    // NCBI computes the lookup segments (BLAST_ComplementMaskLocations, which skips the
    // contexts marked invalid) before its ungapped blocks mark contexts invalid
    // (blast_setup.c:633-653), so every context of a query with letters counts.
    for (ctx, masks) in contexts.iter().zip(context_masks.iter()) {
        if ctx.len == 0 {
            continue;
        }
        let ctx_offset = ctx.query_offset.max(0) as usize;
        let ranges = build_unmasked_ranges(ctx.len, masks);
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1081-1115
        // A context masked from its first letter to its last gets the range
        // (end_offset + 1, end_offset): one entry less, and its end in `max_off`.
        if ranges.is_empty() {
            approx_entries -= 1;
            max_q_off = max_q_off.max(ctx_offset + ctx.len - 1);
            continue;
        }
        for (start, end) in ranges {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_util.c:190-203
            // NCBI's `right` is inclusive (the last unmasked letter), LOSAT's `end` is not.
            approx_entries += (end - start) as i64 - 1;
            max_q_off = max_q_off.max(ctx_offset + end - 1);
        }
    }

    // A negative estimate (every context masked) selects the smallest table, as 0 does.
    (approx_entries.max(0) as usize, max_q_off)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/unit_tests/api/ntscan_unit_test.cpp:739-777
// NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1647-1674
// The ranges of a soft-masked subject (bl2seq `-lcase_masking`) come from the inclusive
// (from, to) pairs of its lowercase masks.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:700-707
// The pairs are stored from `_data[1]`, but `get_data` returns `_data`, which
// `SetupSubjects_OMF` passes as `size + 1` ranges.
// NCBI reference (598d8ae6): c++/include/objtools/blast/seqdb_reader/seqdb.hpp:280-282
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:819-826
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:211-213
// So range k is {`to` of mask k-1, `from` of mask k}: it starts at the last lowercase
// letter of the previous mask (a word may start there) and ends before the next mask.
fn build_subject_seq_ranges_from_masks(
    masks: &[MaskedInterval],
    subject_len: usize,
) -> Vec<(i32, i32)> {
    let mut ranges = Vec::with_capacity(masks.len() + 1);
    let mut left = 0i32;
    for mask in masks.iter().filter(|mask| mask.start < subject_len) {
        ranges.push((left, mask.start as i32));
        left = (mask.end.min(subject_len) - 1) as i32;
    }
    ranges.push((left, subject_len as i32));
    ranges
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2049-2067
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-892
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:493-555
struct BlastnTiming {
    scan_ns: std::sync::atomic::AtomicU64,
    scan_calls: std::sync::atomic::AtomicU64,
    ungapped_ns: std::sync::atomic::AtomicU64,
    ungapped_calls: std::sync::atomic::AtomicU64,
    gapped_ns: std::sync::atomic::AtomicU64,
    gapped_calls: std::sync::atomic::AtomicU64,
    traceback_ns: std::sync::atomic::AtomicU64,
    traceback_calls: std::sync::atomic::AtomicU64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:358-613
    traceback_sort_prelim_ns: std::sync::atomic::AtomicU64,
    traceback_tree_precheck_ns: std::sync::atomic::AtomicU64,
    traceback_start_offsets_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_dp_ns: std::sync::atomic::AtomicU64,
    traceback_alignment_greedy_ns: std::sync::atomic::AtomicU64,
    traceback_hsp_build_ns: std::sync::atomic::AtomicU64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:633-692
    purge_common_endpoint_pass1_ns: std::sync::atomic::AtomicU64,
    reevaluate_ambiguities_ns: std::sync::atomic::AtomicU64,
    identity_length_test_ns: std::sync::atomic::AtomicU64,
    purge_common_endpoint_pass2_ns: std::sync::atomic::AtomicU64,
    traceback_score_resort_ns: std::sync::atomic::AtomicU64,
    traceback_final_tree_contains_ns: std::sync::atomic::AtomicU64,
    traceback_final_tree_add_ns: std::sync::atomic::AtomicU64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:836-892
    traceback_postprocess_ns: std::sync::atomic::AtomicU64,
    traceback_prelim_hsps: std::sync::atomic::AtomicU64,
    traceback_tree_precheck_skipped_hsps: std::sync::atomic::AtomicU64,
    traceback_full_traceback_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_evalue_cutoff_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_reevaluation_hsps: std::sync::atomic::AtomicU64,
    traceback_deleted_identity_length_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_endpoint_pass1_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_endpoint_pass2_hsps: std::sync::atomic::AtomicU64,
    traceback_removed_final_tree_hsps: std::sync::atomic::AtomicU64,
    traceback_edit_script_lengths: std::sync::Mutex<Vec<u32>>,
    traceback_alignment_lengths: std::sync::Mutex<Vec<u32>>,
    format_ns: std::sync::atomic::AtomicU64,
    format_calls: std::sync::atomic::AtomicU64,
}

impl BlastnTiming {
    fn new() -> Self {
        Self {
            scan_ns: std::sync::atomic::AtomicU64::new(0),
            scan_calls: std::sync::atomic::AtomicU64::new(0),
            ungapped_ns: std::sync::atomic::AtomicU64::new(0),
            ungapped_calls: std::sync::atomic::AtomicU64::new(0),
            gapped_ns: std::sync::atomic::AtomicU64::new(0),
            gapped_calls: std::sync::atomic::AtomicU64::new(0),
            traceback_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_calls: std::sync::atomic::AtomicU64::new(0),
            traceback_sort_prelim_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_tree_precheck_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_start_offsets_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_dp_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_alignment_greedy_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_hsp_build_ns: std::sync::atomic::AtomicU64::new(0),
            purge_common_endpoint_pass1_ns: std::sync::atomic::AtomicU64::new(0),
            reevaluate_ambiguities_ns: std::sync::atomic::AtomicU64::new(0),
            identity_length_test_ns: std::sync::atomic::AtomicU64::new(0),
            purge_common_endpoint_pass2_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_score_resort_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_final_tree_contains_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_final_tree_add_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_postprocess_ns: std::sync::atomic::AtomicU64::new(0),
            traceback_prelim_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_tree_precheck_skipped_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_full_traceback_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_evalue_cutoff_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_reevaluation_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_deleted_identity_length_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_endpoint_pass1_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_endpoint_pass2_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_removed_final_tree_hsps: std::sync::atomic::AtomicU64::new(0),
            traceback_edit_script_lengths: std::sync::Mutex::new(Vec::new()),
            traceback_alignment_lengths: std::sync::Mutex::new(Vec::new()),
            format_ns: std::sync::atomic::AtomicU64::new(0),
            format_calls: std::sync::atomic::AtomicU64::new(0),
        }
    }

    fn record_duration(counter: &std::sync::atomic::AtomicU64, start: std::time::Instant) {
        counter.fetch_add(
            start.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }

    fn record_ns(counter: &std::sync::atomic::AtomicU64, elapsed_ns: u64) {
        counter.fetch_add(elapsed_ns, std::sync::atomic::Ordering::Relaxed);
    }

    fn record_count(counter: &std::sync::atomic::AtomicU64, count: u64) {
        counter.fetch_add(count, std::sync::atomic::Ordering::Relaxed);
    }

    fn record_traceback_lengths(&self, edit_script_lengths: Vec<u32>, alignment_lengths: Vec<u32>) {
        if !edit_script_lengths.is_empty() {
            if let Ok(mut samples) = self.traceback_edit_script_lengths.lock() {
                samples.extend(edit_script_lengths);
            }
        }
        if !alignment_lengths.is_empty() {
            if let Ok(mut samples) = self.traceback_alignment_lengths.lock() {
                samples.extend(alignment_lengths);
            }
        }
    }
}

fn timing_seconds(counter: &std::sync::atomic::AtomicU64) -> f64 {
    counter.load(std::sync::atomic::Ordering::Relaxed) as f64 / 1e9
}

fn timing_count(counter: &std::sync::atomic::AtomicU64) -> u64 {
    counter.load(std::sync::atomic::Ordering::Relaxed)
}

fn timing_length_stats(samples: &std::sync::Mutex<Vec<u32>>) -> (usize, u32, f64) {
    let Ok(samples) = samples.lock() else {
        return (0, 0, 0.0);
    };
    if samples.is_empty() {
        return (0, 0, 0.0);
    }
    let mut sorted = samples.clone();
    sorted.sort_unstable();
    let count = sorted.len();
    let max = sorted[count - 1];
    let median = if count % 2 == 0 {
        let upper = sorted[count / 2] as f64;
        let lower = sorted[(count / 2) - 1] as f64;
        (lower + upper) / 2.0
    } else {
        sorted[count / 2] as f64
    };
    (count, max, median)
}

// No NCBI counterpart: LOSAT_TIMING lines of one query batch (diagnostics only). The batch
// returns them as text, which `run_in_pool` writes to standard error when it takes the batch.
fn blastn_timing_text(
    timing: &BlastnTiming,
    t_search_start: Option<std::time::Instant>,
    t_total: Option<std::time::Instant>,
) -> String {
    use std::fmt::Write as _;
    let mut text = String::new();
    // One line of the block, as `eprintln!` wrote it before.
    macro_rules! timing_line {
        ($($arg:tt)*) => {{
            let _ = writeln!(text, $($arg)*);
        }};
    }
    let scan_s = timing_seconds(&timing.scan_ns);
    let scan_n = timing_count(&timing.scan_calls);
    let ungapped_s = timing_seconds(&timing.ungapped_ns);
    let ungapped_n = timing_count(&timing.ungapped_calls);
    let gapped_s = timing_seconds(&timing.gapped_ns);
    let gapped_n = timing_count(&timing.gapped_calls);
    let traceback_s = timing_seconds(&timing.traceback_ns);
    let traceback_n = timing_count(&timing.traceback_calls);
    let format_s = timing_seconds(&timing.format_ns);
    let format_n = timing_count(&timing.format_calls);

    timing_line!("[TIMING] scan_lookup: {:.3}s (calls={})", scan_s, scan_n);
    timing_line!(
        "[TIMING] ungapped_extend: {:.3}s (calls={})",
        ungapped_s,
        ungapped_n
    );
    timing_line!(
        "[TIMING] gapped_extend: {:.3}s (calls={})",
        gapped_s,
        gapped_n
    );
    timing_line!(
        "[TIMING] traceback_prune: {:.3}s (calls={})",
        traceback_s,
        traceback_n
    );
    timing_line!(
        "[TIMING] traceback_sort_prelim: {:.3}s",
        timing_seconds(&timing.traceback_sort_prelim_ns)
    );
    timing_line!(
        "[TIMING] traceback_tree_precheck: {:.3}s",
        timing_seconds(&timing.traceback_tree_precheck_ns)
    );
    timing_line!(
        "[TIMING] traceback_start_offsets: {:.3}s",
        timing_seconds(&timing.traceback_start_offsets_ns)
    );
    timing_line!(
        "[TIMING] traceback_alignment: {:.3}s (dp={:.3}s greedy={:.3}s)",
        timing_seconds(&timing.traceback_alignment_ns),
        timing_seconds(&timing.traceback_alignment_dp_ns),
        timing_seconds(&timing.traceback_alignment_greedy_ns)
    );
    timing_line!(
        "[TIMING] traceback_hsp_build: {:.3}s",
        timing_seconds(&timing.traceback_hsp_build_ns)
    );
    timing_line!(
        "[TIMING] purge_common_endpoint_pass1: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass1_ns)
    );
    timing_line!(
        "[TIMING] reevaluate_ambiguities: {:.3}s",
        timing_seconds(&timing.reevaluate_ambiguities_ns)
    );
    timing_line!(
        "[TIMING] identity_length_test: {:.3}s",
        timing_seconds(&timing.identity_length_test_ns)
    );
    timing_line!(
        "[TIMING] purge_common_endpoint_pass2: {:.3}s",
        timing_seconds(&timing.purge_common_endpoint_pass2_ns)
    );
    timing_line!(
        "[TIMING] traceback_score_resort: {:.3}s",
        timing_seconds(&timing.traceback_score_resort_ns)
    );
    timing_line!(
        "[TIMING] traceback_final_tree_contains: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_contains_ns)
    );
    timing_line!(
        "[TIMING] traceback_final_tree_add: {:.3}s",
        timing_seconds(&timing.traceback_final_tree_add_ns)
    );
    timing_line!(
        "[TIMING] traceback_postprocess: {:.3}s",
        timing_seconds(&timing.traceback_postprocess_ns)
    );
    timing_line!(
        "[TIMING] traceback_counts: prelim={} tree_skipped={} full_traceback={} deleted_evalue_cutoff={} deleted_reevaluation={} deleted_identity_length={} endpoint_pass1_removed={} endpoint_pass2_removed={} final_tree_removed={}",
        timing_count(&timing.traceback_prelim_hsps),
        timing_count(&timing.traceback_tree_precheck_skipped_hsps),
        timing_count(&timing.traceback_full_traceback_hsps),
        timing_count(&timing.traceback_deleted_evalue_cutoff_hsps),
        timing_count(&timing.traceback_deleted_reevaluation_hsps),
        timing_count(&timing.traceback_deleted_identity_length_hsps),
        timing_count(&timing.traceback_removed_endpoint_pass1_hsps),
        timing_count(&timing.traceback_removed_endpoint_pass2_hsps),
        timing_count(&timing.traceback_removed_final_tree_hsps)
    );
    let (edit_count, edit_max, edit_median) =
        timing_length_stats(&timing.traceback_edit_script_lengths);
    let (aln_count, aln_max, aln_median) = timing_length_stats(&timing.traceback_alignment_lengths);
    timing_line!(
        "[TIMING] traceback_lengths: edit_script_samples={} edit_script_max={} edit_script_median={:.1} alignment_samples={} alignment_max={} alignment_median={:.1}",
        edit_count, edit_max, edit_median, aln_count, aln_max, aln_median
    );
    timing_line!(
        "[TIMING] format_output: {:.3}s (calls={})",
        format_s,
        format_n
    );
    if let Some(t_search_start) = t_search_start {
        let search_total = t_search_start.elapsed();
        // Diagnostic (LOSAT_TIMING): the same interval feeds the common line (utils/xtiming.rs):
        // one query batch's search, from after its lookup table to its end (the batches' sum
        // over all threads). No NCBI counterpart.
        crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::Search, search_total);
        timing_line!("[TIMING] search_total: {:.3}s", search_total.as_secs_f64());
    }
    if let Some(t_total) = t_total {
        timing_line!("[TIMING] total: {:.3}s", t_total.elapsed().as_secs_f64());
    }
    text
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:121-161
/// A subject's final HSPs (`SubjectPasses::search_subject`'s traceback).
#[derive(Debug)]
enum SubjectHits {
    /// The HSPs of every query in one list, in any order.
    Mixed(Vec<BlastnHsp>),
    /// No NCBI counterpart (LOSAT_X_QGROUPPAR's traceback): the lists of its queries, one query
    /// each, in query order, each sorted by score (`sort_hsps_by_score`). These are the lists that
    /// `update_hitlists_with_subject_hits` makes of the subject's list sorted by score: that sort
    /// is stable and those lists are sorted by the same comparison, so splitting the joined, sorted
    /// list by query gives each query's list back unchanged. No HSP is moved to make them.
    ByQuery(Vec<Vec<BlastnHsp>>),
}

impl SubjectHits {
    fn len(&self) -> usize {
        match self {
            SubjectHits::Mixed(hits) => hits.len(),
            SubjectHits::ByQuery(lists) => lists.iter().map(Vec::len).sum(),
        }
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// No NCBI counterpart (memory placement only): each list keeps only its length. A list
    /// leaves the traceback with its preliminary count as its capacity, and the purges, the final
    /// tree and the e-value reap drop HSPs from it; a list kept past its subject holds its count.
    /// A shrinking `realloc` (glibc shrinks the block in place; another allocator may copy the
    /// kept HSPs).
    fn fit(&mut self) {
        match self {
            SubjectHits::Mixed(hits) => hits.shrink_to_fit(),
            SubjectHits::ByQuery(lists) => lists.iter_mut().for_each(Vec::shrink_to_fit),
        }
    }
}

fn update_hitlists_with_subject_hits(
    hit_lists: &mut [Option<BlastnHitList>],
    hits: SubjectHits,
    prelim_hitlist_size: usize,
) {
    let hits = match hits {
        SubjectHits::Mixed(hits) => hits,
        SubjectHits::ByQuery(lists) => {
            for mut hsps in lists {
                let Some(first) = hsps.first() else {
                    continue;
                };
                let (q_idx, oid) = (first.q_idx as usize, first.s_idx);
                debug_assert!(hsps.iter().all(|hsp| hsp.q_idx as usize == q_idx));
                if q_idx >= hit_lists.len() {
                    continue;
                }
                // The kept list holds its exact count (the final tree and the e-value reap drop
                // HSPs from the preliminary count the list was given).
                hsps.shrink_to_fit();
                sort_hsps_by_score(&mut hsps);
                let hsp_list = BlastnHspList {
                    oid,
                    query_index: q_idx as u32,
                    hsps,
                    best_evalue: i32::MAX as f64,
                };
                hit_lists[q_idx]
                    .get_or_insert_with(|| BlastnHitList::new(prelim_hitlist_size))
                    .update(hsp_list);
            }
            return;
        }
    };
    if hits.is_empty() {
        return;
    }

    // Each query's HSPs get a vector of their exact count (LOSAT's memory: a vector grown by
    // `push` holds four HSPs from the first one, and most subjects give a query one HSP).
    let mut counts = vec![0usize; hit_lists.len()];
    for hit in &hits {
        if let Some(count) = counts.get_mut(hit.q_idx as usize) {
            *count += 1;
        }
    }
    let mut per_query: Vec<Vec<BlastnHsp>> = counts.into_iter().map(Vec::with_capacity).collect();
    for hit in hits {
        let q_idx = hit.q_idx as usize;
        if q_idx < per_query.len() {
            per_query[q_idx].push(hit);
        }
    }

    for (q_idx, mut hsps) in per_query.into_iter().enumerate() {
        if hsps.is_empty() {
            continue;
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:555
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1374-1382
        sort_hsps_by_score(&mut hsps);

        let oid = hsps[0].s_idx;
        let hsp_list = BlastnHspList {
            oid,
            query_index: q_idx as u32,
            hsps,
            best_evalue: i32::MAX as f64,
        };

        let hit_list =
            hit_lists[q_idx].get_or_insert_with(|| BlastnHitList::new(prelim_hitlist_size));
        hit_list.update(hsp_list);
    }
}

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
    /// The records whose letters the rows read (`queries`, and with `-query_loc` their input
    /// records, whole).
    query_records: crate::algorithm::blastn::pairwise::QueryRecords<'a>,
    /// The texts of the subjects that each query's report reads one by one (its batch's
    /// prefetch off), by query, which the first format's writer places where it reads them
    /// (`QueryWarnings::subject_read`).
    query_reads: &'a [crate::report::query_warnings::QueryReads],
    /// The subject records and their letters (theirs, or their packed subjects': with
    /// LOSAT_X_BNSTORE a subject set's, with `-db` the reported sequences' stored data).
    subjects: SubjectLetters<'a>,
    /// Where the searched letters lie in the records (`-query_loc`, `-subject_loc`).
    ranges: &'a ReportRanges<'a>,
    /// Per query, plus strand: DUST and, with `-lcase_masking`, the input lowercase.
    query_masks: &'a [Vec<MaskedInterval>],
    /// The DUST-only masks of every query (outfmt 11's `LCaseMask`), or with
    /// `-window_masker_db` the WindowMasker masks; empty without outfmt 11.
    archive_dust: &'a [Vec<MaskedInterval>],
    /// With `-filtering_db`, the repeat hits of every query that its batch's set-up found
    /// (outfmt 11's `LCaseMask`, `BatchResults::repeat_hits`); empty without outfmt 11 or without
    /// the option.
    archive_repeats: &'a [Vec<(usize, usize)>],
    lcase_masking: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`).
    query_eff_searchsp: &'a [i64],
    /// The first valid context of a query, as its offset in the query's two contexts: 1 for
    /// `-strand minus`, else 0 (`context 2 * q + first_context`).
    first_context: usize,
    megablast: bool,
    /// The search used a MegaBLAST database index (`BlastnIntro::indexed_megablast`).
    indexed_megablast: bool,
    /// Matrix-only scoring (reward 0 and penalty 0, `TaskConfig::matrix_only`): the outfmt 0
    /// report shows raw scores only (`BlastnIntro`).
    matrix_only: bool,
    reward: i32,
    penalty: i32,
    gap_open: i32,
    gap_extend: i32,
    /// The task is megablast or blastn (`BlastnPairwiseReport::zero_gap_extension_formula`).
    zero_gap_extension_formula: bool,
    /// The two-hit window of the task, printed in the epilog when not 0.
    window_size: usize,
    /// An ungapped (`-ungapped`) search: one Karlin block, no gap costs in the reports, and the
    /// hit list of the end of the traceback is the preliminary one (`finish_query_hit_list`).
    ungapped: bool,
    /// Per query context (the plus strand of query `q` is context `2 * q`); `None` for
    /// the contexts of an invalid query (not searched by NCBI) and of the strand that
    /// `-strand` does not search.
    context_karlin: &'a [Option<ContextKarlin>],
    /// The numbers of descriptions and alignments of the pairwise report
    /// (`BlastnArgs::report_counts`).
    num_descriptions: usize,
    num_alignments: usize,
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:251-253
    /// `-line_length`, `-sorthits` and `-sorthsps` of the pairwise report
    /// (`BlastnArgs::formatting`).
    formatting: crate::blastinput::report_counts::FormattingArgs,
    /// `-culling_limit` (0: no culling), read when the hit lists are finished.
    culling_limit: i32,
    /// `-best_hit_overhang` and `-best_hit_score_edge` (`None`: neither given), read when the
    /// hit lists are finished.
    best_hit: Option<crate::core::hspfilter_besthit::BestHitParams>,
    db_num_seqs: i32,
    /// NCBI's `Int8` (`SubjectMetadata::db_len_total`).
    db_len_total: i64,
    /// The `Database:` lines of a database search's outfmt 0 report.
    database: Option<&'a PairwiseDatabase>,
    /// The tabular subject ids of a database search's reported sequences, by rank
    /// (`database::ReportedSubjects::tabular_subject_ids`); `None` for a subject set.
    database_subject_ids: Option<&'a [crate::algorithm::common::tabular::SubjectId]>,
    /// The ordinal ID of each reported sequence of a database search, by rank
    /// (`database::ReportedSubjects::ordinals`); `None` for a subject set.
    database_ordinals: Option<&'a [u32]>,
    /// The opened database of a database search, which the XML2 and JSON formats read the taxids
    /// of the reported sequences from and whose OIDs name the sequences without parsed Seq-ids
    /// (`oid_map`); `None` for a subject set.
    database_reader: Option<&'a database::BlastnDatabase>,
    /// The algorithm of a database search with `-db_soft_mask` or `-db_hard_mask`
    /// (`ResolvedMask::id`), whose masks outfmt 0 shows in lowercase (`database::BlastnDatabase::display_masks`); `None` otherwise.
    database_mask: Option<i32>,
    /// The reported sequences with parsed Seq-ids of a database search, by rank
    /// (`database::ReportedSubjects::parsed`; empty without them): their ids and titles come from
    /// their definition lines. Of a subject set, the records whose ids `-parse_deflines` read, by
    /// subject index (`db_subject::parsed_records`).
    database_parsed: &'a [Option<Arc<crate::report::db_subject::ParsedSubject>>],
    /// What the XML2 and JSON reports of a database search do at each subject's `Hit`, by rank
    /// (`database::ReportedSubjects::loads`; every subject loaded for a subject set).
    database_loads: &'a crate::report::db_subject::HitLoads,
    /// The sources from which NCBI makes the empty stored titles of the other reported sequences
    /// of a database search, by rank (`database::ReportedSubjects::sources`; empty without them
    /// and for a subject set).
    database_sources: &'a [crate::report::defline::SourceTaxname],
    /// Whether the hit lists were finished after each query batch's traceback (a database
    /// search, `run_in_pool`), not by `post_process_hits_and_write`.
    finished_hit_lists: bool,
    /// Whether NCBI's epilog (`PrintEpilog`: the outfmt 0 footer, outfmt 7's last line)
    /// follows the reports; an error in a later batch skips it.
    epilog: bool,
    /// `-show_gis` (the labels of the query-anchored sections, outfmt 1 to 4).
    show_gis: bool,
    /// The `-html` links of the outfmt 0 to 4 reports (`blastinput::app::html_links`); None
    /// without `-html` or any of those formats.
    html: Option<Arc<crate::report::html_ids::HtmlLinks>>,
    /// `-html` (`BlastnArgs::html`, NCBI's `m_IsHTML`): the HTML form of the outfmt 0 to 4 and
    /// 18 reports (the organism report takes no links).
    is_html: bool,
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
    /// The database of `-filtering_db` (the filter string's `R -d <name>;`).
    repeat_filter: Option<&'a str>,
    /// The WindowMasker option (the filter string's `W -d <name>;` or `W -t <taxid>;`).
    window_masker: Option<crate::report::blast_xml::filter::WindowMasker<'a>>,
    /// `-soft_masking` (the filter string's `m;`).
    soft_masking: bool,
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
    /// Whether `taxdb.bti` is on NCBI's search path, for the custom tabular reports
    /// (`ShownSubjects::tabular_taxdb_found`).
    tabular_taxdb_found: bool,
    /// The taxonomy database of a `-subject` search's organism report (`BlastnArgs::subject_taxdb`).
    subject_taxdb: Option<&'a crate::algorithm::common::database::SubjectTaxdb>,
}

impl BlastnReportInputs<'_> {
    /// The number of sequences the reports print: the database's own (with a taxid list, the
    /// whole database's: `SeqDb::printed_num_seqs`), or the subject set's; NCBI's `Int8`
    /// (`PairwiseDatabase::num_sequences`).
    fn shown_num_seqs(&self) -> i64 {
        self.database
            .map_or(i64::from(self.db_num_seqs), |database| {
                database.num_sequences
            })
    }

    /// The total length the reports print (as `shown_num_seqs`), NCBI's `Int8`
    /// (`PairwiseDatabase::total_letters`).
    fn shown_total_letters(&self) -> i64 {
        self.database
            .map_or(self.db_len_total, |database| database.total_letters)
    }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-104
    // The first valid context of a query is its minus context with `-strand minus`.
    let first_context = report.first_context;
    let queries = report
        .queries
        .iter()
        .enumerate()
        .map(|(q_idx, query)| BlastnPairwiseQuery {
            query_name: query_titles[q_idx].to_vec(),
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:742-744
            // The length of the query's record, also where `-query_loc` searches part of it.
            query_length: report
                .ranges
                .queries
                .placements
                .length(q_idx, query.seq().len()),
            karlin: context_karlin[2 * q_idx + first_context]
                .map(|blocks| (blocks.ungapped, blocks.gapped)),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-104
            // An invalid query has no valid context, and its search space stays 0.
            effective_search_space: if context_karlin[2 * q_idx + first_context].is_some() {
                report.query_eff_searchsp[2 * q_idx + first_context]
            } else {
                0
            },
        })
        .collect();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2910-2928
    // The defaults are 500 descriptions and 250 alignments (format_flags.cpp:219,221).
    let (num_descriptions, num_alignments) = (report.num_descriptions, report.num_alignments);
    let pairwise_report = BlastnPairwiseReport {
        version: crate::wording::REPORT_VERSION.to_string(),
        intro: BlastnIntro {
            megablast: report.megablast,
            indexed_megablast: report.indexed_megablast,
            matrix_only: report.matrix_only,
        },
        database_name: subject_title.as_bytes().to_vec(),
        database_num_sequences: report.shown_num_seqs(),
        database_total_letters: report.shown_total_letters(),
        database: report.database.cloned(),
        reward: report.reward,
        penalty: report.penalty,
        gap_open: report.gap_open,
        gap_extend: report.gap_extend,
        zero_gap_extension_formula: report.zero_gap_extension_formula,
        window_size: report.window_size,
        gapped: !report.ungapped,
        num_descriptions,
        num_alignments,
        unsearched,
        epilog: report.epilog,
        prolog: report.prolog,
        parsed_subjects: report.database_parsed.to_vec(),
        subject_sources: report.database_sources.to_vec(),
    };
    Ok((queries, pairwise_report))
}

/// The queries of the call's batches as the archive takes them (`archive::for_each_batch`).
fn call_batches<'a>(
    report: &'a BlastnReportInputs<'a>,
    query_lengths: &'a [usize],
    unsearched: &'a [bool],
) -> archive::CallBatches<'a> {
    archive::CallBatches {
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
        query_masks: report.query_masks,
        dust_masks: report.archive_dust,
        repeat_hits: report.archive_repeats,
        unsearched,
        lcase_masking: report.lcase_masking,
        window_masker: report.window_masker.is_some(),
        both_strands: report.both_strands,
        minus_strand: report.minus_strand,
        spans: report.batch_spans,
        gapped: !report.ungapped,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
/// The batch's context 0 as the `EffectiveSearchSpace` option of the archive and the search
/// strategy reads it (`archive::RequestLength::space`), for a batch whose first query is `first`
/// (an index of `query_eff_searchsp`'s and `context_karlin`'s queries, two contexts each) and of
/// `query_length` letters: the search space of its context 0, the option
/// `SplitQuery_SetEffectiveSearchSpace` sets, and that context's query length. Context 0 is the
/// plus strand (`-strand minus` leaves it invalid: no Karlin block, length 0). The archive's
/// (`archive_first_context`) and the search strategy's (`run_in_pool`,
/// `StrategyExport::note_searched`).
fn context0_searched(
    query_eff_searchsp: &[i64],
    context_karlin: &[Option<ContextKarlin>],
    query_length: usize,
    first: usize,
) -> crate::report::archive::SearchedContext0 {
    crate::report::archive::SearchedContext0 {
        space: query_eff_searchsp[2 * first],
        query_length: if context_karlin[2 * first].is_some() {
            query_length as i64
        } else {
            0
        },
    }
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
    let context0 = context0_searched(
        report.query_eff_searchsp,
        report.context_karlin,
        query_lengths[first],
        first,
    );
    let Some(blocks) = report.context_karlin[2 * first + report.first_context] else {
        return Ok(archive::FirstContext {
            karlin: None,
            search_space: 0,
            length_adjustment: 0,
            context0_space: context0.space,
            context0_length: context0.query_length,
        });
    };
    let search_space = report.query_eff_searchsp[2 * first + report.first_context];
    // `stats::setup::length_adjustment`: the setup's (`calc_eff_lengths`) for this context.
    let length_adjustment = i64::from(crate::stats::setup::length_adjustment(
        query_lengths[first] as i64,
        report.eff_db_length,
        report.db_num_seqs,
        &blocks.gapped,
    ));
    Ok(archive::FirstContext {
        karlin: Some(blocks),
        search_space,
        length_adjustment,
        context0_space: context0.space,
        context0_length: context0.query_length,
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
            // `stats::setup::length_adjustment`: the setup's (`calc_eff_lengths`) for this
            // context.
            let length_adjustment = i64::from(crate::stats::setup::length_adjustment(
                query_lengths[q_idx] as i64,
                report.eff_db_length,
                report.db_num_seqs,
                &blocks.gapped,
            ));
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
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1520-1534
/// NCBI reference (598d8ae6): c++/src/objmgr/util/seq_align_util.cpp:71-78
/// A subject interval has the strand both, so its mapping is the shift by its start.
///
/// An HSP that NCBI keeps outside its query context with a Seq-align query start below 0
/// (`BlastnHsp::seqalign_query_start`) reaches the reports only where the shift moves its
/// Dense-diag's start to 0 or past it (`BlastnHsp::record_query_start`,
/// `left_of_query_refusal`): its offsets, which `adjust_blastn_offsets` kept at 0, are those of
/// the moved start (`Blast_HSPGetAdjustedOffsets` of its `Int4` offsets plus the shift).
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
                if hsp.seqalign_query_start() < 0 {
                    let start = hsp.record_query_start(q_shift, true);
                    let length =
                        i64::from(hsp.internal_q_end_0) - i64::from(hsp.internal_q_offset_0);
                    hsp.q_start = (start + 1) as usize;
                    hsp.q_end = (start + length) as usize;
                } else {
                    hsp.q_start += q_shift;
                    hsp.q_end += q_shift;
                }
                hsp.s_start += s_shift;
                hsp.s_end += s_shift;
            }
        }
    }
}

/// The subject ID of the tabular formats (`sseqid`, `sacc`, `saccver`) of every subject, by
/// subject index (`report::defline::tabular_subject_id`: the title's first word, or the
/// local ID, and for `lcl|Subject_...` the first word of `GenerateDefline`); `None` where
/// `GenerateDefline` throws (`common::tabular::SubjectId`) (a `Subject_` first word and a title in no encoding that
/// `CUtf8::GuessEncoding` recognises), whose rows NCBI writes with its catch's values
/// (`report::tabular::SubjectLookup::Failed`) or stops at (`check_failed_lookups`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:828-841
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

/// The check of the rows of `hit_lists` whose subjects' lookups fail where NCBI has no result
/// (`report::tabular::check_failed_lookups`), for each tabular format of `output_formats` (the
/// formats of `outputs`, whose refusal each is: `FormatRefusal::check`), before any format is
/// written; `subject_ids` are the formats' subject ids (`tabular_subject_ids`, or a database's).
/// Out of line and called only for a tabular format whose subjects' lookups may fail
/// (`common::tabular::lookups_may_fail`), so that the other reports do not change for it (no NCBI
/// counterpart: code placement only).
#[inline(never)]
fn check_tabular_lookups(
    outputs: &mut [crate::api::local_blast::FormatOutput<'_>],
    output_formats: &[BlastnOutputFormat],
    hit_lists: &[Option<BlastnHitList>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    report: &BlastnReportInputs<'_>,
) -> Result<()> {
    let database = report.database_ordinals.zip(report.database_reader);
    for (output, format) in outputs
        .iter_mut()
        .zip(output_formats)
        .filter(|(_, format)| format.is_tabular())
    {
        let spec = format.tabular_spec();
        if !crate::algorithm::common::tabular::checks_failed_lookups(
            spec,
            subject_ids,
            report.database_parsed,
            report.database_loads,
        ) {
            continue;
        }
        let checked = crate::report::tabular::check_failed_lookups(
            spec,
            hit_lists.iter().flatten().map(|hit_list| {
                hit_list
                    .hsplist_array
                    .iter()
                    .filter_map(|hsp_list| hsp_list.hsps.first())
                    .map(|hsp| hsp.s_idx as usize)
            }),
            |s_idx| {
                crate::algorithm::common::tabular::subject_lookup_failure(
                    spec,
                    subject_ids,
                    report.database_parsed,
                    report.database_loads,
                    database.is_some(),
                    false,
                    s_idx,
                )
            },
            |s_idx| match database {
                Some((ordinals, reader)) => crate::report::tabular::failed_database_subject(
                    reader.oid_map().oid(ordinals[s_idx]),
                    report.database_loads.thrown(s_idx),
                    report.database_parsed.get(s_idx).and_then(Option::as_deref),
                ),
                None => crate::report::tabular::failed_record_title(
                    s_idx,
                    report.subjects.records()[s_idx].shown_id(),
                ),
            },
            "BLASTN",
        );
        output.refusal.check(checked)?;
    }
    Ok(())
}

/// The end of the traceback for the final hit list of one query (`hit_list`, of a query of
/// `query_len` letters, the length of each of its contexts): the HSP lists cut to
/// `max_hsps_per_subject`, their HSPs below `query_cov_hsp_perc` percent of the query
/// removed (`-qcov_hsp_perc`; a list left empty is removed) and their subject best hits
/// (`-subject_besthit`), the lists sorted by e-value and cut to `hitlist_size`, and the HSPs
/// of each list sorted for the Seq-aligns (`post_process_hits_and_write`; a database search
/// does it after each query batch's traceback, `run_in_pool`).
fn finish_query_hit_list(
    hit_list: &mut BlastnHitList,
    query_len: usize,
    hitlist_size: usize,
    max_hsps_per_subject: usize,
    query_cov_hsp_perc: f64,
    subject_besthit: bool,
    culling_limit: i32,
    best_hit: Option<crate::core::hspfilter_besthit::BestHitParams>,
    ungapped: bool,
    sorts_results: bool,
) {
    // A database whose `BlastSeqSrcGetTotLen` is not positive keeps the traceback's order, which
    // the post-traceback steps below then read: the traceback saved the final lists in the order
    // of their ordinal ids into a hit list of `hitlist_size` (`update_hitlists_with_subject_hits`,
    // `Blast_HSPResultsInsertHSPList`), the order of NCBI's one thread.
    // NCBI's culling pipe is the first thing the traceback does with the consolidated results
    // (`-culling_limit`; the mask level of 101, the default, is no step).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1721-1723,1763-1766
    if culling_limit > 0 {
        culling::cull_traceback(
            hit_list,
            query_len,
            culling_limit,
            get_prelim_hitlist_size(hitlist_size, false, !ungapped),
        );
    }
    // NCBI's best hit pipe is there in its place (`-best_hit_overhang`, `-best_hit_score_edge`;
    // the two exclude `-culling_limit`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:377-385
    if let Some(params) = best_hit {
        besthit::best_hit_traceback(
            hit_list,
            query_len,
            params,
            get_prelim_hitlist_size(hitlist_size, false, !ungapped),
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:836-867
    for hsp_list in &mut hit_list.hsplist_array {
        if max_hsps_per_subject > 0 {
            trim_by_max_hsps(hsp_list, max_hsps_per_subject);
        }
        // Both contexts of a query have its length, and the offsets are those of the HSP's
        // own context.
        if query_cov_hsp_perc > 0.0
            && reap_by_query_coverage(&mut hsp_list.hsps, query_cov_hsp_perc, |hsp| {
                (
                    i64::from(hsp.internal_q_offset_0),
                    i64::from(hsp.internal_q_end_0),
                    query_len as i64,
                )
            })
        {
            update_best_evalue(hsp_list);
        }
        // NCBI calls Blast_HSPListSubjectBestHit on a list the coverage left empty after
        // freeing it (blast_traceback.c:853-861); LOSAT reads it as the empty list it was.
        if subject_besthit {
            subject_best_hit(&mut hsp_list.hsps, query_len);
        }
    }
    if query_cov_hsp_perc > 0.0 {
        hit_list
            .hsplist_array
            .retain(|hsp_list| !hsp_list.hsps.is_empty());
        hit_list.hsplist_count = hit_list.hsplist_array.len();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3383-3400
    if sorts_results {
        hit_list.sort_by_evalue();
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-892
    hit_list.prune_by_size(hitlist_size);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1571-1577
    // Every format shows the HSPs of a subject in this order. It is the score order,
    // except where the two strands of a query have different gapped blocks.
    for hsp_list in &mut hit_list.hsplist_array {
        sort_hsplist_by_evalue(hsp_list);
    }
}

/// The HSPs of score 0 dropped from the finished hit list `hit_list`, and the HSP lists
/// left without HSPs (`post_process_hits_and_write`).
fn drop_hsps_of_score_zero(hit_list: &mut BlastnHitList) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:672-674
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1485-1490
    // An HSP of score 0 (the traceback reads the ambiguity codes that the preliminary
    // search read as random bases) has no Seq-align, so no format shows it, and a
    // subject without other HSPs has no alignments.
    for hsp_list in &mut hit_list.hsplist_array {
        hsp_list.hsps.retain(|hsp| hsp.raw_score != 0);
    }
    hit_list
        .hsplist_array
        .retain(|hsp_list| !hsp_list.hsps.is_empty());
    hit_list.hsplist_count = hit_list.hsplist_array.len();
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
    query_warnings: &[Posts],
    mut sam: Vec<Option<SamTarget<'_>>>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:183-187
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:836-867
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
                report.best_hit,
                report.ungapped,
                // A FASTA subject set (NCBI's `dbscan` source totals its records): a set
                // without letters finds no hit.
                true,
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:877-893
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

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:3077-3106
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
    // Each requested format prints the same final result without searching again.
    // `unsearched` marks the queries of the batches without a valid context, which NCBI
    // does not search (outfmt 0 and 7 show them without results).
    // Rendered hits are needed by outfmt 0 and by the caller's hit records; both get the
    // same final hit list, built once.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
    // The reports print record coordinates: the HSPs move by the start of the query's and
    // the subject's interval (`shift_to_records`); their internal offsets stay those of
    // the searched letters, from which the shown rows are read.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1619-1636
    // An HSP that NCBI keeps outside its query context reaches here only with a Seq-align query
    // start of at least 0 (`left_of_query_refusal` fails or refuses the batch of one below 0), and
    // the reports show it with query coordinates past the query's end (remaining-gaps X1b-p).
    shift_to_records(&mut hit_lists, report.ranges);
    let wants_rows =
        outputs.hits.is_some() || output_formats.iter().any(|format| format.is_pairwise());
    // EXPERIMENT (LOSAT_X_RPTSTREAM): the same masks, made once for every batch.
    let subject_masks = (wants_rows && report.lcase_masking).then(|| {
        // The records as read with `-subject_loc` (never packed), else the searched ones.
        let subject_records = match report.ranges.subjects.records {
            Some(records) => SubjectLetters::new(records, None),
            None => report.subjects,
        };
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1580-1585
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:154-157
    // The stored masks of the reported sequences of a masked database search, by rank: the
    // displayed rows show the letters of those that meet an HSP in lowercase (`pairwise::Rows`).
    let database_masks = match (
        wants_rows,
        report.database_mask,
        report.database_reader,
        report.database_ordinals,
    ) {
        (true, Some(algo_id), Some(reader), Some(ordinals)) => {
            let mut masks = crate::seqdb::MaskBlock::default();
            Some(
                ordinals
                    .iter()
                    .map(|&ordinal| reader.display_masks(ordinal, algo_id, &mut masks))
                    .collect::<Result<Vec<_>>>()?,
            )
        }
        _ => None,
    };
    // The shown masks of every query and subject, made once for every window of queries.
    let rows = wants_rows.then(|| {
        PairwiseRows::new(
            report.query_records,
            report.subjects,
            &DisplayMasks {
                query: report.query_masks,
                query_offsets: &report.ranges.queries.placements,
                subject: subject_masks.as_deref(),
                database: database_masks.as_deref(),
            },
        )
    });
    // The rendered hits of the queries `queries`, under `-html` with their rows' masks in
    // `html_masks` (`pairwise_hits`).
    let rows_of = |queries: std::ops::Range<usize>,
                   html_masks: Option<&mut Vec<crate::report::query_anchored::HspMasks>>|
     -> Result<Vec<PairwiseHit>> {
        let rows = rows
            .as_ref()
            .expect("the rows' inputs are made for outfmt 0");
        let mut hits = pairwise_hits(
            &hit_lists[queries],
            rows,
            report.x_pairwise_pool,
            html_masks,
        )?;
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2479
        // The length of the subject's record, also where `-subject_loc` searches part
        // of it.
        for hit in &mut hits {
            let s_idx = hit.hit.s_idx as usize;
            hit.subject_length = hit
                .subject_length
                .map(|searched| report.ranges.subjects.placements.length(s_idx, searched));
        }
        Ok(hits)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1585
    // NCBI makes and prints the alignments of one query's result set at a time; outfmt 0 makes
    // the rows of a window of queries, writes them and drops them before the next window's
    // (`pairwise::pairwise_windows`, line DBID, ID4f). The rows' errors are LOSAT's internal
    // invariants (`pairwise_hits`), so making them while the report is written changes no
    // run's output. The caller's hit records and a format observer (whose HSP numbers count
    // over the whole report) read every query's rows: those are made once, before any format
    // is written.
    // Under `-html` the rows' masks of those hits (`HtmlRowMasks`, the spans of outfmt 0).
    let mut html_masks = report.html.is_some().then(Vec::new);
    let pairwise_hits = if outputs.hits.is_some()
        || (outputs.observer.is_some() && output_formats.iter().any(|format| format.is_pairwise()))
    {
        Some(Arc::new(rows_of(0..hit_lists.len(), html_masks.as_mut())?))
    } else {
        None
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // The caller receives the same final result that every formatter prints.
    if let (Some(hits_sink), Some(hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref()) {
        hits_sink(Arc::clone(hits));
    }
    // Built before any format is written, so that a report that cannot be made fails the
    // run without output.
    let pairwise_data = if output_formats.iter().any(|format| format.is_pairwise()) {
        Some(blastn_pairwise_report(
            report,
            unsearched.clone(),
            query_titles,
            subject_title,
        )?)
    } else {
        None
    };
    // EXPERIMENT (LOSAT_X_RPTSTREAM): a subject set's ids made once for every batch where LOSAT
    // refuses none of them (`x_rpt_stream::SubjectCache`), so they are the ids made here.
    let tabular = output_formats.iter().any(|format| format.is_tabular());
    let subject_ids: std::borrow::Cow<'_, [crate::algorithm::common::tabular::SubjectId]> =
        if tabular {
            match (
                report.database_subject_ids,
                report.subject_cache.and_then(|cache| cache.ids.as_deref()),
            ) {
                (Some(ids), _) => std::borrow::Cow::Borrowed(ids),
                (None, Some(ids)) => std::borrow::Cow::Borrowed(ids),
                (None, None) => {
                    std::borrow::Cow::Owned(tabular_subject_ids(report.subjects.records()))
                }
            }
        } else {
            std::borrow::Cow::Owned(Vec::new())
        };
    // The rows of subjects whose lookups fail where NCBI has no result, before any format is
    // written.
    if tabular
        && crate::algorithm::common::tabular::lookups_may_fail(&subject_ids, report.database_parsed)
    {
        check_tabular_lookups(
            &mut outputs.formats,
            output_formats,
            &hit_lists,
            &subject_ids,
            report,
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:497-498
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2273
    // The outfmt 0 descriptions and headings are made from the subjects' title bytes
    // (`report::defline::generate_defline`), by subject index.
    // EXPERIMENT (LOSAT_X_RPTSTREAM): the same titles, made once for every batch.
    let subject_titles: std::borrow::Cow<'_, [Arc<[u8]>]> =
        if output_formats.iter().any(|format| format.is_pairwise()) {
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
    // The query-anchored formats (outfmt 1 to 4) refuse the call before any format of it is written.
    report_write::refuse_anchored_reports(
        output_formats,
        report,
        &hit_lists,
        pairwise_hits.as_deref().map(Vec::as_slice),
        &subject_titles,
        // The check reads the hits' segments only, not the rows' masks of `-html`'s spans.
        &|queries| rows_of(queries, None),
    )?;
    report_write::refuse_described_records(
        output_formats,
        &mut outputs.formats,
        report,
        &hit_lists,
        pairwise_hits.as_deref().map(Vec::as_slice),
        &|queries| rows_of(queries, None),
    )?;
    // The queries' result warnings that a structured first format does not write are each later
    // format's that writes them, to the caller's receiver of that format's
    // (`ReportOutputs::later_warnings`).
    // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    let keeps_later = outputs.keeps_later_query_warnings();
    let observer = &mut outputs.observer;
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
    // The statistics of each query for the observer, the values of the XML reports: `q_idx` is
    // that of the hit lists (a caller with an observer gets one report of all queries,
    // `x_rpt_stream`).
    if let Some(observer) = observer.as_deref_mut() {
        let stats = blastn_xml_stats(report, query_lengths)?;
        for (q_idx, stats) in stats.into_iter().enumerate() {
            observer.query_stats(q_idx, stats);
        }
    }
    // The warnings are written once, with the first format (`QueryWarnings`).
    let mut warning_lines = Some(query_warnings);
    // Every structured format of the call is XML2 or JSON: none writes the report that a
    // `-subject` document holds back (`blast_xml::write_blastn_xml`).
    let xml_only_xml2 = output_formats
        .iter()
        .filter(|format| format.is_structured_xml())
        .all(|format| format.xml2_kind().is_some());
    for (format_index, (format, &output_format)) in
        outputs.formats.iter_mut().zip(output_formats).enumerate()
    {
        // A format that failed alone writes nothing more; the next takes the warnings.
        if format.refusal.refused() {
            continue;
        }
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
        // So do those of the formatters of a tabular field list that shows a taxonomy name
        // (`x_CheckTaxDB`, `algorithm::blastn::tabular`), and the Seq-annots' posts of the bytes
        // that NCBI's stream replaces (`seqalign::write_query_to`).
        // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
        let later = format.later_warnings(keeps_later, &mut outputs.later_warnings, format_index);
        let mut format_warnings = (warning_lines.is_some()
            || later.is_some()
            || crate::report::query_warnings::report_posts(
                format.report_format(),
                report.html.is_some(),
            )
            || output_format
                .tabular_spec()
                .is_some_and(|spec| spec.looks_up_every_taxid()))
        .then(|| {
            let first = warning_lines.take();
            let warnings = crate::report::query_warnings::QueryWarnings::of_posts(
                first.unwrap_or(&[]),
                report.archive,
                &mut *outputs.diagnostics,
            )
            .with_later(query_warnings, later)
            .for_report(format_index)
            .with_html_posts(
                outputs
                    .html_posts
                    .as_deref_mut()
                    .map(|posts| posts as &mut dyn std::io::Write),
            );
            // The first format's writer also writes the texts of the subjects that the reports
            // read one by one (`BlastnReportInputs::query_reads`).
            match first {
                Some(_) => warnings.with_subject_reads(report.query_reads),
                None => warnings,
            }
        });
        // The caller's receiver of the side files of outfmt 13 and 14, if it gives one.
        let side_sink = format.sink.side_files().cloned();
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
        // The format's report, whose error goes to its refusal (below).
        let written = (|| -> anyhow::Result<()> {
            let mut writer = crate::cli::NcbiOstream::sink(&mut format.sink)?;
            if let BlastnOutputFormat::Pairwise(view) = output_format {
                let (queries, pairwise_report) = pairwise_data
                    .as_ref()
                    .expect("the pairwise report data is built for outfmt 0");
                // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:251-253
                let config = PairwiseConfig {
                    program: "blastn".to_string(),
                    show_frame: false,
                    line_length: report.formatting.line_length(),
                    // `-sorthits` for every view, `-sorthsps` for outfmt 0 only (NCBI's warning for the
                    // others, `PairwiseView::report_sort`).
                    sort: view.report_sort(report.formatting.report_sort()),
                    view,
                    html: report.html.as_ref().map(|links| {
                        crate::report::html_ids::HtmlReport::new(
                            links,
                            report
                                .database_ordinals
                                .zip(report.database_reader)
                                .into_iter()
                                .flat_map(|(ordinals, reader)| {
                                    ordinals.iter().map(move |&o| reader.oid_map().oid(o))
                                }),
                        )
                    }),
                    ..PairwiseConfig::default()
                };
                // The sequences of the query-anchored sections of outfmt 1 to 4 (line RG T3b).
                let anchored =
                    (view != crate::report::pairwise::PairwiseView::Pairwise).then(|| {
                        crate::algorithm::blastn::anchored::AnchoredSequences {
                            queries: report.queries,
                            query_records: report.query_records,
                            subjects: report.subjects,
                            subject_placements: report.ranges.subjects.placements,
                            parsed: report.database_parsed,
                            listed_ordinals: report
                                .database_ordinals
                                .zip(report.database_reader.map(|reader| reader.oid_map())),
                            show_gis: report.show_gis,
                            rows: rows
                                .as_ref()
                                .expect("the rows' inputs are made for outfmt 0 to 4"),
                            hit_lists: &hit_lists,
                        }
                    });
                let mut write = |hits: &[PairwiseHit],
                                 masks: Option<&[crate::report::query_anchored::HspMasks]>,
                                 range: std::ops::Range<usize>| {
                    write_blastn_pairwise_report(
                        hits,
                        &mut writer,
                        &config,
                        queries,
                        range,
                        &subject_titles,
                        pairwise_report,
                        probe.as_mut(),
                        format_warnings.as_mut(),
                        report.x_pairwise_pool,
                        anchored.as_ref().map(|sequences| {
                            sequences as &dyn crate::report::query_anchored::ReportSequences
                        }),
                        masks,
                        report.database_loads,
                    )
                    // A subject whose load throws may end the run there
                    // (`x_pairwise_par::alignment_end`).
                    .map_err(crate::cli::NativeError::of_writer)
                };
                match pairwise_hits.as_deref() {
                    Some(hits) => write(hits, html_masks.as_deref(), 0..queries.len())?,
                    None => {
                        // The windows of the parallel writer (`x_pairwise_par`).
                        let windows = pairwise_windows(
                            &hit_lists,
                            queries.len(),
                            report.x_pairwise_pool,
                            config.line_length,
                        );
                        if windows.is_empty() {
                            // No query: the prolog and the epilog alone.
                            write(&[], None, 0..0)?;
                        }
                        for range in windows {
                            let lists =
                                range.start.min(hit_lists.len())..range.end.min(hit_lists.len());
                            let mut masks = report.html.is_some().then(Vec::new);
                            let hits = rows_of(lists, masks.as_mut())?;
                            write(&hits, masks.as_deref(), range)?;
                        }
                    }
                }
            } else if let BlastnOutputFormat::Seqalign(serial) = output_format {
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:690-698
                // The Seq-annot of each query with alignments (`seqalign`).
                let subjects = report_write::subject_names_of(report);
                write_blastn_seqalign(
                    serial,
                    !report.ungapped,
                    &hit_lists,
                    &QueryNames {
                        records: report.queries,
                        searched_lengths: query_lengths,
                        placements: &report.ranges.queries.placements,
                    },
                    subjects,
                    report.ranges.subjects.placements,
                    &mut writer,
                    format_warnings
                        .as_mut()
                        .expect("a Seq-annot writes its posts to the diagnostics"),
                )?;
            } else if output_format == BlastnOutputFormat::Archive {
                // The BLAST archive (outfmt 11): one document per query batch of the call.
                let handle = report
                    .archive
                    .expect("a BLAST archive request has its handle");
                let subjects = report_write::subject_names_of(report);
                let call = archive::ArchiveCall {
                    handle,
                    batches: call_batches(report, query_lengths, &unsearched),
                    subjects,
                    subject_placements: report.ranges.subjects.placements,
                    hit_lists: &hit_lists,
                };
                archive::write_blastn_archives(
                    &call,
                    format_warnings
                        .as_mut()
                        .expect("the BLAST archive writes its posts to the diagnostics"),
                    &|first| archive_first_context(report, query_lengths, first),
                    &mut writer,
                )?;
            } else if let BlastnOutputFormat::Sam(_) = output_format {
                report_write::write_sam(
                    &mut sam[format_index],
                    report,
                    &hit_lists,
                    query_lengths,
                    format_warnings.as_mut(),
                    &mut writer,
                )?;
            } else if output_format.is_structured_xml() {
                report_write::write_structured_xml(
                    report_write::XmlDocument {
                        kind: output_format.xml2_kind(),
                        only_xml2: xml_only_xml2,
                        side_sink: side_sink.as_ref(),
                        prefetch: database::prefetch_format(output_formats),
                    },
                    report,
                    &hit_lists,
                    query_lengths,
                    &unsearched,
                    subject_title,
                    &mut writer,
                    format_warnings.as_mut(),
                )?;
            } else if let Some(spec) = output_format.tabular_spec() {
                // The tabular reports with a field list or a delimiter, and outfmt 10 and 20
                // (`algorithm::blastn::tabular`); never the default 6 and 7 below.
                write_custom_tabular(
                    spec,
                    report,
                    &hit_lists,
                    query_ids,
                    query_titles,
                    query_lengths,
                    subject_title,
                    &subject_ids,
                    &unsearched,
                    &mut writer,
                    format_warnings.as_mut(),
                    probe.as_mut(),
                )?;
            } else {
                crate::algorithm::blastn::x_tab_par::write_output(
                    report.x_tab_pool,
                    &hit_lists,
                    &mut writer,
                    query_ids,
                    &subject_ids,
                    report.database_loads,
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
            std::io::Write::flush(&mut writer)?;
            Ok(())
        })();
        // An error of the format's report fails the call, or this format alone (`FormatRefusal`):
        // the formats after it are written then.
        format.refusal.check(written)?;
    }
    // Without an output format the warnings are written all at once.
    if let Some(handle) = report.archive {
        handle.set_muted(false);
    }
    if let Some(before) = warning_lines {
        for query_warnings in before {
            query_warnings.write(&mut *outputs.diagnostics, report.archive)?;
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:770-782
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

/// The custom tabular report of the queries of `hit_lists` (`algorithm::blastn::tabular`): the
/// subjects' values as the report shows them (a database's stored lengths and, where a field
/// shows them, its taxids), outfmt 20's header line with the run's first report and outfmt 7's
/// last line with the epilog. Kept out of line: only a field list reaches it.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-846
/// `slen` is the length of the subject's Bioseq: the whole record with `-subject_loc`, as the
/// outfmt 0 `Length=` line.
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn write_custom_tabular(
    spec: &crate::report::tabular::TabularSpec,
    report: &BlastnReportInputs<'_>,
    hit_lists: &[Option<BlastnHitList>],
    query_ids: &[Arc<[u8]>],
    query_titles: &[Arc<[u8]>],
    query_lengths: &[usize],
    subject_title: &str,
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    unsearched: &[bool],
    writer: &mut crate::cli::NcbiOstream<'_>,
    warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    probe: Option<&mut FormatProbe<'_>>,
) -> Result<()> {
    use super::super::tabular::{write_blastn_tabular, BlastnTabular};
    use crate::algorithm::common::tabular::TabularSubjects;
    let needs = spec.needs();
    let record_length = |s_idx: usize| {
        report
            .ranges
            .subjects
            .placements
            .length(s_idx, report.subjects.len(s_idx))
    };
    let (lengths, taxonomy) = match (report.database_ordinals, report.database_reader) {
        (Some(ordinals), Some(database)) => (
            if spec.has(crate::report::tabular::TabularField::SubjectLength) {
                database.sam_lengths(ordinals)?
            } else {
                Vec::new()
            },
            if needs.taxonomy {
                database.tabular_taxonomy(
                    ordinals,
                    report.database_parsed,
                    hit_lists,
                    report.database_loads,
                )?
            } else {
                Vec::new()
            },
        ),
        _ => (Vec::new(), Vec::new()),
    };
    // A `-subject` search's organism report looks names up in the run's database.
    let names = match report.database_reader {
        Some(database) => database.tabular_names(spec, &taxonomy)?,
        None => report
            .subject_taxdb
            .map(|taxdb| taxdb.report_names("BLASTN")),
    };
    let subjects = match (report.database_ordinals, report.database_reader) {
        (Some(ordinals), Some(reader)) => TabularSubjects::Database {
            ids: subject_ids,
            ordinals,
            oids: reader.oid_map(),
            records: report.subjects.records(),
            parsed: report.database_parsed,
            sources: report.database_sources,
            lengths: &lengths,
            taxonomy: &taxonomy,
            loads: report.database_loads,
        },
        _ => TabularSubjects::Records {
            ids: subject_ids,
            length: &record_length,
            parsed: report.database_parsed,
            records: report.subjects.records(),
        },
    };
    // Outfmt 18, the organism report: the `Query=` block shows the query record's length, also
    // with `-query_loc` (as outfmt 0's, `QueryNames`).
    if spec.kind == crate::report::tabular::TabularKind::TaxReport {
        let queries: Vec<crate::report::tax_report::TaxQuery<'_>> = hit_lists
            .iter()
            .enumerate()
            .map(|(q_idx, hit_list)| crate::report::tax_report::TaxQuery {
                title: query_titles
                    .get(q_idx)
                    .map_or(&b""[..], |title| title.as_ref()),
                length: report
                    .ranges
                    .queries
                    .placements
                    .length(q_idx, query_lengths.get(q_idx).copied().unwrap_or_default()),
                // A query of a batch that NCBI did not search has no Seq-align-set.
                hsps: (!unsearched.get(q_idx).copied().unwrap_or(false)).then(|| {
                    hit_list
                        .iter()
                        .flat_map(|hit_list| &hit_list.hsplist_array)
                        .flat_map(|hsp_list| &hsp_list.hsps)
                        .map(|hsp| crate::report::tax_report::TaxHsp {
                            subject: hsp.s_idx as usize,
                            evalue: hsp.e_value,
                            bit_score: hsp.bit_score,
                        })
                        .collect()
                }),
            })
            .collect();
        return crate::algorithm::common::tabular::write_tax_report(
            &subjects,
            report.subjects.records(),
            false,
            names.as_ref(),
            "BLASTN",
            &queries,
            writer,
            warnings,
            report.is_html,
        );
    }
    let tab = BlastnTabular {
        spec,
        queries: report.queries,
        query_records: report.query_records,
        query_ids,
        query_titles,
        searched_lengths: query_lengths,
        range_length: report.ranges.queries.range_length,
        ungapped: report.ungapped,
        database: subject_title,
        subjects,
        letters: report.subjects,
        unsearched,
        prolog: report.queries_before == 0,
        epilog: report
            .epilog
            .then_some(report.queries_before + hit_lists.len()),
        taxdb_found: report.tabular_taxdb_found,
        names: names.as_ref(),
    };
    write_blastn_tabular(&tab, hit_lists, writer, warnings, probe)
}

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::ApplicationSettings`): the data loaders decide whether
/// the readers of the subjects and the queries try a first line as a Seq-id
/// (`ReaderConfig`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
pub fn run(
    mut args: BlastnArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.id_scope = args.parse_deflines.then(Default::default);
    // The command line's progress bar (`BlastnArgs::progress_bar`; no NCBI counterpart). The
    // library API and the web entries come through `run_local` and leave it off.
    args.progress_bar = true;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:128-133
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too (`output_formats_of`). `GetExportSearchStrategyStream` comes first: it creates
    // the search strategy file (`strategy::open_strategy_file`), which an error after it
    // leaves empty, and with it a value that does not parse is left to the formatting handler
    // (`app::early_format_choice`).
    let strategy_file = crate::report::archive::strategy::open_strategy_file(
        args.export_search_strategy.as_deref(),
        args.out.as_deref(),
        "BLASTN",
    )?;
    let choice = app::early_format_choice(&args.outfmt, strategy_file.is_some())?;
    let output_formats = output_formats_of(&choice, args.out.as_deref())?;
    // The BLAST archive (outfmt 11): its request, before anything is read or written.
    let archive_requested = output_formats.contains(&BlastnOutputFormat::Archive);
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    // The run's options handle, which its archive and its search strategy share.
    let options = crate::report::archive::OptionsHandle::default();
    args.archive = archive::handle_for(
        &args,
        args.dust_filtering.as_deref(),
        archive_requested,
        &options,
    )
    .map(|handle| handle.with_outfmt_text(&args.outfmt));
    // The search strategy's request, with the archive's options handle.
    args.strategy = strategy_file.map(|file| {
        crate::report::archive::strategy::StrategyExport::for_args(
            &args,
            file,
            choice.number(),
            &options,
        )
    });
    // LOSAT's thread capability is checked before any input or output, as before.
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3631-3636
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:63-70
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2553-2557
    // The first handler opens and reads the subjects (an empty subject set fails there),
    // before the query and the output are opened and the options are checked. Files are
    // opened when a handler asks for them, and each once: a named pipe gives its bytes to
    // one reader only.
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2558-2562
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:172-176
    // The handler of the database arguments raises it, before the query and the output
    // are opened (`CStdCmdLineArgs` comes after it).
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    // With `-db` the handler of the database arguments keeps the name and reads the list: the
    // database is opened after every option handler has run (`search_cli`). With `-subject` it
    // reads the subjects here, before the query and the output are opened.
    let archive = args.archive.clone();
    let mut stderr = crate::cli::CliStderr::process();
    let mut archive_tee;
    let diagnostics: &mut (dyn std::io::Write + Send) = match &archive {
        Some(handle) => {
            archive_tee = crate::report::archive::DiagTee::new(&mut stderr, handle.clone());
            &mut archive_tee
        }
        None => &mut stderr,
    };
    // The taxid or seqid list (line DI-8): NCBI's handler of the database arguments
    // (`s_GetTaxIDList`, `CSeqDBFileGiList`) comes before the other option handlers
    // (`m_BlastDbArgs` is among the first of `CBlastpAppArgs`' `m_Args`, blastp_args.cpp:61-65)
    // and before `CStdCmdLineArgs`, which opens the query and the output: its warning before
    // theirs (probe F04 of line DI-8), and a list error leaves no `-out`, only the error archive
    // that NCBI's `PrintErrorArchive` writes for outfmt 11 (probe E of line DI-8).
    let list = match args.database_name() {
        Some(_) => match args
            .db_filter
            .resolve(&settings, "blastn", &mut *diagnostics)
        {
            Ok(list) => list,
            Err(error) => {
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        },
        None => None,
    };
    let cli_subjects = if args.db.is_some() {
        None
    } else {
        match read_cli_subjects(
            &args,
            &settings,
            &mut *diagnostics,
            archive.is_some() || args.strategy.is_some(),
        ) {
            Ok(read) => Some(read),
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        }
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:777-781
    // The output file is created here, before the options are checked, so a run that
    // stops at a check leaves an empty file; `-` is standard output.
    check_utf8_file_name(&args.query, "query", "BLASTN")?;
    // NCBI opens the query file inside `Run` too (`CStdCmdLineArgs` in its handler of the
    // arguments does not): its error is saved and an error archive follows.
    let query_file = match open_input(&args.query, "query", "BLASTN") {
        Ok(file) => file,
        Err(error) => {
            write_run_error_archive(&args, archive.as_ref(), &error);
            return Err(error);
        }
    };
    // The output file is opened once, as NCBI's stream (a named pipe gives one reader one
    // end of file).
    let out_file = open_out_file(&args)?;
    let outfmt = args.outfmt.clone();
    let pairwise = output_formats.iter().any(|format| format.is_pairwise());
    let tax_report = output_formats.iter().any(|format| {
        format
            .tabular_spec()
            .is_some_and(|spec| spec.kind == crate::report::tabular::TabularKind::TaxReport)
    });
    let mut report = crate::cli::ReportStream::of_cli(out_file, tax_report, args.strategy.as_ref());
    let strategy = args.strategy.clone();
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut report), &mut *diagnostics);
        search_cli(
            args,
            list,
            query_file,
            cli_subjects.as_ref().map(|read| {
                (
                    read.subjects(),
                    SubjectRanges {
                        placements: &read.placements,
                        records: read.records.as_deref(),
                    },
                )
            }),
            &mut outputs,
            output_formats,
            settings,
        )
    };
    // A post of the search strategy left unwritten is LOSAT's internal error (`end_run`).
    let result = crate::report::archive::strategy::end_run(strategy.as_ref(), result);
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status:
    // blastn_app.cpp:329-332).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                if let Some(lines) = archive::error_lines(error, false) {
                    let _ = crate::report::archive::write_error_archive_after_post(
                        handle,
                        error,
                        &lines,
                        &mut report,
                    );
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut report);
                }
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
    // The outfmt 0 formatter's stream throws when a write fails (oracle: `-out /dev/full`,
    // exit 6; `cli::ReportStream::finish`).
    report.finish(result, pairwise || tax_report)
}

/// The subjects of `-subject` as the handler of the database arguments reads them (`run`):
/// the searched records (their intervals with `-subject_loc`), where their letters lie in
/// the records, and the records as read (with `-subject_loc`).
/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastnArgs,
    archive: Option<&crate::report::archive::ArchiveHandle>,
    error: &anyhow::Error,
) {
    let Some(handle) = archive else {
        return;
    };
    // `-outfmt "11 qseqid"`: NCBI's `PrintErrorArchive` stops at the number before it opens `-out`.
    if !handle.error_archive_enabled() {
        return;
    }
    let Some(lines) = archive::error_lines(error, true) else {
        return;
    };
    if let Ok(out_file) = open_out_file(args) {
        let mut report = crate::cli::ReportStream::of_cli(out_file, false, None);
        let _ = crate::report::archive::write_error_archive_after_post(
            handle,
            error,
            &lines,
            &mut report,
        );
        let _ = std::io::Write::flush(&mut report);
    }
}

/// The `-out` file, created (`-`: standard output).
fn open_out_file(args: &BlastnArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    args.out
        .as_deref()
        .filter(|path| path.as_os_str() != "-")
        .map(|path| crate::cli::create_out_file(path, Some("BLASTN")).map(std::io::BufWriter::new))
        .transpose()
}

struct CliSubjects {
    subjects: Vec<FastaRecord>,
    placements: Placements,
    records: Option<Vec<FastaRecord>>,
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of `subjects`, packed (`x_read_subject_stores`;
    /// the records keep none).
    stores: Option<Vec<PackedNucleotide>>,
}

impl CliSubjects {
    /// The searched subjects.
    fn subjects(&self) -> Subjects<'_> {
        Subjects::Records(SubjectLetters::new(&self.subjects, self.stores.as_deref()))
    }
}

/// Reads the subjects of `-subject` where NCBI's handler of the database arguments reads
/// them (`run`; NCBI references there), with NCBI's error when `-subject` is missing (a
/// search without `-db`).
fn read_cli_subjects(
    args: &BlastnArgs,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut (dyn std::io::Write + Send),
    keep_deflines: bool,
) -> Result<CliSubjects> {
    let Some(subject_path) = args.subject.as_deref() else {
        // See docs/product_decisions/PD-LOSAT-CLI-NONSEARCH-DIFFERENCES.md exception 12.
        return Err(crate::cli::NativeError {
            exit: 1,
            message: format!(
                "{QUERY_OPTIONS_ERROR} Either a BLAST database or subject sequence(s) must be specified\n{}\n",
                crate::wording::USER_MANUAL_FOOTER
            ).into(),
        }
        .into());
    };
    check_utf8_file_name(subject_path, "subject", "BLASTN")?;
    let subject_file = open_input(subject_path, "subject", "BLASTN")?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range = parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // NCBI reads the records one at a time, writes each one's messages as it reads it, and
    // checks each one's range after reading it: a range that starts past the end of a
    // record stops the reading there (`read_subjects`).
    let mut subject_source = FastaInputSource::from_file(
        subject_file,
        ReaderConfig::subject("BLASTN", false, settings.data_loaders)
            .with_defline(keep_deflines)
            .with_parse_deflines(args.parse_deflines),
    );
    subject_source.set_id_scope(args.id_scope.clone());
    // EXPERIMENT (LOSAT_X_BNSTORE): without a subject range, each record's letters are packed as
    // soon as it has been read, and the records keep none (`x_read_subject_stores`).
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs); with
    // LOSAT_X_BNSTORE the packing of each record is in this interval too. No NCBI counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let (read_subjects, stores) = {
        let mut warn = |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message);
        // The archive writes the subjects' letters from the records, not from a packed store.
        if x_bnstore() && subject_range.is_none() && !keep_deflines {
            x_read_subject_stores(&mut subject_source, &mut warn)?
        } else {
            (
                read_subjects(&mut subject_source, subject_range.as_ref(), &mut warn)?,
                None,
            )
        }
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    if x_bnstore() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNSTORE", stores.is_some());
    }
    drop(subject_source);
    // NCBI reads records without residues (and intervals without letters) without a message
    // and reports them when it sets up the search (`search`).
    let ranged_subjects = cut_subjects(&read_subjects, subject_range.as_ref())?;
    let (subjects, subject_placements, full_subjects) = match ranged_subjects {
        Some((cut, placements)) => (cut, placements, Some(read_subjects)),
        None => (read_subjects, Placements::default(), None),
    };
    check_subjects_not_empty(&subjects)?;
    Ok(CliSubjects {
        subjects,
        placements: subject_placements,
        records: full_subjects,
        stores,
    })
}

/// What `search_cli` and `run_local_database` do after the processing of the options and before
/// the query: the opening of the database of `-db` (`list`: its taxid or seqid list, `None` for a
/// search without `-db`), and the two questions that its reports ask of it.
fn open_search_database(
    args: &mut BlastnArgs,
    list: Option<&crate::seqdb::DbList>,
    output_formats: &[BlastnOutputFormat],
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Option<database::BlastnDatabase>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-219
    // The database of `-db` is opened here: after the option handlers (the query file is
    // open and the output file exists) and before the test for an empty query.
    let database = match args.database_name() {
        Some(name) => {
            let database = database::BlastnDatabase::open_restricted(
                name,
                settings,
                list,
                args.db_mask.option().as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis)
            .after_calls(args.first_query_ordinal > 0);
            crate::algorithm::common::database::note_archive_database_filter(
                args.archive.as_ref(),
                args.strategy.as_ref(),
                &args.db_filter,
                database.list_taxids(),
                &database.statistics(),
                false,
                args.searchsp.is_some_and(|space| space > 0),
                args.dbsize.is_some_and(|size| size > 0),
                "BLASTN",
            )?;
            crate::algorithm::common::database::check_formats_report_sums(
                || database.pairwise_database(),
                outputs
                    .formats
                    .iter_mut()
                    .map(|format| (format.outfmt, &mut format.refusal)),
                "BLASTN",
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-133
    // A tabular field list that looks every taxid up (the plural names of the taxonomy database,
    // the organism report) asks whether `taxdb.bti` is on the search path: the plural names'
    // formatters warn without it (`algorithm::blastn::tabular`). A `-db` search asks its
    // database, a `-subject` search the search path of NCBI's databases (`subject_taxonomy`).
    match database.as_ref() {
        Some(database) => {
            if output_formats
                .iter()
                .filter_map(|format| format.tabular_spec())
                .any(|spec| spec.looks_up_every_taxid())
            {
                args.tabular_taxdb = Some(database.taxdb_found());
            }
        }
        None => subject_taxonomy(
            args,
            Some(&crate::algorithm::common::database::TaxdbSetting::of(
                settings,
            )),
            output_formats,
        )?,
    }
    // LOSAT's refusal of a database name for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for serial in output_formats.iter().filter_map(|format| format.serial()) {
            crate::report::serial::seqalign::check_database(
                serial,
                database.name(),
                Program::Blastn,
            )?;
        }
    }
    Ok(database)
}

/// What the reports `output_formats` of a `-subject` search ask of the taxonomy database,
/// through `setting` (`common::database::subject_taxonomy`): asked by the command line
/// (`open_search_database`) and by a search of records (`run_local_formats`) at the same point,
/// after the options are processed.
fn subject_taxonomy(
    args: &mut BlastnArgs,
    setting: Option<&crate::algorithm::common::database::TaxdbSetting>,
    output_formats: &[BlastnOutputFormat],
) -> Result<()> {
    let taxonomy = crate::algorithm::common::database::subject_taxonomy(
        setting,
        "BLASTN",
        output_formats
            .iter()
            .filter_map(|format| format.tabular_spec()),
    )?;
    args.subject_taxdb = taxonomy.taxdb;
    args.tabular_taxdb = taxonomy.found;
    Ok(())
}

/// The part of `run` after the output is opened: NCBI's filtering handler (`-dust`,
/// `-window_masker_db`), the
/// processing of the options, the opening of the database of `-db`, `Query is Empty!`, the
/// reading of the queries and the search. `list` is the taxid or seqid list of `-db`, read
/// before the output is opened; `subjects` are those of `-subject` (`None` with `-db`).
fn search_cli(
    mut args: BlastnArgs,
    list: Option<crate::seqdb::DbList>,
    query_file: crate::blastinput::input_files::InputFile,
    subjects: Option<(Subjects<'_>, SubjectRanges<'_>)>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.resolve_dust()?;
    args.window_masker_registry = Some(settings.window_masker_registry.clone());
    args.resolve_window_masker()?;
    args.resolve_repeat_filter(&settings);
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    let database = open_search_database(
        &mut args,
        list.as_ref(),
        &output_formats,
        &settings,
        outputs,
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:201-214
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // The query's data loaders are those of the subjects (`SDataLoaderConfig`; the BLAST
    // database of a database search is added only to the loaders' databases). The position
    // is taken on the opened file, before it is read (`stream_is_empty`). When the subjects
    // were read from standard input too, `cin` has reached its end (a failed stream), so
    // NCBI gets no position for the query either, and its reader is at its end: no query
    // batch is read.
    let shared_standard_input =
        args.query.as_os_str() == "-" && args.subject_path().as_os_str() == "-";
    // EXPERIMENT (LOSAT_X_QSTREAM): what its rule reads (`x_query_stream::Rule`): the largest query
    // batch.
    let x_rule = x_query_stream::Rule {
        batch_letters: x_query_batch_letters(&args.task),
    };
    // EXPERIMENT (LOSAT_X_QSTREAM): whether the query file is read a second time, batch by batch
    // (a second handle on it, at the position its reader starts from), where the reports may be
    // written batch by batch (LOSAT_X_RPTSTREAM, unless it keeps them after the loop for a reason
    // known before the reading).
    let x_reports_streamed = !args.parse_deflines
        && x_rpt_stream::refused_before_reading(
            x_rpt_stream::switch_on(),
            &args,
            outputs,
            &output_formats,
        )
        .is_none();
    let mut x_opening = x_query_stream::Opening::new(
        &args.query,
        &query_file,
        shared_standard_input,
        x_reports_streamed,
        x_rule,
    );
    // EXPERIMENT (LOSAT_X_BNQPACK): whether the records read before the search move into a store
    // that the batch loop takes them out of (`x_query_drop::Rule`), decided before the reading so
    // that their letters are packed as they are read.
    let x_store_rule = x_query_drop::Rule {
        formats: &output_formats,
        query_range: query_range.is_some(),
        hit_records: outputs.hits.is_some(),
        strategy: args.strategy.is_some(),
    }
    .acts();
    let x_store_acts = !x_opening.acts() && x_query_drop::on() && x_store_rule;
    let query_config = ReaderConfig::query("BLASTN", false, settings.data_loaders)
        .with_defline(args.archive.is_some() || args.strategy.is_some())
        .with_parse_deflines(args.parse_deflines)
        .with_strand(&args.strand);
    let mut query_source = FastaInputSource::from_file(query_file, query_config);
    query_source.set_id_scope(args.id_scope.clone());
    if !shared_standard_input && query_source.stream_is_empty() {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // The records are read here, and their messages are written batch by batch
    // (`run_in_pool`). What ends the reading (the end of the input, an `eEOF`, a reader
    // error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads it in,
    // after the reports of the batches before (`run_in_pool`).
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    // EXPERIMENT (LOSAT_X_QSTREAM): with the switch acting this reading keeps only each record's
    // length and how the reading ends; the batches read their records again (`x_query_stream`).
    let mut x_packed = Vec::new();
    // EXPERIMENT (LOSAT_X_BNQPARTS): the parts of the store, and whether they were read in parts.
    let mut x_parts = None;
    let mut x_parts_read_in_parts = false;
    let (records, end, x_store) = if x_opening.acts() {
        let (store, end) = x_query_stream::Store::read(
            x_opening,
            &mut query_source,
            query_config,
            query_range.as_ref(),
            args.num_threads,
        )?;
        (Vec::new(), end, store)
    } else if x_store_acts && x_query_drop::parts_on() {
        // EXPERIMENT (LOSAT_X_BNQPARTS): the records read into the parts of the store, on the
        // search's threads where the file can be cut into parts (each part read to its end: the
        // reading ended at the end), else by the search's reader.
        let (parts, end) = match x_opening.store_parts(query_config, args.num_threads) {
            Some(parts) => {
                x_parts_read_in_parts = true;
                (parts, QueryEnd::Input)
            }
            None => x_query_drop::read_parts(&mut query_source, query_config.id_prefix),
        };
        x_parts = Some(parts);
        (Vec::new(), end, x_query_stream::Store::off(x_opening))
    } else if x_store_acts {
        let (QueryRecords { records, end }, packed) =
            crate::blastinput::fasta_reader::x_read_queries_packed(&mut query_source);
        x_packed = packed;
        (records, end, x_query_stream::Store::off(x_opening))
    } else {
        let QueryRecords { records, end } = read_queries(&mut query_source);
        (records, end, x_query_stream::Store::off(x_opening))
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    // A database has no subject ranges (`-subject_loc` excludes `-db`).
    let no_placements = Placements::default();
    let (subjects, subject_ranges) = match (subjects, database.as_ref()) {
        (Some((subjects, ranges)), _) => (subjects, ranges),
        (None, Some(database)) => (
            Subjects::Database(database),
            SubjectRanges {
                placements: &no_placements,
                records: None,
            },
        ),
        (None, None) => unreachable!("`run` reads the subjects of a search without -db"),
    };
    // EXPERIMENT (LOSAT_X_BNQPACK): where no report reads the queries' letters, the records read
    // before the search move into a store that the batch loop takes them out of, so a batch can drop
    // its letters once it has encoded them (`x_query_drop.rs`). The reach of the switch is recorded
    // here, after the reading, for the searches that reach it.
    if !x_store.acts() && x_query_drop::on() {
        crate::utils::xreach::x_reach!("LOSAT_X_BNQPACK", x_store_rule);
    }
    let x_taken_store;
    let queries = if x_store.acts() {
        SearchQueries::Streamed(&x_store)
    } else if let Some(parts) = x_parts {
        x_taken_store =
            x_query_drop::Store::from_parts(parts, query_config.id_prefix, x_parts_read_in_parts);
        SearchQueries::Taken(&x_taken_store)
    } else if x_store_acts {
        x_taken_store = x_query_drop::Store::new(records, x_packed);
        SearchQueries::Taken(&x_taken_store)
    } else {
        SearchQueries::Records(&records)
    };
    search(
        args,
        queries,
        reading,
        Ranges {
            query: query_range,
            subjects: subject_ranges,
        },
        subjects,
        outputs,
        output_formats,
        None,
    )
}

/// The query range, read by NCBI's query options handler: after the filtering handler
/// (`-dust`), before the formatting handler and the check of the options.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blastn_args.cpp:100-102
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:1995-1999
fn parse_query_range(args: &BlastnArgs) -> Result<Option<SequenceRange>> {
    parse_optional_range(args.query_loc.as_deref(), RangeRole::Query)
}

/// The ranges of a search: the query range as given, and the subjects' ranges.
struct Ranges<'a> {
    query: Option<SequenceRange>,
    subjects: SubjectRanges<'a>,
}

/// Where each subject's searched letters lie in its record, and the records as read
/// (with `-subject_loc`; the report shows their lowercase letters).
#[derive(Clone, Copy)]
struct SubjectRanges<'a> {
    placements: &'a Placements,
    records: Option<&'a [FastaRecord]>,
}

/// The ranges as the report uses them: the query input and the subjects' ranges.
struct ReportRanges<'a> {
    queries: &'a QueryInput,
    subjects: SubjectRanges<'a>,
}

/// NCBI's error for a subject file without records, raised when it reads the subjects
/// (the first option handler, `run`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:225-227
pub(crate) fn check_subjects_not_empty<R>(subjects: &[R]) -> Result<()> {
    if subjects.is_empty() {
        return Err(crate::cli::NativeError {
            exit: 3,
            message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
        }
        .into());
    }
    Ok(())
}

/// The report of the CLI's `-outfmt`, read as NCBI reads it before its option handlers run
/// (`run`; `blastinput/app.rs` `app::early_format_choice`, `report_format`). A format that
/// NCBI itself rejects later has no report: the run stops at NCBI's error, the formatting
/// handler's in `process_options` (AIRR rearrangement, FASTA) or the XInclude formats' in
/// `search` (13 and 14 to standard output). Nor has a value that does not parse where the argv
/// names a search strategy file (`strategy`): the formatting handler reads it again
/// (`process_options`) and fails there, after the query options.
fn cli_output_formats(
    outfmt: &str,
    out: Option<&std::path::Path>,
    strategy: bool,
) -> Result<Vec<BlastnOutputFormat>> {
    output_formats_of(&app::early_format_choice(outfmt, strategy)?, out)
}

/// The report of the choice `choice` that `app::early_format_choice` read (`cli_output_formats`):
/// none for one that NCBI itself rejects later, nor for one it left to the formatting handler.
fn output_formats_of(
    choice: &app::EarlyFormat,
    out: Option<&std::path::Path>,
) -> Result<Vec<BlastnOutputFormat>> {
    let Some(choice) = choice.parsed() else {
        return Ok(Vec::new());
    };
    app::report_format(choice, Program::Blastn, out)?
        .map(|format| BlastnOutputFormat::of(format, choice))
        .into_iter()
        .collect()
}

/// The requested output formats of `run_local`, each read as the CLI reads its single
/// `-outfmt` (`cli_output_formats`; `strategy`: the argv names a search strategy file): NCBI's
/// parse and LOSAT's refusals, before `-out` is opened. A format without a report
/// (`app::report_format`: one that NCBI rejects later) is left out, as on the command line: NCBI's
/// formatting handler rejects it after `-out` is opened (`process_options`), or `search` (the
/// XInclude formats without a file). `out`: the file of `api::run_to_file`, or the `-out` of the
/// caller's receiver of the side files of a format (`FormatOutput::side_out`).
fn parse_output_formats<'a>(
    formats: impl IntoIterator<Item = (&'a str, Option<&'a std::path::Path>)>,
    strategy: bool,
) -> Result<Vec<BlastnOutputFormat>> {
    let mut output_formats = Vec::new();
    for (outfmt, out) in formats {
        output_formats.extend(cli_output_formats(outfmt, out, strategy)?);
    }
    Ok(output_formats)
}

/// NCBI's processing of the options after the query and the output are opened
/// (`SetOptions`): the formatting handler (its errors for the formats of other programs,
/// then its few-matches warning), and the check of the options.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2874-2886
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2975-2977
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3636
/// The hit list size is the -max_target_seqs value, or 500 when it is omitted.
fn process_options(args: &mut BlastnArgs, outputs: &mut ReportOutputs<'_>) -> Result<()> {
    // The formatting handler comes after the query options handler (blastn_args.cpp:100-105).
    let mut outfmts = Vec::with_capacity(outputs.formats.len());
    for format in &outputs.formats {
        let choice = app::parse_formatting_string(format.outfmt)?;
        app::formatting_handler_check(&choice, Program::Blastn)?;
        outfmts.push(choice.number);
    }
    args.html_links = app::html_links(args.html, Program::Blastn, args.database_name(), &outfmts)?;
    // The hit list size and the counts of the pairwise report, with the formatting
    // handler's warnings (`ExtractAlgorithmOptions`), once for the run.
    let (counts, warnings) = crate::algorithm::blastn::scoring::report_counts_of(args, &outfmts)?;
    args.report_counts = counts;
    outputs.diagnostics.write_all(&warnings)?;
    check_scoring_options(args)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
// NCBI formats one result set without searching again; several requested formats
// are several CBlastFormat printers over the same result set.
/// Runs one BLASTN search over already parsed records and writes every requested
/// output format from the same result (the shared entry of the CLI, web ABI v1 and
/// v2).
///
/// Each requested `-outfmt` is parsed as on the command line, before the search starts.
/// No BLASTN search option depends on the output formats it supports (6 and 7). The
/// `-query` and `-subject` values of `args` are used only as display names. The hit
/// records (`ReportOutputs::hits`) are not produced yet.
///
/// The BLAST archive (outfmt 11) writes each record's DefLine as the reader kept it
/// (`ReaderConfig::with_defline(true)`). Records without a kept defline get `>` and their
/// title, which differs from NCBI for a defline with extra spaces or control bytes.
pub fn run_local(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    crate::algorithm::blastn::args::BLASTN_TASK.check(&mut args.task)?;
    // The BLAST archive (outfmt 11) among the formats: its request and a diagnostics stream
    // that saves the messages.
    let archive_requested = outputs.formats.iter().any(|format| {
        app::parse_formatting_string(format.outfmt)
            .is_ok_and(|choice| choice.number == app::output_format::ARCHIVE_FORMAT)
    });
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    if archive_requested && args.db.is_none() {
        QueryStrand::parse(&args.strand)?;
        args.archive = archive::handle_for(
            &args,
            args.dust_filtering.as_deref(),
            true,
            &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
        )
        .map(|handle| outputs.with_archive_outfmt(handle));
    }
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        run_local_formats(args, query_records, query_end, subject_records, outputs)
    })
}

fn run_local_formats(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    subject_records: &[FastaRecord],
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // `run_local` searches subject records (the library, `api::run`; the WebAssembly adapter,
    // `web/adapter`; ABI v1, which refuses `-db` earlier: `parse_blastn_args` in `web_api.rs`);
    // a `-db` search runs through `run_local_database` (or the command line's `run`).
    if args.db.is_some() {
        anyhow::bail!(crate::api::local_blast::RUN_LOCAL_DB_REFUSAL);
    }
    // The counts file of `-window_masker_db` is read through the caller's counts source
    // (`BlastnArgs::window_masker_counts`: `api::run` gives `FileCounts`; `web/adapter` gives
    // its registered counts files); a caller of `run_local` that sets none has no counts
    // source; ABI v1 refuses `-window_masker_db` earlier (`parse_blastn_args` in `web_api.rs`).
    if args.window_masker_db.is_some() && args.window_masker_counts.is_none() {
        anyhow::bail!(
            "-window_masker_db (a WindowMasker counts file) needs its counts source set before calling run_local (BlastnArgs::window_masker_counts)"
        );
    }
    // The counts file of `-window_masker_taxid` (non-zero) is looked up through the caller's
    // registry (`BlastnArgs::window_masker_registry`, from the run's settings) and counts
    // source (`BlastnArgs::window_masker_counts`: `api::run` gives `FileCounts`; `web/adapter`
    // gives its registered counts files); a caller of `run_local` that omits either cannot
    // resolve the taxid; ABI v1 refuses `-window_masker_taxid` earlier (`parse_blastn_args` in
    // `web_api.rs`).
    if args.window_masker_taxid.is_some_and(|taxid| taxid != 0)
        && (args.window_masker_registry.is_none() || args.window_masker_counts.is_none())
    {
        anyhow::bail!(
            "-window_masker_taxid (a directory of WindowMasker counts files) needs its registry and counts source set before calling run_local (BlastnArgs::window_masker_registry, BlastnArgs::window_masker_counts)"
        );
    }
    // The database of `-filtering_db` is looked up with the caller's settings (the library gives
    // them, `api::run`; the WebAssembly adapter, `web_adapter::resolve_blastn_repeat_filter`);
    // a caller of `run_local` that has not resolved it (`BlastnArgs::resolve_repeat_filter`) has
    // no repeat search; ABI v1 refuses the option earlier (`parse_blastn_args` in `web_api.rs`).
    if args.filtering_db.is_some() && args.repeat_filter.is_none() {
        anyhow::bail!(
            "-filtering_db (a BLAST database of repeats) needs its repeat filter resolved with the run's settings (BlastnArgs::resolve_repeat_filter)"
        );
    }
    // The command line has checked `-strand` when it read it (NCBI's argument reader); a caller of
    // `run_local` has not.
    QueryStrand::parse(&args.strand)?;
    let output_formats = parse_output_formats(
        outputs
            .formats
            .iter()
            .map(|format| (format.outfmt, format.side_out(args.out.as_deref()))),
        args.export_search_strategy.is_some(),
    )?;
    // The thread count is checked before the inputs are read, as on the command line (`run`).
    crate::utils::threading::validate_threads(args.num_threads)?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)?;
    for record in &subject_records[..subjects_read(subject_records, subject_range.as_ref())] {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let ranged_subjects = cut_subjects(subject_records, subject_range.as_ref())?;
    let (searched_subjects, subject_placements) = match &ranged_subjects {
        Some((records, placements)) => (records.as_slice(), placements.clone()),
        None => (subject_records, Placements::default()),
    };
    check_subjects_not_empty(searched_subjects)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    args.resolve_dust()?;
    args.resolve_window_masker()?;
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    // The taxonomy database of the reports, where the command line asks (`open_search_database`).
    let setting = args.taxdb_setting.take();
    subject_taxonomy(&mut args, setting.as_ref(), &output_formats)?;
    search(
        args,
        SearchQueries::Records(query_records),
        query_end.reading(),
        Ranges {
            query: query_range,
            subjects: SubjectRanges {
                placements: &subject_placements,
                records: ranged_subjects.is_some().then_some(subject_records),
            },
        },
        Subjects::Records(SubjectLetters::new(searched_subjects, None)),
        outputs,
        output_formats,
        None,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:191-197
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:327-332
/// Runs one BLASTN search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `run`'s and `search_cli`'s, in their order; what differs is what they
/// read: the records are read already, and the outputs are the caller's. An archive (outfmt 11)
/// ends as the command line's does: an error archive follows a failure and messages that no
/// archive carried.
pub fn run_local_database(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    crate::algorithm::blastn::args::BLASTN_TASK.check(&mut args.task)?;
    // The command line has checked `-strand` when it read it (NCBI's argument reader); a caller of
    // `run_local` has not.
    QueryStrand::parse(&args.strand)?;
    let output_formats = parse_output_formats(
        outputs
            .formats
            .iter()
            .map(|format| (format.outfmt, format.side_out(args.out.as_deref()))),
        args.export_search_strategy.is_some(),
    )?;
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    args.archive = archive::handle_for(
        &args,
        args.dust_filtering.as_deref(),
        output_formats.contains(&BlastnOutputFormat::Archive),
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads)?;
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input();
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        search_database(
            args,
            query_records,
            query_end,
            output_formats,
            settings,
            outputs,
            None,
        )
    })
}

/// The repeat search of blastn's `-filtering_db` (`repeat_filter`): `args` (the search's options,
/// `RepeatFilter::search_args`) over `query_records` against the filter's database, in `pool`, the
/// pool of the search whose queries it masks, without a report: the hit lists, as the search hands
/// them to a library caller (`ReportOutputs::hits`). Its messages are dropped, as NCBI reads only
/// the alignments of its search of the repeats, and nothing that only a formatter reads is read
/// (`BlastnArgs::hits_only`: NCBI's search of the repeats has no formatter). An exception of NCBI's
/// search is `RepeatError::Posted`, its message: an error of the database reader where the
/// database is opened, and the stop of a batch that the search's set-up finds without a valid query
/// (`BLAST engine error: `, `CBlastException` eSetup).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/repeats_filter_cxx.cpp:368-390
pub(crate) fn repeat_search_hits(
    args: BlastnArgs,
    query_records: &[FastaRecord],
    filter: &super::super::repeat_filter::RepeatFilter,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> std::result::Result<Vec<Arc<Vec<PairwiseHit>>>, super::super::repeat_filter::RepeatError> {
    use super::super::repeat_filter::RepeatError;
    use crate::algorithm::common::database::DatabaseOpenError;
    // The search's `CSeqDB`, without a list, meets the errors of the alias files' counts where it
    // opens the database (`common::database::restrict_database`).
    let opened = crate::algorithm::common::database::open_search_database(
        filter.database(),
        crate::seqdb::SeqType::Nucleotide,
        filter.settings(),
    );
    let failed = match opened {
        Ok(database) => database.check_alias_values().err(),
        Err(DatabaseOpenError::Reader(error)) => Some(error),
        Err(DatabaseOpenError::Registry(_)) => None,
    };
    if let Some(error) = failed {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-639
        // The set-up of each query posts the exception's `GetMsg()`.
        match (error.ncbi_exception(), error.ncbi_message()) {
            (Some(exception), _) => return Err(RepeatError::Posted(exception.get_msg().to_vec())),
            (None, Some(text)) => return Err(RepeatError::Posted(text)),
            (None, None) => {}
        }
    }
    // The search writes a tabular report, which nothing reads (as the library's batches do,
    // `api::batches`): the hits come from the hit lists.
    let mut diagnostics: Vec<u8> = Vec::new();
    let mut report = std::io::sink();
    let mut collected: Vec<Arc<Vec<PairwiseHit>>> = Vec::new();
    let result = {
        let mut on_hits = |hits| collected.push(hits);
        let mut outputs = ReportOutputs {
            other_outfmts: &[],
            html_posts: None,
            formats: vec![crate::api::local_blast::FormatOutput {
                outfmt: "6",
                refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                sink: crate::api::local_blast::OutputSink::Writer(&mut report),
            }],
            diagnostics: &mut diagnostics,
            later_warnings: None,
            hits: Some(&mut on_hits),
            observer: None,
        };
        parse_output_formats(std::iter::once(("6", None)), false).and_then(|output_formats| {
            search_database(
                args,
                query_records,
                // The repeat search's queries are the batch's records, read whole.
                &crate::api::QueryEnd::Input,
                output_formats,
                filter.settings().clone(),
                &mut outputs,
                Some(pool),
            )
        })
    };
    result.map_err(|error| match error.downcast::<crate::cli::NativeError>() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:650-651
        // The search's set-up that finds no valid query stops with its messages, which the
        // command line writes after its prefix (`ENGINE_ERROR`, NCBI's `BLAST engine error:`).
        Ok(native) => match native
            .message
            .strip_prefix(ENGINE_ERROR.as_bytes())
            .and_then(|text| text.strip_prefix(b" "))
            .and_then(|text| text.strip_suffix(b"\n"))
        {
            Some(text) if native.exit == 3 => RepeatError::Posted(text.to_vec()),
            _ => RepeatError::Search(native.into()),
        },
        Err(error) => RepeatError::Search(error),
    })?;
    Ok(collected)
}

/// `run_local_database` from the taxid or seqid list on, in the pool `within` when one is given
/// (`repeat_search_hits`).
fn search_database(
    mut args: BlastnArgs,
    query_records: &[FastaRecord],
    query_end: &crate::api::QueryEnd,
    output_formats: Vec<BlastnOutputFormat>,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
    within: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    // The handler of the database arguments reads the list before the other handlers (`run`).
    let list = args
        .db_filter
        .resolve(&settings, "blastn", &mut *outputs.diagnostics)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    args.resolve_dust()?;
    args.window_masker_registry = Some(settings.window_masker_registry.clone());
    args.resolve_window_masker()?;
    args.resolve_repeat_filter(&settings);
    let query_range = parse_query_range(&args)?;
    process_options(&mut args, outputs)?;
    let database = open_search_database(
        &mut args,
        list.as_ref(),
        &output_formats,
        &settings,
        outputs,
    )?
    .ok_or_else(app::missing_subject_error)?;
    let no_placements = Placements::default();
    search(
        args,
        SearchQueries::Records(query_records),
        query_end.reading(),
        Ranges {
            query: query_range,
            subjects: SubjectRanges {
                placements: &no_placements,
                records: None,
            },
        },
        Subjects::Database(&database),
        outputs,
        output_formats,
        within,
    )
}

/// The queries of a search: records read before it, or EXPERIMENT (LOSAT_X_QSTREAM) a query file
/// whose batches read their records when the batch loop needs them (`x_query_stream`).
#[derive(Clone, Copy)]
enum SearchQueries<'a> {
    Records(&'a [FastaRecord]),
    Streamed(&'a x_query_stream::Store),
    /// EXPERIMENT (LOSAT_X_BNQPACK): records read before the search that the batch loop takes out
    /// of their store (`x_query_drop.rs`).
    Taken(&'a x_query_drop::Store),
}

#[allow(clippy::too_many_arguments)]
fn search(
    mut args: BlastnArgs,
    queries: SearchQueries<'_>,
    reading: QueryReading,
    ranges: Ranges<'_>,
    subjects: Subjects<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    within: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    // EXPERIMENT (LOSAT_X_QSTREAM): a streamed search has no records here (its reading is never
    // `Records`), and its query input comes from the lengths of the first reading; EXPERIMENT
    // (LOSAT_X_BNQPACK) a taken one has its records in its store, with its query input.
    let query_records: &[FastaRecord] = match queries {
        SearchQueries::Records(records) => records,
        SearchQueries::Streamed(_) | SearchQueries::Taken(_) => &[],
    };
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        outputs
            .diagnostics
            .write_all(b"Warning: [blastn] Query is Empty!\n")?;
        return Ok(());
    }
    // The database's MegaBLAST index that `-use_index true` asks for, after `Query is Empty!`
    // and before the formatter (`mb_index::megablast_index`); LOSAT's refusal of an index comes
    // with its other limits below.
    let megablast_index = match subjects {
        Subjects::Database(database) => {
            mb_index::megablast_index(&args, database, outputs.diagnostics)?
        }
        Subjects::Records(_) => mb_index::MegablastIndex::NotUsed,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
    // `blastinput::db_mask::DbMaskOption::resolve`)
    // The formatter resolves the database masking (`-db_soft_mask`, `-db_hard_mask`): its warning
    // (one diagnostic, which the archive saves as one message) or its database error, before its
    // other work; the archive and the search strategy write the result's program options.
    if let (Subjects::Database(database), Some(option)) = (subjects, args.db_mask.option()) {
        let archive = args.archive.clone();
        let diagnostics = &mut *outputs.diagnostics;
        let resolved = database.resolve_mask_option(&option, &mut |text| match &archive {
            Some(handle) => handle.write_warning(diagnostics, text),
            None => diagnostics.write_all(text),
        })?;
        crate::algorithm::common::database::note_archive_database_masking(
            args.archive.as_ref(),
            args.strategy.as_ref(),
            resolved
                .as_ref()
                .map(crate::blastinput::db_mask::ResolvedMask::archive_params)
                .unwrap_or_default(),
            option.strategy_params(resolved.as_ref()),
        )?;
        args.db_mask_resolved = resolved;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-136
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // The formatter, made after `Query is Empty!` and before the outfmt 0 prolog, sets up
    // the subjects: a subject without letters (a record without residues, or an interval
    // that starts just past its record's end) gets its warning, stays in the database
    // statistics with no letters, and never hits. (A database is not set up so.)
    if let Subjects::Records(subject_set) = subjects {
        write_empty_subject_warnings_by_length(
            subject_set.records(),
            |index| subject_set.len(index),
            "blastn",
            outputs.diagnostics,
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:255-258
    // The XInclude formats written to standard output fail after the formatter is made (a
    // format 13 or 14 to a file takes the file name as the base of its side files).
    outputs.xinclude_check(args.out.as_deref())?;
    // NCBI reads the queries after `Query is Empty!`, one batch at a time, with its
    // reader's warnings (`run_in_pool`).
    // NCBI's x_CleanAndCompress reads past the end of some titles of punctuation (NCBI
    // crashes when it writes the title of such a subject with hits); LOSAT stops at the end
    // of the title (`report/defline.rs`, approved exception 2 of PD-LOSAT-NCBI-DEFECTS).
    // LOSAT's limits come where NCBI starts the search, after its checks and its
    // `Query is Empty!` success.
    check_losat_limits(&args)?;
    args.megablast_index_loaded = megablast_index.check()?;
    if args.megablast_index_loaded.is_some() {
        mb_index::note_index_use(args.archive.as_ref(), args.strategy.as_ref())?;
    }
    if let (Subjects::Database(database), Some(mask)) = (subjects, args.db_mask_resolved.as_ref()) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1625-1633
        // The words that the scan reads: the discontiguous template's, else the word size.
        let template_length = super::super::coordination::determine_template(&args).1;
        let scanned_word = if template_length > 0 {
            usize::from(template_length)
        } else {
            super::super::coordination::determine_effective_word_size(&args)
        };
        database.check_mask_limits(mask, scanned_word)?;
    }
    app::check_search_environment("BLASTN", !matches!(subjects, Subjects::Database(_)))?;
    // The variables' error stops the run after the outfmt 0 prolog (`run_in_pool`).
    let batching = query_batching_from_environment(
        super::super::coordination::task_uses_megablast_chunks(&args.task),
    )
    .map(|(mut batching, overlap_error)| {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:263-264
        // The caller's records are one batch (`BlastnArgs::one_query_batch`): a size that the
        // records do not reach, fixed as `BATCH_SIZE` fixes it.
        if args.one_query_batch {
            batching.fixed_batch_size = i32::MAX;
            batching.batch_size = i32::MAX as u32;
        }
        (batching, overlap_error)
    });
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`cut_queries`): a record is searched cut, and one whose interval starts just past
    // its end is a query without data (`run_in_pool`).
    let ranged_queries = ranges
        .query
        .as_ref()
        .filter(|_| matches!(queries, SearchQueries::Records(_)))
        .map(|range| cut_queries(query_records, range));
    let input_records = query_records;
    let (query_records, query_input) = match (ranged_queries, queries) {
        (Some(ranged), _) => (
            std::borrow::Cow::Owned(ranged.records),
            std::borrow::Cow::Owned(ranged.input),
        ),
        (None, SearchQueries::Records(_)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Owned(QueryInput::whole(query_records)),
        ),
        (None, SearchQueries::Streamed(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
        (None, SearchQueries::Taken(store)) => (
            std::borrow::Cow::Borrowed(query_records),
            std::borrow::Cow::Borrowed(store.query_input()),
        ),
    };
    let queries = match queries {
        SearchQueries::Records(_) => x_query_stream::Queries::Records {
            input: input_records,
            searched: &query_records,
        },
        SearchQueries::Streamed(store) => x_query_stream::Queries::Streamed(store),
        SearchQueries::Taken(store) => x_query_stream::Queries::Taken(store),
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-138
    // The report and the query batches (`GetDbTotalLength`) take the total through an
    // `int`, which LOSAT does not reproduce beyond its range. (A database's total is the
    // `Int8` of `GetBlastDbInfo`, blast_format.cpp:140-154.)
    if let Subjects::Records(subject_set) = subjects {
        let subject_letters: usize = (0..subject_set.records().len())
            .map(|index| subject_set.len(index))
            .sum();
        if subject_letters > i32::MAX as usize {
            anyhow::bail!(
                "the subjects have {subject_letters} letters; NCBI BLAST+ reports a total of 2^31 letters or more through a 32-bit int, which is not supported by LOSAT's BLASTN"
            );
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // The repeat search of `-filtering_db` runs in the pool of the search whose queries it masks
    // (`repeat_search_hits`); every other search has its own.
    let requested = args.num_threads;
    let run = |pool: &crate::utils::threading::SearchPool<'_>| {
        run_in_pool(
            args,
            queries,
            reading,
            &query_input,
            subjects,
            ranges.subjects,
            outputs,
            output_formats,
            pool,
            batching,
        )
    };
    let x_result = match within {
        Some(pool) => run(pool),
        None => crate::utils::threading::with_search_pool(requested, "blastn", run),
    };
    // No NCBI counterpart: the PAIRPAR counters with LOSAT_TIMING; nothing without it.
    x_pair_seed::summary();
    x_result
}

/// An error of a query batch after the first that does not depend on the batch's queries
/// (the Karlin-Altschul error of the scoring options, or the average subject length;
/// `search_query_batch`, in `context_statistics`): the first batch with a valid query meets it, and NCBI raises it
/// after the reports of the batches before, whose queries were all invalid.
#[derive(Debug)]
struct ErrorAfterInvalidBatches(anyhow::Error);

impl std::fmt::Display for ErrorAfterInvalidBatches {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl std::error::Error for ErrorAfterInvalidBatches {}

// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:295-312
/// A database error of a query batch's search (`database::search_blocks`,
/// `database::traceback_listed`) as NCBI reports it: a sequence file that the preliminary
/// search cannot read is NCBI's `BLAST engine error` (`database::PreliminaryReadError`,
/// local_blast.cpp:240-251), one that the traceback cannot read its `BLAST Database
/// error`. NCBI raises it in the batch's search, after the outfmt 0 prolog and the reports
/// of the batches before (`ErrorAfterInvalidBatches`; for the first batch, the prolog of a
/// file sink is written here, the others were written before the search). LOSAT's
/// refusals are errors of their own.
fn database_search_error(
    error: anyhow::Error,
    first_batch: bool,
    outputs: BatchOutputs<'_, '_>,
    intro: BlastnIntro,
    shown: &ShownSubjects,
) -> anyhow::Error {
    let error = match error.downcast::<database::PreliminaryReadError>() {
        Ok(read_error) => crate::blastinput::app::engine_error(read_error.0),
        Err(error) => error,
    };
    if !error.is::<crate::cli::NativeError>() {
        return error;
    }
    if !first_batch {
        return ErrorAfterInvalidBatches(error).into();
    }
    let (outputs, output_formats) = outputs.expect(FIRST_BATCH_OUTPUTS);
    match write_pairwise_prologs(outputs, output_formats, intro, shown, false) {
        Ok(()) => error,
        Err(write_error) => write_error,
    }
}

/// NCBI's query batch sizes of a search of `subjects`, one batch at a time, as `run_in_pool`
/// takes them: for a caller that hands over one batch at a time (`api::run_batches`,
/// `BlastnArgs::one_query_batch`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-270
pub(crate) struct QueryBatchSizes {
    batching: QueryBatching,
    mixer: BatchSizeMixer,
}

impl QueryBatchSizes {
    pub(crate) fn new(args: &BlastnArgs, subjects: &[FastaRecord]) -> Result<Self> {
        use crate::blastinput::seq_range::{record_interval, RecordInterval};
        let subject_range = parse_optional_range(args.subject_loc.as_deref(), RangeRole::Subject)?;
        let total_length: i64 = subjects
            .iter()
            .map(|record| match subject_range.as_ref() {
                Some(range) => match record_interval(range, record.seq().len()) {
                    RecordInterval::Letters { from, to } => (to - from + 1) as i64,
                    RecordInterval::Empty | RecordInterval::PastEnd => 0,
                },
                None => record.seq().len() as i64,
            })
            .sum();
        Self::of_total_length(args, total_length)
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:338-345
    /// The sizes of a search of the database `name`. A database's total is the `Int8` of its index
    /// files (`GetDbTotalLength` sums `m_DbInfo`, `GetBlastDbInfo`: the whole database's, also with
    /// a list, as `run_in_pool` reads it). A database that does not open has no total here: the
    /// search of the first batch fails with NCBI's error.
    pub(crate) fn of_database(
        args: &BlastnArgs,
        name: &str,
        settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    ) -> Result<Self> {
        let total_length = database::BlastnDatabase::open(name, settings)
            .map_or(0, |database| database.pairwise_database().total_letters);
        Self::of_total_length(args, total_length)
    }

    fn of_total_length(args: &BlastnArgs, total_length: i64) -> Result<Self> {
        let (mut batching, _) = query_batching_from_environment(
            super::super::coordination::task_uses_megablast_chunks(&args.task),
        )?;
        let mut mixer = BatchSizeMixer::new(batching.split.mixer_max_batch_size());
        if batching.fixed_batch_size == 0 {
            if total_length > 0 {
                mixer.set_target_hits((total_length / 3000) as i32);
            }
            batching.batch_size = mixer.batch_size(None) as u32;
        }
        Ok(Self { batching, mixer })
    }

    /// The size of the next batch (`TSeqPos`).
    pub(crate) fn batch_size(&self) -> u32 {
        self.batching.batch_size
    }

    /// The batch just searched made `extensions` successful initial extensions.
    pub(crate) fn searched(&mut self, extensions: i32) {
        if self.batching.fixed_batch_size == 0 {
            self.batching.batch_size = self.mixer.batch_size(Some(extensions)) as u32;
        }
    }
}

/// NCBI's query batch size (`GetQueryBatchSize`, 0 for the adaptive `CBatchSizeMixer`) and
/// query chunk sizes, as a query batch of the search uses them.
#[derive(Clone, Copy, Debug)]
struct QueryBatching {
    /// The `BATCH_SIZE` value; 0 when it is not set.
    fixed_batch_size: i32,
    /// The size of the current batch (`TSeqPos`).
    batch_size: u32,
    split: SplitSizes,
}

/// Reads the environment variables that set NCBI's query batches and chunks, in NCBI's
/// order: `CHUNK_SIZE` (the `CBatchSizeMixer` maximum), `BATCH_SIZE`, then
/// `OVERLAP_CHUNK_SIZE` (the query splitter of the first batch).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-262
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:55-59
/// A `CHUNK_SIZE` or `BATCH_SIZE` that `NStr::StringToInt` cannot convert is the `Err`: its
/// `CStringException` (`app::environment_int`), which stops NCBI after the outfmt 0 prolog and
/// before the first batch is read (`run_in_pool`). `OVERLAP_CHUNK_SIZE` is read by the query
/// splitter of each batch, after the batch is read: one that it cannot convert is the second
/// value, its exception, which stops the first batch there (`run_in_pool`); no batch is split
/// with it, so the sizes count it as unset. The batch sizes alone (`QueryBatchSizes`) do not
/// depend on it. A negative `CHUNK_SIZE` is a `size_t` near 2^64; LOSAT rejects it where it splits
/// a query batch (`search_query_chunks`), and otherwise searches as NCBI.
fn query_batching_from_environment(megablast: bool) -> Result<(QueryBatching, Result<()>)> {
    use crate::blastinput::app::read_query_split_variable;
    let chunk_size = read_query_split_variable("CHUNK_SIZE")?;
    let fixed_batch_size = crate::blastinput::query_batch::get_query_batch_size(
        std::env::var_os("BATCH_SIZE").as_deref(),
    )
    .map_err(|e| crate::cli::NativeError::exception(&e))?;
    let (overlap, overlap_error) = match read_query_split_variable("OVERLAP_CHUNK_SIZE") {
        Ok(overlap) => (overlap, Ok(())),
        Err(error) => (None, Err(error)),
    };
    let split = SplitSizes::new(megablast, chunk_size, overlap);
    Ok((
        QueryBatching {
            fixed_batch_size,
            batch_size: fixed_batch_size as u32,
            split,
        },
        overlap_error,
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
/// The subjects of a search, read and encoded once for all query batches, as NCBI's
/// `CLocalDbAdapter` of the subject set serves every batch. Only the packed (ncbi2na)
/// letters that the preliminary stage reads are kept for every subject; the traceback
/// makes the blastna letters of a subject with preliminary HSPs from its record.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
struct PreparedSubjects<'a> {
    records: &'a [FastaRecord],
    /// EXPERIMENT (LOSAT_X_BNSTORE): the letters of the records, packed (the records keep
    /// none); the preliminary stage reads their packed bytes.
    stores: Option<&'a [PackedNucleotide]>,
    metadata: SubjectMetadata,
    /// The packed letters of each record; empty when `stores` holds them.
    packed: Vec<Vec<u8>>,
}

impl<'a> PreparedSubjects<'a> {
    /// The records and their letters.
    fn letters(&self) -> SubjectLetters<'a> {
        SubjectLetters::new(self.records, self.stores)
    }

    /// Subject `index` as the per-subject kernel reads it.
    fn subject(&self, index: usize) -> BlastnSubjectRef<'_> {
        match self.stores {
            Some(stores) => BlastnSubjectRef::packed_record(&stores[index]),
            None => BlastnSubjectRef::record(&self.records[index], &self.packed[index]),
        }
    }
}

/// The subjects of a search as its caller gives them.
#[derive(Clone, Copy)]
enum Subjects<'a> {
    /// A FASTA subject set (`-subject`, or the records of `run_local`): the searched
    /// records (their intervals with `-subject_loc`), and their letters (with
    /// LOSAT_X_BNSTORE, packed).
    Records(SubjectLetters<'a>),
    /// An NCBI BLAST database (`-db`).
    Database(&'a database::BlastnDatabase),
}

/// The steps of NCBI's preliminary stage that finish the HSP list of one subject of an ungapped
/// (`-ungapped`) search (`Blast_HSPListReevaluateUngapped`, the e-values, the reap by the
/// preliminary e-value and by query coverage): the HSPs are re-evaluated with the ambiguity
/// codes of the subject (`s_blastna`, whole, without sentinel bytes), the ones below the
/// word cutoff of their context or the `-perc_identity` are deleted, the list is sorted by
/// score again, and the ones above the e-value (`evalue_of`) are deleted. The reap by the
/// query coverage (`-qcov_hsp_perc`) follows it, where the subject's HSPs are given to the
/// collector (`ungapped_reap_query_coverage`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2609-2737
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2694-2697
/// The HSP is read in the concatenated query (`query_concat`) from NCBI's start, `query_start` plus
/// `hsp->query.offset`, which lies outside the context for an HSP whose ungapped hit passes an end
/// of it (a negative `PrelimHit::prelim_qs`, or a `prelim_qe` past the context's length); the
/// context's letters are the same bytes. Where the re-evaluation keeps a piece outside the
/// context, NCBI's HSP keeps those offsets, and so does LOSAT's through the later stages, which may
/// drop it; a report that reaches it is refused (`post_process_hits_and_write`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2708-2716
/// The identities are counted in the concatenated letters without masks (`query_concat_nomask`,
/// NCBI's `sequence_nomask`) from the same start (`&query[q_off]`, c++/src/algo/blast/core/blast_hits.c:764).
#[allow(clippy::too_many_arguments)]
fn ungapped_subject_stage(
    hits: &mut Vec<PrelimHit>,
    s_blastna: &[u8],
    query_concat: &[u8],
    query_concat_nomask: &[u8],
    score_matrix: &[i32; 256],
    word_cutoffs: &[i32],
    percent_identity: f64,
    evalue_of: impl Fn(&PrelimHit) -> f64,
    reap: PrelimReap,
) {
    hits.retain_mut(|hit| {
        let context_idx = hit.context_idx as usize;
        let context_start = hit.query_context_offset as isize;
        let Some(reevaluated) = ungapped::reevaluate_with_ambiguities_ungapped(
            query_concat,
            s_blastna,
            (context_start + hit.prelim_qs as isize) as usize,
            hit.prelim_ss,
            (hit.prelim_qe - hit.prelim_qs) as usize,
            word_cutoffs[context_idx],
            score_matrix,
        ) else {
            return false;
        };
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:453-456
        // A kept piece before the context (the part there reaches the cutoff and the context's
        // letters score no more; the sentinel before the context starts the sum over,
        // blast_stat.c:1122-1127) has negative offsets (remaining-gaps X1b-o, 2026-11-17).
        hit.prelim_score = reevaluated.score;
        hit.prelim_qs = (reevaluated.q_offset as isize - context_start) as i32;
        hit.prelim_qe = (reevaluated.q_end as isize - context_start) as i32;
        hit.prelim_ss = reevaluated.s_offset;
        hit.prelim_se = reevaluated.s_end;
        let align_length = (hit.prelim_qe - hit.prelim_qs) as usize;
        let identities = ungapped::count_identities_ungapped(
            query_concat_nomask,
            s_blastna,
            reevaluated.q_offset,
            hit.prelim_ss,
            align_length,
        );
        // `Blast_HSPTest` (cited above); NCBI's `min_hit_length` is always 0
        // (`core::blast_hits::hsp_test`).
        !crate::core::blast_hits::hsp_test(
            identities as i64,
            align_length as i64,
            percent_identity,
            0,
        )
    });
    sort_prelim_hits_by_score(hits);
    reap.reap(hits, evalue_of);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1542-1547
// NCBI reference (598d8ae6): c++/src/objects/seq/seq_align_mapper_base.cpp:1288-1294, 2901-2904
/// LOSAT's refusal (`hsp_left_of_query`) of the first query of a batch whose final hit list keeps
/// an HSP that NCBI keeps outside its query context (`BlastnHsp::internal_q_offset_0`) with a
/// Seq-align query start below 0 in the query's record (`BlastnHsp::starts_left_of_record`: a
/// plus-strand HSP left of its context, a minus-strand one right of it, that the `-query_loc`
/// interval's start does not move to 0 or past it); `hit_lists`, `lengths`, `offsets` (the
/// intervals' starts, `Placements::offset`) and `records` (the searched records) are the
/// batch's, by query of the batch. A list that holds such an HSP is copied and made final as the
/// reports get it: a subject set's as `post_process_hits_and_write` finishes it, a database's
/// (finished after the traceback) without its HSPs of score 0; no list is copied in a search
/// without one (an ungapped X-drop of `INT_MIN`, `-xdrop_ungap 1e10`, gives them, remaining-gaps
/// X1b-o).
///
/// NCBI moves the query row of each Seq-align by the interval's start (`RemapToQueryLoc`,
/// blast_seqalign.cpp:1634 and 1641, before the subject's mapping): an ungapped search's
/// Dense-diag whatever its start, a gapped search's Dense-seg only at its starts at 0 or past it.
/// Of a subject set NCBI then maps each Seq-align of the batch to the subject's location in the
/// search (`s_RemapToSubjectLoc`, blast_seqalign.cpp:1630-1636): the query row of such an HSP maps
/// to nothing, which throws `CAnnotMapperException` (exit status 255) at
/// seq_align_mapper_base.cpp:1292 or 2902 (both seen with `-ungapped`, by task; a gapped HSP's
/// Dense-seg not reached).
/// Its start at 0 or past it is shown (`algorithm::blastn::pairwise::context_letters`;
/// remaining-gaps X1b-p).
fn left_of_query_refusal(
    args: &BlastnArgs,
    source: SubjectSource<'_>,
    hit_lists: &[Option<BlastnHitList>],
    lengths: &[usize],
    offsets: &[usize],
    records: &[FastaRecord],
) -> Option<anyhow::Error> {
    let database_search = matches!(source, SubjectSource::Database(_));
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1609-1628
    // An ungapped search's HSPs are Dense-diags (`x_UngappedHSPToDenseDiag`).
    let dense_diag = args.ungapped;
    let keeps_one = |hit_list: &BlastnHitList, offset: usize| {
        hit_list
            .hsplist_array
            .iter()
            .flat_map(|hsp_list| &hsp_list.hsps)
            .any(|hsp| hsp.starts_left_of_record(offset, dense_diag))
    };
    let q_idx = hit_lists.iter().enumerate().position(|(q_idx, hit_list)| {
        let offset = offsets.get(q_idx).copied().unwrap_or(0);
        let Some(hit_list) = hit_list
            .as_ref()
            .filter(|hit_list| keeps_one(hit_list, offset))
        else {
            return false;
        };
        let mut finished = hit_list.clone();
        if !database_search {
            finish_query_hit_list(
                &mut finished,
                lengths[q_idx],
                args.report_counts.hitlist_size,
                args.max_hsps_per_subject.unwrap_or(0),
                args.qcov_hsp_perc.unwrap_or(0.0),
                args.subject_besthit,
                args.culling_limit_value(),
                args.best_hit(),
                args.ungapped,
                true,
            );
        }
        drop_hsps_of_score_zero(&mut finished);
        keeps_one(&finished, offset)
    })?;
    let id = records.get(q_idx).map_or(&[][..], FastaRecord::shown_id);
    Some(hsp_left_of_query(id, database_search))
}

/// LOSAT's refusal of the query `query_id`, whose final hit list keeps an HSP with a Seq-align
/// query start below 0 (`left_of_query_refusal`). Of a subject set NCBI stops with an exception
/// there, after what it wrote before (`fail_batch!` in `run_in_pool`). Of a database
/// (`database_search`) NCBI maps nothing, and its reports read the query row of the Seq-align
/// each its own way: some stop with an exception after the rows before, others leave the HSP out
/// or show it; LOSAT does not reproduce them yet (remaining-gaps X1b-p, decision X1b-p-1).
#[cold]
#[inline(never)]
fn hsp_left_of_query(query_id: &[u8], database_search: bool) -> anyhow::Error {
    let ncbi = if database_search {
        "whose reports from a database NCBI BLAST+ writes each its own way (some stop with an exception), which LOSAT does not reproduce yet"
    } else {
        "where NCBI BLAST+ stops with an exception when it maps the alignment to the subject (exit status 255); LOSAT refuses it after what NCBI BLAST+ writes before"
    };
    anyhow::anyhow!(
        "BLASTN: query {}: an HSP that NCBI BLAST+ keeps outside its query context (from an ungapped extension past an end of the context, such as with -xdrop_ungap 1e10) starts before the query's first letter, {ncbi}",
        String::from_utf8_lossy(query_id)
    )
}

/// The last step of the preliminary stage's finish of a subject's HSP list in an ungapped
/// (`-ungapped`) search: the HSPs below `-qcov_hsp_perc` percent of their query are deleted
/// (`Blast_HSPListReapByQueryCoverage`). It does not depend on the subject, so the collector's
/// caller makes it, where the subject's HSPs are given over.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1482-1541
fn ungapped_reap_query_coverage(
    hits: &mut Vec<PrelimHit>,
    query_cov_hsp_perc: f64,
    context_length: impl Fn(usize) -> usize,
) {
    reap_by_query_coverage(hits, query_cov_hsp_perc, |hit| {
        (
            hit.prelim_qs as i64,
            hit.prelim_qe as i64,
            context_length(hit.context_idx as usize) as i64,
        )
    });
}

/// LOSAT's refusal of a word cutoff whose arithmetic NCBI does in a way that overflows (memory
/// "no NCBI UB / overflow"): a context of 2^30 letters or more, or a search space beyond `Int8`.
/// The word cutoff of an ungapped search and of matrix-only scoring (`batch_setup::batch_cutoffs`).
fn ungapped_cutoff_limit_error(limit: ungapped::CutoffLimit) -> anyhow::Error {
    let what = match limit {
        ungapped::CutoffLimit::QueryLength => {
            "a query of 2^30 letters or more, whose doubled length overflows NCBI BLAST+'s 32-bit int in the ungapped word cutoff"
        }
        ungapped::CutoffLimit::SearchSpace => {
            "an average subject length whose ungapped word cutoff search space overflows NCBI BLAST+'s 64-bit int"
        }
    };
    anyhow::anyhow!(
        "{what} is not supported by LOSAT's BLASTN with -ungapped or with reward 0 and penalty 0 (matrix-only scoring)"
    )
}

/// What the query batches search (`search_query_batch`): the prepared subject set, or the
/// database, which is read in blocks for each batch (`database::search_blocks`).
#[derive(Clone, Copy)]
enum SubjectSource<'a> {
    Prepared(&'a PreparedSubjects<'a>),
    Database(&'a database::BlastnDatabase),
}

impl<'a> SubjectSource<'a> {
    /// The number of sequences and the total length (the real ones, also with `-dbsize`).
    fn metadata(self) -> &'a SubjectMetadata {
        match self {
            SubjectSource::Prepared(prepared) => &prepared.metadata,
            SubjectSource::Database(database) => database.metadata(),
        }
    }

    // No NCBI counterpart: the schedules (`subjects_parallel`, the progress bar); it changes no
    // value NCBI computes.
    /// The number of subjects the search reads: the records, or the database's included
    /// sequences (`SeqDb::num_seqs`; the statistics' number, `metadata`, can be another).
    pub(super) fn num_subjects(self) -> usize {
        match self {
            SubjectSource::Prepared(prepared) => prepared.metadata.db_num_seqs as usize,
            SubjectSource::Database(database) => database.num_subjects(),
        }
    }

    /// The prepared subject set; none for a database.
    fn prepared(self) -> Option<&'a PreparedSubjects<'a>> {
        match self {
            SubjectSource::Prepared(prepared) => Some(prepared),
            SubjectSource::Database(_) => None,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:129-153
/// What the reports show of the subjects (NCBI's `m_DbInfo` and `m_DbName`): the name of
/// outfmt 7's `# Database:` line and of the outfmt 0 `Database:` lines of a subject set
/// (`User specified sequence set (Input: PATH)`; the `-db` value for a database), the
/// counts, and for a database the lines of `GetBlastDbInfo` (`PairwiseDatabase`).
struct ShownSubjects {
    name: String,
    /// The statistics' number of sequences, NCBI's `Int4` (`SubjectMetadata::db_num_seqs`).
    num_seqs: i32,
    /// The statistics' total length, NCBI's `Int8` (`SubjectMetadata::db_len_total`).
    total_letters: i64,
    database: Option<PairwiseDatabase>,
    /// Whether `taxdb.bti` is on NCBI's search path, for the custom tabular reports
    /// (`BlastnArgs::tabular_taxdb`; not found where it was not asked).
    tabular_taxdb_found: bool,
    /// The run's `-html` data for the outfmt 0 to 4 prolog (`BlastnArgs::html_links`,
    /// `write_blastn_pairwise_prolog`); `None`: the plain prolog.
    html_links: Option<std::sync::Arc<crate::report::html_ids::HtmlLinks>>,
}

impl ShownSubjects {
    fn new(args: &BlastnArgs, source: SubjectSource<'_>) -> Self {
        let metadata = source.metadata();
        match source {
            SubjectSource::Prepared(_) => ShownSubjects {
                name: format!(
                    "User specified sequence set (Input: {})",
                    args.subject_path().display()
                ),
                num_seqs: metadata.db_num_seqs,
                total_letters: metadata.db_len_total,
                database: None,
                tabular_taxdb_found: args.tabular_taxdb.unwrap_or(false),
                html_links: args.html_links.clone(),
            },
            SubjectSource::Database(database) => ShownSubjects {
                name: database.name().to_string(),
                num_seqs: metadata.db_num_seqs,
                total_letters: metadata.db_len_total,
                database: Some(database.pairwise_database_masked(args.db_mask_resolved.as_ref())),
                tabular_taxdb_found: args.tabular_taxdb.unwrap_or(false),
                html_links: args.html_links.clone(),
            },
        }
    }
}

/// One subject as the per-subject kernel (`SubjectPasses::search_subject`)
/// reads it: a record of a subject set, or a sequence of a database.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
#[derive(Clone, Copy)]
struct BlastnSubjectRef<'a> {
    /// The number of bases.
    length: usize,
    /// The packed ncbi2na letters that the preliminary stage scans (four bases per byte,
    /// the count of the last byte's bases in its two low bits): the subject set's
    /// (`encode_subject_ncbi2na_packed`) or the database's stored data.
    packed: &'a [u8],
    /// The letters of a subject set's record: its lowercase masks (`-lcase_masking`) and
    /// the traceback's BLASTNA letters come from them. `None` for a database sequence.
    letters: Option<&'a [u8]>,
    /// The BLASTNA letters of a database sequence in the traceback
    /// (`database::traceback_listed`), without sentinel bytes; `None` otherwise.
    blastna: Option<&'a [u8]>,
    /// EXPERIMENT (LOSAT_X_BNSTORE): the packed subject of a subject set's record that keeps
    /// no letters: its lowercase masks and the traceback's BLASTNA letters come from it.
    store: Option<&'a PackedNucleotide>,
    /// EXPERIMENT (LOSAT_X_BNTBRANGE): with `blastna` decoded in ranges only, set when a
    /// traceback alignment met a fence (`alignment::x_fence`): the traceback stops, and
    /// `database::traceback_listed` traces the sequence again whole. `None` otherwise.
    x_fenced: Option<&'a std::sync::atomic::AtomicBool>,
    /// The block of the preliminary stage that a database sequence was read in
    /// (`database::block_subject`), whose masks a search with `-db_soft_mask` reads once for all
    /// its sequences (`database::BlastnDatabase::seq_ranges_into`). `None` otherwise.
    block: Option<&'a crate::seqdb::BlockRange>,
}

impl<'a> BlastnSubjectRef<'a> {
    /// The record `record` of a subject set, with its packed letters `packed`.
    fn record(record: &'a FastaRecord, packed: &'a [u8]) -> Self {
        BlastnSubjectRef {
            length: record.seq().len(),
            packed,
            letters: Some(record.seq()),
            blastna: None,
            store: None,
            x_fenced: None,
            block: None,
        }
    }

    /// EXPERIMENT (LOSAT_X_BNSTORE): a record of a subject set held as the packed subject
    /// `store`; the preliminary stage scans its own packed bytes.
    fn packed_record(store: &'a PackedNucleotide) -> Self {
        BlastnSubjectRef {
            length: store.len(),
            packed: store.packed(),
            letters: None,
            blastna: None,
            store: Some(store),
            x_fenced: None,
            block: None,
        }
    }
}

/// What both stages of `search_subject` read of one call: the subject's index, trace id and
/// length. Built once per call, after its early return without a valid context, for
/// `prelim_subject::prelim_subject` or `traceback_subject::traceback_subject`; the batch's values
/// (the queries, their contexts, the effective search spaces, the cutoffs, the LOSAT_TIMING
/// counters) are read from the stage's batch (`PrelimBatch`, `TraceBatch`), built once per batch,
/// or with the statistics of each subject from `SubjectPasses::special_subject`'s copy of it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1412-1413,1446
/// (`BLAST_PreliminarySearchEngine`'s loop on the subjects)
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1521
/// (`BLAST_ComputeTraceback_MT`'s loop on the subjects' batches)
/// The trace id has no NCBI counterpart.
struct SubjectCtx<'a> {
    s_idx: usize,
    s_id: &'a str,
    s_len_full: usize,
}

impl SubjectCtx<'_> {
    /// LOSAT_DEBUG_BLASTN: the subject's cutoffs (the batch's `cutoff_scores`), printed by each
    /// stage of `search_subject` before it searches the subject (`prelim_subject::prelim_subject`,
    /// `traceback_subject::traceback_subject`). No NCBI counterpart (a LOSAT diagnostic line).
    fn debug_cutoffs(&self, cutoff_scores: &[i32]) {
        eprintln!(
            "[BLASTN_DEBUG] Subject {} (len={}): cutoff_scores={:?}",
            self.s_id, self.s_len_full, cutoff_scores
        );
    }
}

/// The results of one query batch, with the values of its query contexts that the reports
/// use.
struct QueryBatch {
    hit_lists: Vec<Option<BlastnHitList>>,
    good_init_extends: u64,
    searched: bool,
    query_masks: Vec<Vec<MaskedInterval>>,
    context_karlin: Vec<Option<ContextKarlin>>,
    query_eff_searchsp: Vec<i64>,
    /// The preliminary hit lists of a query chunk, one per query part
    /// (`BatchStage::ChunkPrelim`).
    chunk_prelim_lists: PrelimHitLists,
    /// The batch's `LOSAT_TIMING` lines (`blastn_timing_text`; empty without it), which
    /// `run_in_pool` writes to standard error when it takes the batch, where the search wrote
    /// them before.
    timing: String,
    /// The first of the batch's fetches with `check_oid_exclusion` (the traceback's, or the
    /// ungapped re-evaluation's) that read a header with a UTF8String where a VisibleString is
    /// expected: NCBI's warning comes there, before its Seq-aligns' reads (`read_titles`).
    fetched_utf8: Option<crate::seqdb::SeqLocation>,
}

/// Which of NCBI's stages `search_query_batch` runs for its queries.
#[derive(Clone, Copy)]
enum BatchStage<'a> {
    /// The preliminary search and the traceback of a query batch (`CLocalBlast::Run`). A
    /// batch that NCBI splits has its preliminary search run in query chunks
    /// (`search_query_chunks`), which need where each query's searched letters start in
    /// its record (`-query_loc`).
    ///
    /// `query_masks` are the masks of the batch's queries when the loop has made them
    /// (`-window_masker_db`, whose errors stop the batch before its search); `None`: made here.
    /// `lower_case_before`: whether the letter before each query's interval is lower case, which
    /// the chunks restrict too (empty without a range, without `-lcase_masking` and with
    /// `-window_masker_db`). `masks_outside`: with `-filtering_db`, the masks next to each query's
    /// interval that the batch's masks joined (`repeat_filter::MasksOutside`; empty without them).
    ///
    /// `last`: the batch loop reads no batch after this one (`batches_follow`), so the loop's
    /// keep of LOSAT_X_MBLUTKEEP frees its table after the batch's preliminary stage
    /// (`x_lut_drop::after_prelim`), and its `qs` (LOSAT_X_CSRQSKEEP) once the batch's tables are
    /// built (`BatchStage::x_last_tables`).
    ///
    /// `setup_only`: the set-up of the batch alone, up to its checks of the subjects
    /// (`context_statistics`), which ends with `SetUpOnly` where they pass: the set-up of a batch
    /// that LOSAT refuses (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20).
    Search {
        query_offsets: &'a [usize],
        query_masks: Option<&'a [Vec<MaskedInterval>]>,
        lower_case_before: &'a [bool],
        masks_outside: &'a [MasksOutside],
        last: bool,
        setup_only: bool,
    },
    /// The preliminary search of one query chunk of a split batch
    /// (`SplitQuery_CreateChunkData`), with the masks of its query parts and the search
    /// spaces of the batch's contexts.
    ///
    /// EXPERIMENT (LOSAT_X_CSRQSKEEP) `x_last`: the last chunk of the batch loop's last batch, after
    /// whose tables no build follows (`BatchStage::x_last_tables`), so the loop's keep frees the
    /// `qs` it holds once they are built and takes none back when they drop
    /// (`x_mb_lut_keep::Keep::close_qs`).
    ChunkPrelim {
        query_masks: &'a [Vec<MaskedInterval>],
        batch_eff_searchsp: &'a [i64],
        x_last: bool,
    },
}

impl BatchStage<'_> {
    /// Whether the batch loop builds no lookup table after this call's (EXPERIMENT
    /// LOSAT_X_MBLUTKEEP, LOSAT_X_CSRQSKEEP; no NCBI counterpart, memory only): the loop's last
    /// batch (`Search`'s `last`; a split one builds none of its own, after its chunks') or the last
    /// chunk of a split last batch (`ChunkPrelim`'s `x_last`). The one test of it for the loop's
    /// keep: it closes the keep to `qs` once these tables are built and frees its table after the
    /// last batch's preliminary stage.
    fn x_last_tables(self) -> bool {
        matches!(
            self,
            BatchStage::Search { last: true, .. } | BatchStage::ChunkPrelim { x_last: true, .. }
        )
    }
}

/// The end of a batch's set-up alone (`BatchStage::Search::setup_only`) whose checks pass: no
/// error of NCBI's. Never shown: the refusal of the batch takes its place.
#[derive(Debug)]
struct SetUpOnly;

impl std::fmt::Display for SetUpOnly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("the set-up of the batch passes")
    }
}

impl std::error::Error for SetUpOnly {}

/// Whether `run_in_pool`'s batch loop reads a batch from input record `input_start` of
/// `records`: records are left, or what ends the reading after them (`reading`: blank lines or
/// an error, `QueryReading::reads_past_records`) has not come yet (`reading_ended`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277
fn batches_follow(
    input_start: usize,
    records: usize,
    reading: &QueryReading,
    reading_ended: bool,
) -> bool {
    input_start < records || (reading.reads_past_records() && !reading_ended)
}

/// Searches the queries in NCBI's query batches and writes the results of all of them.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:261-300
/// `GetQueryBatchSize` is the `BATCH_SIZE` value (0 when it is not set) and the chunk size
/// that of `CHUNK_SIZE` or the task (`query_batching_from_environment`).
#[allow(clippy::too_many_arguments)]
fn run_in_pool(
    args: BlastnArgs,
    queries: x_query_stream::Queries<'_>,
    reading: QueryReading,
    query_input: &QueryInput,
    subjects: Subjects<'_>,
    subject_ranges: SubjectRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: Vec<BlastnOutputFormat>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: Result<(QueryBatching, Result<()>)>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:258-261
    // What the run writes before its query batches ends here.
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.queries_begin();
    }
    let prepared: PreparedSubjects<'_>;
    let source = match subjects {
        Subjects::Records(subject_set) => {
            let subject_records = subject_set.records();
            let metadata =
                subject_metadata_from_records(subject_records, |index| subject_set.len(index));
            if metadata.db_num_seqs == 0 {
                // EXPERIMENT (LOSAT_X_QSTREAM): the records are read again for their warnings.
                queries.write_input_warnings(outputs.diagnostics)?;
                return Ok(());
            }
            // Diagnostic (LOSAT_TIMING): the subjects' packing for the common line
            // (utils/xtiming.rs). No NCBI counterpart.
            let t_encode = crate::utils::xtiming::clock();
            // The packed subjects hold the packed letters (LOSAT_X_BNSTORE).
            let packed = if subject_set.stores().is_some() {
                Vec::new()
            } else {
                subject_records
                    .iter()
                    .map(|record| encode_subject_ncbi2na_packed(record.seq()))
                    .collect()
            };
            crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
            prepared = PreparedSubjects {
                records: subject_records,
                stores: subject_set.stores(),
                packed,
                metadata,
            };
            SubjectSource::Prepared(&prepared)
        }
        // A database without sequences fails in the first searched batch
        // (`search_query_batch`, in `context_statistics`).
        Subjects::Database(database) => SubjectSource::Database(database),
    };
    // The subjects as the reports show them, and the stored titles of the database
    // sequences that the final hit lists of the batches name (`database::read_titles`).
    let shown = ShownSubjects::new(&args, source);
    let mut titles = database::DatabaseTitles::new();
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`), added as the
    // loop cuts the batches.
    // With outfmt 5 against a database, every query's ids from the start (`with_whole_run_of`):
    // outfmt 5 reads the queries' letters, which keeps the records read before the search
    // (`x_query_drop::Rule`), and `-parse_deflines` reads no query as a stream.
    let subject_records = matches!(source, SubjectSource::Prepared(_));
    let xml = output_formats.contains(&BlastnOutputFormat::Xml);
    let mut query_scopes = crate::blastinput::query_batch::QueryScopes::of_search(
        args.parse_deflines,
        subject_records,
        xml,
        "BLASTN",
    )
    .map(|scopes| match &queries {
        x_query_stream::Queries::Records { input, .. } if xml && !subject_records => {
            scopes.with_whole_run_of(input)
        }
        _ => scopes,
    })
    .map(|scopes| scopes.with_refusal(args.scope_refusal.clone()));
    let intro = blastn_intro(&args);
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:258-261
    // NCBI writes the outfmt 0 prolog before it reads the first query batch, and its stream
    // throws at the prolog's first flush when the write fails (blast_format.cpp:118-119),
    // before any query is read or searched.
    write_pairwise_prologs(outputs, &output_formats, intro, &shown, true)?;
    // The mixer's `CHUNK_SIZE` and then `BATCH_SIZE` are read here, after the prolog and before
    // the first batch: a value that NCBI cannot convert stops the run with its exception
    // (`query_batching_from_environment`), as a batch size of 0 does below.
    let (mut batching, mut overlap_error) = match batching {
        Ok(batching) => batching,
        Err(error) => {
            write_pairwise_prologs(outputs, &output_formats, intro, &shown, false)?;
            return Err(error);
        }
    };
    let mut mixer = BatchSizeMixer::new(batching.split.mixer_max_batch_size());
    if batching.fixed_batch_size == 0 {
        // A database's total is the `Int8` of its index files (`GetDbTotalLength` sums
        // `m_DbInfo`, `GetBlastDbInfo`: the whole database's, also with a list).
        let total_length = shown
            .database
            .as_ref()
            .map_or(source.metadata().db_len_total, |database| {
                database.total_letters
            });
        if total_length > 0 {
            mixer.set_target_hits((total_length / 3000) as i32);
        }
        // `SetBatchSize` of the input, given `GetBatchSize` of the mixer, converts the `Int4` to `TSeqPos`.
        batching.batch_size = mixer.batch_size(None) as u32;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-280
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // A reader at the end of its input reads no batch (its `End()` is true): the report of no
    // query, the prolog and the epilog (an empty pipe, or `-query -` after `-subject -`).
    if queries.input_count() == 0 && matches!(reading, QueryReading::End) {
        return write_batch_reports(
            &args,
            x_query_stream::ReportQueries::default(),
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &ReportRanges {
                queries: query_input,
                subjects: subject_ranges,
            },
            outputs,
            &output_formats,
            BatchResults {
                hit_lists: Vec::new(),
                warnings: &[],
                reads: &[],
                lengths: &[],
                unsearched: Vec::new(),
                query_masks: &[],
                repeat_hits: &[],
                context_karlin: &[],
                query_eff_searchsp: &[],
                batches: &[],
            },
            true,
            None,
        );
    }
    // A batch size of 0 (`CHUNK_SIZE` 1000 without `BATCH_SIZE`: the mixer's maximum is
    // 0) reads no query, so the first batch fails after the prolog (outfmt 0).
    if batching.batch_size == 0 {
        write_pairwise_prologs(outputs, &output_formats, intro, &shown, false)?;
        return Err(crate::cli::NativeError {
            exit: 3,
            message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
        }
        .into());
    }

    // The lengths of the searched queries (their intervals with `-query_loc`; with
    // LOSAT_X_QSTREAM from the first reading's lengths).
    let lengths: Vec<usize> = queries.searched_lengths();
    let report_ranges = ReportRanges {
        queries: query_input,
        subjects: subject_ranges,
    };
    let query_offsets: Vec<usize> = (0..lengths.len())
        .map(|index| query_input.placements.offset(index))
        .collect();
    // EXPERIMENT (LOSAT_X_RPTSTREAM): whether the loop writes each batch's reports when it takes
    // the batch (`x_rpt_stream.rs`); the summary lines come when it is dropped, at the end.
    let mut stream = x_rpt_stream::Stream::start(
        &args,
        outputs,
        &output_formats,
        source,
        subject_ranges,
        query_input,
        &lengths,
        batching,
    );
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): the megablast table this loop keeps between its batches
    // (`lookup::x_mb_lut_keep`), whose table is freed after the last batch's preliminary stage.
    let x_lut_keep = x_mb_lut_keep::Keep::for_loop();
    // No NCBI counterpart (memory placement only; NCBI allocates the buffer of each batch,
    // `x_query_sentinels`): the query buffer this loop keeps between its batches, freed when the
    // loop ends.
    let mut query_keep: Vec<u8> = Vec::new();
    let mut hit_lists: Vec<Option<BlastnHitList>> = Vec::with_capacity(lengths.len());
    let mut batch_spans: Vec<std::ops::Range<usize>> = Vec::new();
    let mut unsearched = Vec::with_capacity(lengths.len());
    let mut query_masks = Vec::with_capacity(lengths.len());
    // With `-filtering_db` and the BLAST archive, the repeat hits of each query that its batch's
    // set-up found (`repeat_batch_masks`), which the archive's masks join
    // (`BatchResults::repeat_hits`).
    let keeps_repeat_hits =
        args.repeat_filter.is_some() && output_formats.contains(&BlastnOutputFormat::Archive);
    let mut repeat_hits: Vec<Vec<(usize, usize)>> = Vec::new();
    let mut context_karlin = Vec::with_capacity(2 * lengths.len());
    // The warnings written before each query's report (`QueryWarnings`).
    let mut warnings: Vec<Posts> = vec![Posts::new(); lengths.len()];
    // The texts of the subjects that each query's report reads one by one, by query, where the
    // first format places them (`BlastnReportInputs::query_reads`).
    let mut query_reads: Vec<QueryReads> = vec![Vec::new(); lengths.len()];
    // The handle of the BLAST archive, which saves each whole post of `warnings` as one message
    // (`Posts::write`).
    let archive = args.archive.clone();
    let archive = archive.as_ref();
    // The search strategy's requests of the batches (`-export_search_strategy`), with every
    // `-subject` record.
    let requests = args.strategy.as_ref().map(|_| {
        archive::requests(
            match source {
                SubjectSource::Prepared(prepared) => {
                    Some(subject_ranges.records.unwrap_or(prepared.records))
                }
                SubjectSource::Database(_) => None,
            },
            args.lcase_masking,
        )
    });
    let mut query_eff_searchsp = Vec::with_capacity(2 * lengths.len());
    // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQPACK): the searched records of the batches taken,
    // kept where a writer reads them after their batch: the reports after the loop and the
    // error paths that write the reports of the batches before (LOSAT_X_RPTSTREAM not acting).
    // Without the switches the records read before the search are those.
    let keep = queries.per_batch() && !stream.acting();
    let mut kept = x_query_stream::Kept::default();
    // A batch is read from the input records (`input_start..input_end`); its searched
    // queries are `start..end` (all of them without `-query_loc`). A batch size counts the
    // whole records, and a record whose range starts past its end is skipped
    // (`next_ranged_batch`).
    let mut reading = reading;
    // Whether what ended the reading after the records (`reading`) has come.
    let mut reading_ended = false;
    let mut input_start = 0;
    // NCBI's reader's count of the records into the batches (`BatchCount`).
    let mut batch_count = BatchCount::of(query_scopes.as_ref());
    // Where a query's report of a batch whose prefetch is off ends the run at a subject whose
    // load throws, no report and no batch's prefetch after it loads anything
    // (`database::QueryFetches`).
    let run_loads = crate::report::db_subject::RunLoads::new(database::after_throw(
        &output_formats,
        outputs.other_outfmts,
        args.html_links.is_some(),
    ));
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // The reader is at the end of its input after the last record, unless an `eEOF` or an
    // error ends the reading there; that comes in the batch being read when the records
    // run out before the batch reaches its size, and otherwise in a batch of its own.
    while batches_follow(
        input_start,
        query_input.input_lengths.len(),
        &reading,
        reading_ended,
    ) {
        let (input_end, reached_size) = BatchCount::next_batch(
            batch_count.as_mut(),
            &|record| queries.believed(record),
            &query_input.input_lengths,
            &query_input.skipped,
            input_start,
            batching.batch_size,
        );
        let std::ops::Range { start, end } = query_input.searched_in(input_start..input_end);
        let reading_ends_here = input_end == query_input.input_lengths.len()
            && !reached_size
            && reading.reads_past_records();
        // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
        // `CFastaReader` writes its messages about a batch's queries (their lines and
        // titles) when it reads them, after the report of the batch before: before the
        // report of the batch's first query (`QueryWarnings`). Skipped records are read too.
        // EXPERIMENT (LOSAT_X_QSTREAM): the batch's records, read now from the query file when
        // the switch acts (an error here is the file changed between the two readings).
        let mut batch_records =
            queries.batch_in(input_start..input_end, start..end, parallel_pool)?;
        let mut batch_title_warnings = Posts::of_lines(
            batch_records
                .input()
                .iter()
                .flat_map(|record| record.warnings.iter().copied())
                .collect(),
        );
        // The reports of the batches before a batch at which the run stops (for the first batch,
        // the prologs that come with them), before the batch's messages and its error.
        macro_rules! write_reports_before {
            () => {{
                // EXPERIMENT (LOSAT_X_CSRQSKEEP): the loop builds no table after this, so its keep
                // frees the `qs` it holds before the reports are written (`Keep::close_qs`;
                // nothing with the switch off). No NCBI counterpart (memory only).
                if let Some(keep) = &x_lut_keep {
                    keep.close_qs();
                }
                // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                // written when the loop took them.
                let before = |outputs: &mut ReportOutputs<'_>| {
                    write_before_failed_batch(
                        &args,
                        queries.before(&kept, start),
                        ReportSubjects {
                            source,
                            shown: &shown,
                            titles: &titles,
                        },
                        &report_ranges,
                        outputs,
                        &output_formats,
                        BatchResults {
                            hit_lists,
                            warnings: &warnings[..start],
                            reads: &query_reads[..start],
                            lengths: &lengths[..start],
                            unsearched,
                            query_masks: &query_masks,
                            repeat_hits: &repeat_hits,
                            context_karlin: &context_karlin,
                            query_eff_searchsp: &query_eff_searchsp,
                            batches: &batch_spans,
                        },
                    )
                };
                // In the first batch nothing was streamed: the call of the off path writes the
                // prologs that come with the reports (outfmt 20's line of the field names).
                if stream.acting() && start > 0 {
                    stream.finish_sam(outputs)?;
                } else {
                    before(outputs)?;
                }
            }};
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input.cpp:145-152
        // An `eEOF` ends the batch with the records read so far. A reader error (or LOSAT's
        // rejection of a Seq-id line) stops the run while the batch is read: after the
        // reports of the batches before (without the epilog) and the messages of the
        // batch's records, which are neither searched nor reported.
        if reading_ends_here {
            reading_ended = true;
            if let QueryReading::Error {
                error,
                warnings: error_warnings,
            } = std::mem::replace(&mut reading, QueryReading::End)
            {
                write_reports_before!();
                batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
                outputs.diagnostics.write_all(&error_warnings)?;
                return Err(error.into_app_error());
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
        // A batch of skipped records only (with `-query_loc`), or a batch that an `eEOF`
        // ends before its first record, fails after the reports of the batches before,
        // without the epilog.
        if start == end {
            write_reports_before!();
            batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
            return Err(crate::cli::NativeError {
                exit: 3,
                message: format!("{ENGINE_ERROR} Empty CBlastQueryVector\n").into(),
            }
            .into());
        }
        // LOSAT's refusal of two queries of the batch's scope with one id, in the one form of every
        // program (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20): with
        // `-export_search_strategy` before the batch's request, which NCBI writes with the earlier
        // record's letters for the repeated id (no request of the batch is written); else after the
        // batch's checks before its set-up, where the set-up of the batch's records that are first
        // in the scope gives NCBI's error of the subjects or the refusal comes, before the search.
        // (Round 4 refused after the search, and wrote the batch's request.)
        let mut scope_refusal = None;
        if let Some(scopes) = query_scopes.as_mut() {
            if let Err(refusal) = scopes.add_batch(
                batch_records.input(),
                &query_input.skipped[input_start..input_end],
                start,
            ) {
                scope_refusal = Some(refusal);
            }
        }
        if let Some(refusal) = scope_refusal.take_if(|_| args.strategy.is_some()) {
            write_reports_before!();
            batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
            return Err(refusal);
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:280-282
        // The batch's search strategy request (`-export_search_strategy`): once its query factory
        // is made, before its set-up (the masks of `-window_masker_db` first) and search, so every
        // error of the batch from here on follows it (the file is apart from the reports). Its
        // `EffectiveSearchSpace` option is the one the batch before left (`note_searched` below).
        // NCBI's posts of the bytes it replaces follow the reader's messages of the batch, before
        // anything of its set-up, search and reports (`batch_title_warnings`).
        if let (Some(strategy), Some(requests)) = (args.strategy.as_ref(), requests.as_ref()) {
            requests.write(
                strategy,
                start,
                batch_records.searched(),
                query_input.placements.get(start),
            )?;
            batch_title_warnings.append(&strategy.take_all_posts());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.hpp:104-106
        // The query splitter of the batch reads `CHUNK_SIZE` (converted already by the mixer) and
        // `OVERLAP_CHUNK_SIZE` before the batch's queries are set up (their masks below): a value
        // that NCBI cannot convert stops the first batch here, after its reader's messages
        // (`query_batching_from_environment`).
        if input_start == 0 {
            if let Err(error) = std::mem::replace(&mut overlap_error, Ok(())) {
                write_before_failed_batch(
                    &args,
                    queries.before(&kept, start),
                    ReportSubjects {
                        source,
                        shown: &shown,
                        titles: &titles,
                    },
                    &report_ranges,
                    outputs,
                    &output_formats,
                    BatchResults {
                        hit_lists,
                        warnings: &warnings[..start],
                        reads: &query_reads[..start],
                        lengths: &lengths[..start],
                        unsearched,
                        query_masks: &query_masks,
                        repeat_hits: &repeat_hits,
                        context_karlin: &context_karlin,
                        query_eff_searchsp: &query_eff_searchsp,
                        batches: &batch_spans,
                    },
                )?;
                batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
                return Err(error);
            }
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:518-530
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:149-152
        // With `-window_masker_db` the masks come first in each query's set-up: a counts
        // file NCBI cannot read, or a bad index of its hash table, fails every query of the
        // batch (`window_masker::batch_error`), also those without letters. The preliminary
        // search sets up the whole batch (`BlastSetupPreliminarySearchEx`, called by
        // `CBlastPrelimSearch::x_Init`) before it splits the batch into chunks, so a split
        // batch fails the same way; its chunks then make their own masks in the search, and
        // an error there (the whole queries' masks were made) is LOSAT's rejection
        // (`window_masker::rejection`). The masks made here are those the search uses.
        let mut window_masker_masks = None;
        // The masks' failure of the batch as NCBI reports it: the messages of its queries' set-ups
        // and the error.
        let mut mask_failure: Option<(Posts, crate::cli::NativeError)> = None;
        if let Some(window_masker) = &args.window_masker {
            // The batch's query source loads the registry of `-window_masker_taxid` before it
            // makes the masks, after the reader's messages of the batch: its posts go with them
            // (`WindowMasker::load_posts`), before every later message of the batch; once for a
            // batch whose masks are made, once for each query when they fail (`batch_error`).
            // The run's first load is that of its first batch: the library's call of a later
            // batch (`api::run_batches`, `first_query_ordinal`) has the loads of those before.
            let first_load = input_start == 0 && args.first_query_ordinal == 0;
            let searched = batch_records.searched();
            let masks = match window_masker.build_error() {
                Some(error) => Err(error.clone()),
                None => window_masker.batch_masks(
                    searched.iter().map(|record| record.seq()),
                    QueryStrand::from_option(&args.strand) == QueryStrand::Minus,
                ),
            };
            match masks {
                Ok(masks) => {
                    batch_title_warnings.append(window_masker.load_posts(first_load));
                    window_masker_masks = Some(masks);
                }
                Err(error) => match window_masker.batch_error(&error, searched.len(), first_load) {
                    Some(failure) => mask_failure = Some(failure),
                    None => return Err(window_masker_rejection(&args, &error)),
                },
            }
        }
        // With `-filtering_db` the batch's masks are made here too (`repeat_batch_masks`), where
        // NCBI sets the batch's queries up; the masks next to each query's `-query_loc` interval
        // that they join go on to the chunks of a split batch (`BatchStage::Search`'s
        // `masks_outside`).
        let mut masks_outside: Vec<MasksOutside> = Vec::new();
        // The batch's repeat hits, kept for the BLAST archive (`BatchResults::repeat_hits`).
        let mut batch_repeat_hits: Vec<Vec<(usize, usize)>> = Vec::new();
        if let Some(filter) = &args.repeat_filter {
            let ranged = (args.lcase_masking && query_input.placements.ranged()).then(|| {
                (start..end)
                    .map(|index| {
                        let input = query_input.ordinal(index) - input_start;
                        (&batch_records.input()[input], query_offsets[index])
                    })
                    .collect()
            });
            let searched = batch_records.searched();
            match repeat_batch_masks(&args, filter, searched, ranged, parallel_pool)? {
                Ok((masks, outside, hits)) => {
                    window_masker_masks = Some(masks);
                    masks_outside = outside;
                    if keeps_repeat_hits {
                        batch_repeat_hits = hits;
                    }
                }
                Err(failure) => mask_failure = Some((Posts::new(), failure)),
            }
        }
        if let Some((masker_warnings, failure)) = mask_failure {
            write_reports_before!();
            batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
            masker_warnings.write(&mut *outputs.diagnostics, archive)?;
            return Err(failure.into());
        }
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:632-652
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1013-1025
        // A query without letters (`Sequence contains no data`, blast_setup.hpp:190-198) has
        // no valid context, and the strands of every other BLASTN query are valid when the
        // queries are set up: a batch of such queries only stops NCBI with one message per
        // query (`GetMessage` puts the severity first), after the reports of the batches
        // before (`CATCH_ALL`: `BLAST engine error: `, exit 3).
        if batch_records
            .searched()
            .iter()
            .all(|record| record.seq().is_empty())
        {
            write_reports_before!();
            batch_title_warnings.write(&mut *outputs.diagnostics, archive)?;
            let mut message = format!("{ENGINE_ERROR} ");
            for _ in start..end {
                message.push_str("Warning: ");
                message.push_str(NO_DATA_MESSAGE);
                message.push(' ');
            }
            message.push('\n');
            return Err(crate::cli::NativeError {
                exit: 3,
                message: message.into(),
            }
            .into());
        }
        warnings[start].append(&batch_title_warnings);
        check_batch_window(&args, batch_records.batch_queries().records())?;
        // NCBI restricts the lower-case runs of each whole record to the chunks of a split
        // batch, and of the letters outside a `-query_loc` interval the one before it reaches
        // them (`query_split::restrict_masks_with_letter_before`); WindowMasker replaces the
        // runs. The input records keep that letter.
        let lower_case_before: Vec<bool> = if args.lcase_masking
            && args.window_masker.is_none()
            && query_input.placements.ranged()
        {
            (start..end)
                .map(|index| {
                    let record = &batch_records.input()[query_input.ordinal(index) - input_start];
                    query_offsets[index]
                        .checked_sub(1)
                        .is_some_and(|letter| record.seq()[letter].is_ascii_lowercase())
                })
                .collect()
        } else {
            Vec::new()
        };
        // NCBI fails a batch where it reads a header file that it cannot read (`read_titles`)
        // or where its prefetch looks up a Seq-align id whose lookup fails on a file
        // (`database::batch_prefetch`), after the prologs, the reports of the batches before
        // (without the epilog) and the batch's reader messages: `error` is NCBI's.
        macro_rules! fail_batch {
            ($error:expr) => {{
                write_reports_before!();
                warnings[start].write(&mut *outputs.diagnostics, archive)?;
                return Err($error);
            }};
        }
        // The batch of LOSAT's refusal (`scope_refusal` above, without a request): NCBI's error of
        // the set-up of its records that are first in the scope (`BatchStage::Search::setup_only`,
        // which ends before the search), or the refusal (line WQA D10 round 5, 2026-11-20).
        // BLASTN's formats have no failure of a matrix name (`QueryScopes::with_xml_matrix_failure`):
        // every refusal stops its batch.
        if let Some(refusal) = scope_refusal.take() {
            if let Err(error) =
                crate::blastinput::query_batch::refused_batch_error(refusal, false, |first| {
                    refused_batch_setup_error(
                        &args,
                        batch_records.searched(),
                        (start..end)
                            .filter(|&query| first[query_input.ordinal(query) - input_start]),
                        start,
                        &query_offsets,
                        window_masker_masks.as_deref(),
                        &lower_case_before,
                        &masks_outside,
                        source,
                        Some((&mut *outputs, &output_formats[..])),
                        parallel_pool,
                        batching,
                        input_start == 0,
                        &mut warnings[start],
                        &mut query_keep,
                    )
                    .into()
                })
            {
                // The first batch's set-up wrote the prologs before its error
                // (`database_search_error`).
                if input_start == 0 && error.is::<crate::cli::NativeError>() {
                    warnings[start].write(&mut *outputs.diagnostics, archive)?;
                    return Err(error);
                }
                fail_batch!(error);
            }
        }
        // Whether the loop reads no batch after this one: its keep then frees the table after the
        // batch's preliminary stage (`BatchStage::Search`'s `last`).
        let last = !batches_follow(
            input_end,
            query_input.input_lengths.len(),
            &reading,
            reading_ended,
        );
        // The batch's search, on the whole pool. EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop keeps
        // its megablast table from batch to batch.
        // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:929
        // The instance of the batch search is chosen here, once for the batch.
        let searched = if args.ungapped {
            search_query_batch::<true>(
                &args,
                batch_records.batch_queries(),
                source,
                Some((&mut *outputs, &output_formats)),
                parallel_pool,
                batching,
                input_start == 0,
                BatchStage::Search {
                    query_offsets: &query_offsets[start..end],
                    query_masks: window_masker_masks.as_deref(),
                    lower_case_before: &lower_case_before,
                    masks_outside: &masks_outside,
                    last,
                    setup_only: false,
                },
                &mut warnings[start],
                x_lut_keep.as_ref(),
                &mut query_keep,
            )
        } else {
            search_query_batch::<false>(
                &args,
                batch_records.batch_queries(),
                source,
                Some((&mut *outputs, &output_formats)),
                parallel_pool,
                batching,
                input_start == 0,
                BatchStage::Search {
                    query_offsets: &query_offsets[start..end],
                    query_masks: window_masker_masks.as_deref(),
                    lower_case_before: &lower_case_before,
                    masks_outside: &masks_outside,
                    last,
                    setup_only: false,
                },
                &mut warnings[start],
                x_lut_keep.as_ref(),
                &mut query_keep,
            )
        };
        let batch = match searched {
            Ok(batch) => batch,
            Err(error) => match error.downcast::<ErrorAfterInvalidBatches>() {
                // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:295-318
                // The reports of the batches before are written, then the error, which
                // skips `PrintEpilog` (the outfmt 0 footer and outfmt 7's last line).
                Ok(ErrorAfterInvalidBatches(error)) => {
                    // EXPERIMENT (LOSAT_X_CSRQSKEEP): as in `write_reports_before`, the keep frees
                    // the `qs` it holds first. No NCBI counterpart (memory only).
                    if let Some(keep) = &x_lut_keep {
                        keep.close_qs();
                    }
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): the reports of the batches before were
                    // written when the loop took them.
                    let before = |outputs: &mut ReportOutputs<'_>| {
                        write_batch_reports(
                            &args,
                            queries.before(&kept, start),
                            ReportSubjects {
                                source,
                                shown: &shown,
                                titles: &titles,
                            },
                            &report_ranges,
                            outputs,
                            &output_formats,
                            BatchResults {
                                hit_lists,
                                warnings: &warnings[..start],
                                reads: &query_reads[..start],
                                lengths: &lengths[..start],
                                unsearched,
                                query_masks: &query_masks,
                                repeat_hits: &repeat_hits,
                                context_karlin: &context_karlin,
                                query_eff_searchsp: &query_eff_searchsp,
                                batches: &batch_spans,
                            },
                            false,
                            None,
                        )
                    };
                    // In the first batch nothing was streamed: the call of the off path writes the
                    // prologs that come with the reports (outfmt 20's line of the field names).
                    if stream.acting() && start > 0 {
                        stream.finish_sam(outputs)?;
                    } else {
                        before(outputs)?;
                    }
                    // The failing batch was read (its title warnings) before its search.
                    warnings[start].write(&mut *outputs.diagnostics, archive)?;
                    return Err(error);
                }
                // The warnings of the batches read so far, in NCBI's order (the outfmt 0
                // prolog of a Karlin-Altschul error of the first batch is written).
                // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were written
                // with their reports, which a refusal of LOSAT's in a later batch follows (the
                // traceback's of ambiguity data that it cannot decode, `x_rpt_stream`); the off
                // path writes no report here.
                Err(error) => {
                    let written = if stream.acting() { start } else { 0 };
                    for query_warnings in &warnings[written..=start] {
                        query_warnings.write(&mut *outputs.diagnostics, archive)?;
                    }
                    return Err(error);
                }
            },
        };
        // The batch's set-up has written its `EffectiveSearchSpace` option, which the next batch's
        // search strategy request holds.
        if let Some(strategy) = args.strategy.as_ref() {
            strategy.note_searched(
                start..end,
                batch.searched.then(|| {
                    context0_searched(
                        &batch.query_eff_searchsp,
                        &batch.context_karlin,
                        lengths[start],
                        0,
                    )
                }),
            )?;
        }
        // The batch's LOSAT_TIMING block, where its search wrote it before.
        if !batch.timing.is_empty() {
            eprint!("{}", batch.timing);
        }
        let mut batch_hit_lists = batch.hit_lists;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:264-275
        // A database: the batch's hit lists are finished now (`finish_query_hit_list`, the
        // end of NCBI's traceback), and the definition lines of the sequences they name are
        // read in NCBI's order where NCBI makes the batch's Seq-aligns
        // (`database::read_titles`), before the HSPs of score 0 are dropped. A header file
        // that NCBI cannot read fails this batch, after the reports of the batches before
        // and the messages of this batch's reading.
        // The taxonomy database's messages of each query's own fetches, where the batch's
        // prefetch is off (below), by query of the batch.
        let mut query_fetches: Vec<Vec<u8>> = Vec::new();
        if let SubjectSource::Database(database) = source {
            let hitlist_size = args.report_counts.hitlist_size;
            for (q_idx, hit_list) in batch_hit_lists.iter_mut().enumerate() {
                if let Some(hit_list) = hit_list {
                    finish_query_hit_list(
                        hit_list,
                        lengths[start + q_idx],
                        hitlist_size,
                        args.max_hsps_per_subject.unwrap_or(0),
                        args.qcov_hsp_perc.unwrap_or(0.0),
                        args.subject_besthit,
                        args.culling_limit_value(),
                        args.best_hit(),
                        args.ungapped,
                        database.sorts_results(),
                    );
                }
            }
            let read = (&mut titles, &mut warnings[start]);
            let one_by_one = crate::algorithm::common::database::prefetch_may_be_off(
                database::prefetch_format(&output_formats),
            );
            if let Err(error) =
                database.batch_titles((batch.fetched_utf8, one_by_one), &mut batch_hit_lists, read)
            {
                let error = match error.downcast::<database::HeaderReadError>() {
                    Ok(failure) => failure.error,
                    // LOSAT's refusal: the warnings of the batches read so far, as for the
                    // other errors of a batch (below).
                    // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were
                    // written with their reports, which the refusal follows (`x_rpt_stream`);
                    // the off path writes no report here.
                    Err(error) => {
                        let written = if stream.acting() { start } else { 0 };
                        for query_warnings in &warnings[written..=start] {
                            query_warnings.write(&mut *outputs.diagnostics, archive)?;
                        }
                        return Err(error);
                    }
                };
                fail_batch!(error);
            }
        }
        // An HSP left of its query (`left_of_query_refusal`): of a subject set NCBI's search of the
        // batch throws, after the prologs, the reports of the batches before and the batch's reader
        // messages; of a database LOSAT refuses after the batch's errors and prefetch (below).
        let database_search = matches!(source, SubjectSource::Database(_));
        let left_of_query = match left_of_query_refusal(
            &args,
            source,
            &batch_hit_lists,
            &lengths[start..],
            &query_offsets[start..],
            batch_records.searched(),
        ) {
            Some(error) if !database_search => fail_batch!(error),
            database_refusal => database_refusal,
        };
        // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:304-309
        // The prefetch before the batch's reports reads `PRE_FETCH_SEQS_LIMIT`, after the batch's
        // search (a database's titles are read in its traceback, above): a value that NCBI cannot
        // convert stops the first batch here, after its reader's messages and before its query
        // warnings (`app::check_prefetch_environment`).
        if let Err(error) = crate::blastinput::app::check_prefetch_environment(
            outputs.only_archive() || args.hits_only,
        ) {
            // EXPERIMENT (LOSAT_X_CSRQSKEEP): as in `write_reports_before`, the keep frees the `qs`
            // it holds first. No NCBI counterpart (memory only).
            if let Some(keep) = &x_lut_keep {
                keep.close_qs();
            }
            write_before_failed_batch(
                &args,
                queries.before(&kept, start),
                ReportSubjects {
                    source,
                    shown: &shown,
                    titles: &titles,
                },
                &report_ranges,
                outputs,
                &output_formats,
                BatchResults {
                    hit_lists,
                    warnings: &warnings[..start],
                    reads: &query_reads[..start],
                    lengths: &lengths[..start],
                    unsearched,
                    query_masks: &query_masks,
                    repeat_hits: &repeat_hits,
                    context_karlin: &context_karlin,
                    query_eff_searchsp: &query_eff_searchsp,
                    batches: &batch_spans,
                },
            )?;
            warnings[start].write(&mut *outputs.diagnostics, archive)?;
            return Err(error);
        }
        // LOSAT's refusal of a reported sequence with the id of a query of its scope, after
        // NCBI's own errors of the batch (the prefetch's above): the warnings of the batches read
        // so far, as for a refusal of a title (above).
        if let (SubjectSource::Database(database), Some(scopes)) = (source, query_scopes.as_ref()) {
            let refused = crate::algorithm::common::database::refuse_reported_ids_of_queries(
                &titles,
                batch_hit_lists
                    .iter()
                    .enumerate()
                    .flat_map(|(query, hit_list)| {
                        hit_list.iter().flat_map(move |hit_list| {
                            let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                            hit_list.hsplist_array[..count]
                                .iter()
                                .map(move |hsp_list| (start + query, hsp_list.oid))
                        })
                    }),
                database.oid_map(),
                scopes,
                database::parsed_shown(&output_formats),
                "BLASTN",
            );
            if let Err(error) = refused {
                // EXPERIMENT (LOSAT_X_RPTSTREAM): those of the batches before were written with
                // their reports.
                let written = if stream.acting() { start } else { 0 };
                for query_warnings in &warnings[written..=start] {
                    query_warnings.write(&mut *outputs.diagnostics, archive)?;
                }
                return Err(error);
            }
        }
        if let SubjectSource::Database(database) = source {
            for hit_list in batch_hit_lists.iter_mut().flatten() {
                drop_hsps_of_score_zero(hit_list);
            }
            // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:304-312
            // Before the reports of a batch NCBI fetches every subject that the batch reports
            // (those with Seq-aligns: the HSPs of score 0 are dropped), and with `taxdb.bti`
            // on the search path the taxonomy database writes its messages there
            // (`BlastnDatabase::prefetch_messages`: the texts of opening it at the run's first
            // lookup, a volume's lookup of a taxid once per run): after the reader's messages
            // of the batch and before the warnings of its first query, which are added below.
            // They go with the warnings, where the first format is outfmt 0, 5, 6, 7 or 17, or
            // with `PRE_FETCH_SEQS_LIMIT=2147483647` 8, 9 or 12, which write the reader's messages
            // only (`BlastnOutputFormat::prefetches_subjects`), so every writer of the batch's
            // reports writes them there (the reports after the loop, LOSAT_X_RPTSTREAM's
            // batch writes, the error paths that write the reports of the batches before).
            // Without `taxdb.bti` nothing is read. A refusal (a taxonomy database or a taxid
            // that LOSAT cannot read as NCBI does) or an error of the taxonomy database is
            // written as the titles' refusals above: with LOSAT_X_RPTSTREAM after the reports of
            // the batches before (`x_rpt_stream`).
            // Where the prefetch fails at a subject with parsed Seq-ids, NCBI writes its warning
            // after those messages (`BlastnDatabase::batch_prefetch`), with or without
            // `taxdb.bti`.
            let taxdb_found = database.taxdb_found();
            // Which rows of the first format read the subjects' Bioseqs where a batch's prefetch
            // is off (`TabularSpec::subject_fetch`; every row of a format that is not a field
            // list): what the reports of such a batch load (`shared::batch_loads`).
            let first_rows = match output_formats.first() {
                Some(BlastnOutputFormat::CustomTabular(spec)) => {
                    spec.subject_fetch(crate::blastinput::app::Program::Blastn)
                }
                _ => crate::report::tabular::SubjectFetch::Every,
            };
            // NCBI's scope at the batch's prefetch (`shared::PrefetchQueries`), made only where it
            // stops (both forms of the prefetch below): the batch's input records (those that
            // `-query_loc` drops included, each with its `Query_N` of the run,
            // `FastaRecord::local_id`), in outfmt 5, whose scope the formatter keeps, those of the
            // batches before too (outfmt 5 keeps the records read before the search, as
            // `query_scopes` above).
            let scope = || {
                crate::blastinput::query_batch::QueryScope::of_shown_records(match &queries {
                    x_query_stream::Queries::Records { input, .. } if xml => &input[..input_end],
                    _ => batch_records.input(),
                })
            };
            let prefetch_queries =
                crate::algorithm::common::database::PrefetchQueries::new(start, &scope)
                    .where_reports_unchanged_by_load_throws(
                        database::reports_unchanged_by_load_throws(&output_formats),
                    );
            let prefetch = if run_loads.ended() {
                None
            } else if output_formats
                .first()
                .is_some_and(|format| format.prefetches_subjects())
                && (taxdb_found || database.parsed_seqids())
            {
                // A custom tabular report with a plural name looks every taxid of its rows up,
                // taxid 0 too, which the prefetch does not: the rows of a batch with hits may open
                // the taxonomy database before a later batch's prefetch
                // (`TaxdbRun::open_for_rows`).
                let plural_names = taxdb_found
                    && matches!(
                        output_formats.first(),
                        Some(BlastnOutputFormat::CustomTabular(spec)) if spec.looks_up_every_taxid()
                    );
                // A query's own fetches where the batch's prefetch is off: the texts of the
                // subjects its report reads, kept by subject where the first format writes them
                // where it reads each (`ReportFormat::places_subject_reads`), else returned for the
                // query's warnings; the pairwise report (outfmt 0 to 4) reads those that it shows
                // as its writer reads them (`database::pairwise_reading`).
                let first_report = outputs
                    .formats
                    .first()
                    .and_then(crate::api::local_blast::FormatOutput::report_format);
                let fetches = database::QueryFetches {
                    database,
                    titles: &titles,
                    run: &run_loads,
                    places: first_report
                        .is_some_and(crate::blastinput::app::ReportFormat::places_subject_reads),
                    plural_names,
                    pairwise: database::pairwise_reading(
                        &output_formats,
                        args.formatting.report_sort(),
                        args.report_counts,
                    ),
                    pairwise_first: first_report
                        .and_then(crate::report::pairwise::PairwiseView::of)
                        .is_some(),
                    gapped: !args.ungapped,
                };
                let seq_aligns_of = |hit_lists: &[Option<BlastnHitList>]| -> usize {
                    database::seq_aligns_of(hit_lists, !args.ungapped)
                };
                let seq_aligns = seq_aligns_of(&batch_hit_lists);
                // Where `PRE_FETCH_SEQS_LIMIT` turns the batch's prefetch off
                // (`shared::prefetch_off`), each query's report fetches the subjects it shows,
                // after the query's own warnings (`query_fetches`, added after them below), where
                // the first format reads their Bioseqs (`TabularSpec::subject_fetch`; BLASTN's
                // formatter does not fetch for the identities, so every row reads them or none
                // does). NCBI looks the taxonomy up only where it fetches a subject: in such a
                // report, or in the prefetch of a batch with Seq-aligns (one without them loads
                // no id, below). The run's refusal of NCBI's taxonomy cache
                // (`shared::TaxdbRun::cache_collision`) is asked only there, so a run that fetches
                // no subject is not refused. Where the report stream acts (`x_rpt_stream`), a
                // later batch's refusal comes after the reports of the batches before, as the
                // stream's other refusals of a batch's reported subjects; without it every
                // report is written after the loop, so the refusal comes before any report.
                let prefetch = if crate::algorithm::common::database::prefetch_off(
                    database::prefetch_format(&output_formats),
                    seq_aligns,
                ) {
                    if !taxdb_found && !database.db.scope_utf8_seen() {
                        // No prefetch, and the reports' fetches look nothing up nor post a warning.
                        Ok(Posts::new())
                    } else {
                        fetches
                            .batch(
                                &batch_hit_lists,
                                start,
                                first_rows != crate::report::tabular::SubjectFetch::Never,
                                &mut query_reads,
                            )
                            .map(|messages| {
                                query_fetches = messages;
                                Posts::new()
                            })
                    }
                } else if seq_aligns == 0 {
                    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:768-772
                    // A batch without Seq-aligns prefetches no sequence.
                    fetches
                        .rows_open(&batch_hit_lists, None)
                        .map(|()| Posts::new())
                } else {
                    // Outfmt 5 first: whether NCBI's XML stream cut a report of the batches
                    // before (`shared::XmlCuts`), asked only where the batch's prefetch stops, by
                    // the record of their loads (`LoadedBefore`: a cut report loads nothing) and
                    // the texts that the cut reports left (`shared::XmlReportsBefore`).
                    let xml_cut = |query: usize| {
                        queries
                            .before(&kept, start)
                            .searched
                            .get(query)
                            .is_some_and(|record| {
                                crate::report::blast_xml::names::iteration_cut_before_hits(
                                    record.local_id.as_bytes(),
                                    record.warning_id(),
                                    record.believed(),
                                    &record.title,
                                    false,
                                )
                            })
                    };
                    // The posts before each query's report: those of the prefetch of the
                    // batch it begins (`warnings`, the batches before this one).
                    let prefetched = |query: usize| warnings.get(query);
                    let xml_cuts = crate::algorithm::common::database::XmlCuts::of(
                        matches!(output_formats.first(), Some(BlastnOutputFormat::Xml)),
                        &xml_cut,
                    );
                    let xml_before = xml_cuts.reports_before(&query_reads, &prefetched);
                    database
                        .batch_prefetch(
                            &titles,
                            database::LoadedBefore {
                                hits: &hit_lists,
                                spans: &batch_spans,
                                format: database::prefetch_format(&output_formats),
                                rows: first_rows,
                                cuts: xml_cuts,
                            },
                            &batch_hit_lists,
                            !args.ungapped,
                            &prefetch_queries.after_xml_reports(xml_before.as_ref()),
                        )
                        .and_then(|messages| {
                            fetches.rows_open(&batch_hit_lists, None)?;
                            Ok(messages)
                        })
                };
                Some(prefetch)
            } else if output_formats
                .first()
                .is_some_and(|format| format.xml2_kind().is_some())
                && database.parsed_seqids()
            {
                // The XML2 and JSON formats write the taxonomy texts at each `Hit`: of the
                // batch's prefetch (where it is on and loads a sequence) only its stop
                // (`BlastnDatabase::xml2_prefetch`, with NCBI's scope at the batch's prefetch as
                // above).
                database.xml2_prefetch(
                    &titles,
                    &hit_lists,
                    &batch_spans,
                    (
                        run_loads.answer(),
                        database::prefetch_format(&output_formats),
                    ),
                    &batch_hit_lists,
                    !args.ungapped,
                    &prefetch_queries,
                    &mut query_reads[start..],
                )
            } else {
                None
            };
            let prefetch_error = match prefetch {
                None => None,
                Some(Ok(posts)) => {
                    warnings[start].append(&posts);
                    None
                }
                Some(Err(error)) if error.is::<database::HeaderReadError>() => {
                    let failure = error
                        .downcast::<database::HeaderReadError>()
                        .expect("the error was checked to be a HeaderReadError");
                    fail_batch!(failure.error);
                }
                Some(Err(error)) => Some(error),
            };
            // The prefetch's error, else LOSAT's refusal of a database's HSP left of its query
            // (above) where NCBI's reports read the batch's Seq-aligns, after its prefetch
            // (blastn_app.cpp:304-312, cited above): after the warnings of the batches read so
            // far with the prefetch's posts.
            if let Some(error) = prefetch_error.or(left_of_query) {
                let written = if stream.acting() { start } else { 0 };
                for query_warnings in &warnings[written..=start] {
                    query_warnings.write(&mut *outputs.diagnostics, archive)?;
                }
                return Err(error);
            }
            // EXPERIMENT (LOSAT_X_RPTSTREAM): the loop writes the batch's reports below and drops
            // its lists, so a later batch's stop does not find them among those it holds (the
            // prefetches' `before`): the run's record of what NCBI's data loader loaded for the
            // batch (`BlastnDatabase::note_loaded`; nothing for a database that cannot stop).
            if stream.acting() {
                database.note_loaded(
                    &titles,
                    &batch_hit_lists,
                    start,
                    !args.ungapped,
                    (database::prefetch_format(&output_formats), first_rows),
                );
            }
        }
        // The batch numbers its queries from 0.
        let first = start as u32;
        hit_lists.extend(batch_hit_lists.into_iter().map(|list| {
            list.map(|mut list| {
                for hsp_list in &mut list.hsplist_array {
                    hsp_list.query_index += first;
                    for hsp in &mut hsp_list.hsps {
                        hsp.q_idx += first;
                    }
                }
                list
            })
        }));
        batch_spans.push(start..end);
        unsearched.extend(std::iter::repeat_n(!batch.searched, end - start));
        query_masks.extend(batch.query_masks);
        repeat_hits.extend(batch_repeat_hits);
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
        // The warning of an invalid query comes with its report (`PrintOneResultSet`),
        // before its preamble (`QueryWarnings`): a query without letters has the message of
        // its set-up (`SetupQueries_OMF`, above), any other invalid query that of its
        // Karlin-Altschul parameters.
        // A structured format (outfmt 8, 12) writes no result warnings
        // (`ReportFormat::writes_query_warnings`); the warnings go with the first format, or, where
        // it is structured, they are kept for the later formats that write them, each to the
        // caller's receiver of its own (`Posts::later`, `ReportOutputs::later_warnings`).
        let result_warnings = outputs.first_writes_query_warnings();
        let keeps_later = outputs.keeps_later_query_warnings();
        let filtered_empty = match source {
            SubjectSource::Database(database) => database.filtered_empty(),
            _ => false,
        };
        // The query's first valid context (`-strand minus`: its minus context).
        let first_context = QueryStrand::from_option(&args.strand).first_context();
        for (index, karlin) in batch
            .context_karlin
            .iter()
            .skip(first_context)
            .step_by(2)
            .enumerate()
        {
            if !result_warnings && !keeps_later {
                continue;
            }
            let mut messages = Vec::new();
            if karlin.is_none() {
                messages.push(if lengths[start + index] == 0 {
                    NO_DATA_MESSAGE.to_string()
                } else {
                    crate::report::query_warnings::INVALID_QUERY_MESSAGE.to_string()
                });
            }
            if batch.searched && filtered_empty {
                messages.push(crate::report::query_warnings::EMPTY_DATABASE_MESSAGE.to_string());
            }
            let text = query_warning(
                "blastn",
                args.first_query_ordinal + query_input.ordinal(start + index),
                &batch_records.searched()[index],
                &messages,
            );
            match result_warnings {
                true => warnings[start + index].lines(&text),
                false => warnings[start + index].later(&text),
            }
        }
        for (index, messages) in query_fetches.into_iter().enumerate() {
            warnings[start + index].lines(&messages);
        }
        context_karlin.extend(batch.context_karlin);
        query_eff_searchsp.extend(batch.query_eff_searchsp);
        // NCBI's `Int4` count of the batch's successful initial extensions.
        if batching.fixed_batch_size == 0 {
            batching.batch_size =
                mixer.batch_size(Some(batch.good_init_extends as u32 as i32)) as u32;
        }
        // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:299-300
        // For a caller that sizes the batches itself (`QueryBatchSizes`).
        if let Some(observer) = outputs.observer.as_deref_mut() {
            observer.batch_extensions(batch.good_init_extends as u32 as i32);
        }
        // EXPERIMENT (LOSAT_X_RPTSTREAM): the batch's reports, written now, after the mixer
        // has its gie, as NCBI's loop prints a batch's result sets (`x_rpt_stream.rs`). The
        // batch's lists and values are taken back from the end of what the loop keeps, numbered
        // from 0 again, and written with the batch's part of the query ranges.
        if stream.acting() {
            let count = end - start;
            let at = hit_lists.len() - count;
            fn take<T>(kept: &mut Vec<T>, at: usize) -> Vec<T> {
                kept.drain(at..).collect()
            }
            let mut batch_lists = take(&mut hit_lists, at);
            for hsp_list in batch_lists
                .iter_mut()
                .flatten()
                .flat_map(|list| list.hsplist_array.iter_mut())
            {
                hsp_list.query_index -= first;
                for hsp in &mut hsp_list.hsps {
                    hsp.q_idx -= first;
                }
            }
            let batch_unsearched = take(&mut unsearched, at);
            let batch_masks = take(&mut query_masks, at);
            // (Empty without the archive: nothing kept.)
            let repeats_at = at.min(repeat_hits.len());
            let batch_repeats = take(&mut repeat_hits, repeats_at);
            let batch_karlin = take(&mut context_karlin, 2 * at);
            let batch_searchsp = take(&mut query_eff_searchsp, 2 * at);
            let batch_input = x_rpt_stream::batch_input(query_input, start..end);
            let batch_ranges = ReportRanges {
                queries: &batch_input,
                subjects: subject_ranges,
            };
            // EXPERIMENT (LOSAT_X_PAIRWISEPAR): the pool formats the batch's outfmt 0 report while
            // it is idle (between the loop's searches).
            let x_pool = Some(parallel_pool);
            stream.write(outputs, count, false, |outputs, written| {
                write_reports(
                    &args,
                    batch_records.report_queries(input_start),
                    ReportSubjects {
                        source,
                        shown: &shown,
                        titles: &titles,
                    },
                    &batch_ranges,
                    outputs,
                    &output_formats,
                    BatchResults {
                        hit_lists: batch_lists,
                        warnings: &warnings[start..end],
                        reads: &query_reads[start..end],
                        lengths: &lengths[start..end],
                        unsearched: batch_unsearched,
                        query_masks: &batch_masks,
                        repeat_hits: &batch_repeats,
                        context_karlin: &batch_karlin,
                        query_eff_searchsp: &batch_searchsp,
                        batches: &[],
                    },
                    false,
                    written,
                    x_pool,
                )
            })?;
            // The written queries' warnings are not read again (the error paths write those
            // of the current batch).
            for query_warnings in &mut warnings[start..end] {
                *query_warnings = Posts::new();
            }
            for reads in &mut query_reads[start..end] {
                *reads = Vec::new();
            }
        }
        // EXPERIMENT (LOSAT_X_QSTREAM, LOSAT_X_BNQPACK): the batch's records are kept for a later
        // writer (moved when the loop took them out of their store), or dropped.
        if keep {
            queries.note_kept(kept.keep(batch_records, query_input.placements.ranged()));
        } else {
            drop(batch_records);
        }
        queries.done(input_end);
        input_start = input_end;
    }
    // EXPERIMENT (LOSAT_X_MBLUTKEEP): the loop's kept table is freed with the loop.
    drop(x_lut_keep);
    // So is its query buffer (the reports below do not read it).
    drop(query_keep);

    let reports = |outputs: &mut ReportOutputs<'_>| {
        write_batch_reports(
            &args,
            queries.all(&kept),
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &report_ranges,
            outputs,
            &output_formats,
            BatchResults {
                hit_lists,
                warnings: &warnings,
                reads: &query_reads,
                lengths: &lengths,
                unsearched,
                query_masks: &query_masks,
                repeat_hits: &repeat_hits,
                context_karlin: &context_karlin,
                query_eff_searchsp: &query_eff_searchsp,
                batches: &batch_spans,
            },
            true,
            Some(parallel_pool),
        )
    };
    if !stream.acting() {
        return reports(outputs);
    }
    // EXPERIMENT (LOSAT_X_RPTSTREAM): every batch's reports are written; the epilogs remain
    // (NCBI's `PrintEpilog` after its loop), written by the report code for no query, with
    // the queries of every batch counted (outfmt 7's last line).
    stream.write(outputs, 0, true, |outputs, written| {
        write_reports(
            &args,
            x_query_stream::ReportQueries::default(),
            ReportSubjects {
                source,
                shown: &shown,
                titles: &titles,
            },
            &report_ranges,
            outputs,
            &output_formats,
            BatchResults {
                hit_lists: Vec::new(),
                warnings: &[],
                reads: &[],
                lengths: &[],
                unsearched: Vec::new(),
                query_masks: &[],
                repeat_hits: &[],
                context_karlin: &[],
                query_eff_searchsp: &[],
                batches: &[],
            },
            true,
            written,
            None,
        )
    })?;
    Ok(())
}

// No NCBI counterpart (when LOSAT reads the records only): it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_QSTREAM): the letters of the largest query batch of a search of `task`, as
/// `run_in_pool` reads its batches (`query_batching_from_environment`): `BATCH_SIZE` when it is
/// set (a `TSeqPos`), else the `CBatchSizeMixer` maximum (`CHUNK_SIZE`, or the task's chunk size,
/// less 1000: 999,000 letters for blastn, 4,999,000 for megablast and dc-megablast; as the loop's
/// `TSeqPos`, so a maximum below 0 is a batch of about 4.29e9 letters); 0 when the variables cannot
/// be read (the search stops at them before its first batch).
fn x_query_batch_letters(task: &str) -> u64 {
    match query_batching_from_environment(super::super::coordination::task_uses_megablast_chunks(
        task,
    )) {
        Ok((batching, _)) if batching.fixed_batch_size != 0 => u64::from(batching.batch_size),
        Ok((batching, _)) => u64::from(batching.split.mixer_max_batch_size() as u32),
        Err(_) => 0,
    }
}

// No NCBI counterpart: LOSAT's verbose, debug and trace output (it changes no value NCBI computes).
/// EXPERIMENT (LOSAT_X_RPTSTREAM): whether a search writes verbose, debug or trace output while its
/// batches are searched (`-verbose`, `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`,
/// `LOSAT_DEBUG_COORDS[_START]`, the blastn trace, the thread diagnostics). The switch then keeps
/// the reports written after the loop, where that output comes before them.
fn x_debug_output(args: &BlastnArgs) -> bool {
    args.verbose
        || crate::utils::xdiag::var_os("BLEMIR_DEBUG").is_some()
        || [
            "LOSAT_DEBUG_BLASTN",
            "LOSAT_DEBUG_COORDS",
            "LOSAT_DEBUG_COORDS_START",
        ]
        .iter()
        .any(|name| std::env::var_os(name).is_some())
        || blastn_trace::enabled()
        || crate::utils::threading::diagnostics_enabled()
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
/// LOSAT's refusals of the reports (a title that NCBI crashes on, `names::Refusal`, DECISIONS
/// 20 as K4h left it) come after the search and before the first byte, so the start is written
/// with the first report, or by the path of an NCBI error that stops the run before one.
///
/// Outfmt 20's line of the field names comes with the first report too (or from
/// `write_pairwise_prologs` for a run that stops before one), after LOSAT's refusals of the
/// custom tabular report; the taxonomy database's warning of its formatter, which a run that
/// can refuse never writes, comes before the search (`tabular::write_tabular_prolog_warning`).
fn prolog_with_report(format: &FormatOutput<'_>, output_format: BlastnOutputFormat) -> bool {
    match output_format {
        BlastnOutputFormat::Pairwise(_) => !prolog_before_search(&format.sink),
        BlastnOutputFormat::Xml2(_) => true,
        BlastnOutputFormat::CustomTabular(spec) => {
            spec.kind == crate::report::tabular::TabularKind::CsvWithHeader
        }
        _ => false,
    }
}

/// What the outfmt 0 report of the search introduces (`BlastnIntro`, with its NCBI reference):
/// megablast's reference, and matrix-only scoring from the option values of reward and penalty
/// (NCBI's `m_DisableKAStats`, set for every task of the BLASTN application).
fn blastn_intro(args: &BlastnArgs) -> BlastnIntro {
    let (reward, penalty, _, _) = super::super::coordination::determine_scoring_params(args);
    BlastnIntro {
        megablast: args.task == "megablast",
        indexed_megablast: args.megablast_index_loaded.is_some(),
        matrix_only: super::super::coordination::matrix_only_scoring(reward, penalty),
    }
}

/// Writes and flushes the outfmt 0 prolog of every outfmt 0 sink whose prolog comes before the
/// search (`before_search`) or, for a run that stops before any report, the prolog of the
/// others (`prolog_with_report`: the file sinks of outfmt 0 and every sink of outfmt 15 and 16).
fn write_pairwise_prologs(
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    intro: BlastnIntro,
    shown: &ShownSubjects,
    before_search: bool,
) -> Result<()> {
    for (format, &output_format) in outputs.formats.iter_mut().zip(output_formats) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:348-357
        // Outfmt 20: the warning of its formatter before the search, the line of the field
        // names with the first report or here before an NCBI error (`prolog_with_report`).
        if let (Some(spec), true) = (output_format.tabular_spec(), before_search) {
            super::super::tabular::write_tabular_prolog_warning(
                spec,
                shown.tabular_taxdb_found,
                &mut *outputs.diagnostics,
            )?;
        }
        let writes_prolog = match output_format {
            BlastnOutputFormat::Pairwise(_) => prolog_before_search(&format.sink) == before_search,
            BlastnOutputFormat::Xml2(_) => !before_search,
            BlastnOutputFormat::CustomTabular(_) => {
                !before_search && prolog_with_report(format, output_format)
            }
            _ => false,
        };
        if writes_prolog {
            let posts = crate::report::query_warnings::html_posts_sink(
                outputs.html_posts.as_deref_mut(),
                &mut *outputs.diagnostics,
            );
            write_format_prolog(format, output_format, intro, shown, posts)?;
        }
    }
    Ok(())
}

/// Writes and flushes the prolog of one sink, for the formats that have one (outfmt 0, 15, 16,
/// 20); `diagnostics` takes the posts of a `-html` prolog's references
/// (`query_warnings::html_posts_sink`).
fn write_format_prolog(
    format: &mut FormatOutput<'_>,
    output_format: BlastnOutputFormat,
    intro: BlastnIntro,
    shown: &ShownSubjects,
    diagnostics: &mut dyn std::io::Write,
) -> Result<()> {
    if output_format.is_pairwise() {
        let flush = crate::report::pairwise::flushes_after_prolog(
            shown.html_links.is_some(),
            !format.sink.is_writer(),
        );
        let mut writer = format.sink.open()?;
        write_blastn_pairwise_prolog(
            &mut writer,
            crate::wording::REPORT_VERSION,
            intro,
            shown.name.as_bytes(),
            // Printed for a subject set only (a database prints `database`).
            i64::from(shown.num_seqs),
            shown.total_letters,
            shown.database.as_ref(),
            shown
                .html_links
                .as_deref()
                .map(|links| (links, diagnostics)),
        )?;
        if flush {
            writer.flush()?;
        }
    } else if let Some(kind) = output_format.xml2_kind() {
        let mut writer = format.sink.open()?;
        writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
        writer.flush()?;
    } else if let Some(spec) = output_format.tabular_spec() {
        // Outfmt 20's line of the field names (nothing for the other tabular kinds).
        if spec.kind == crate::report::tabular::TabularKind::CsvWithHeader {
            let mut writer = format.sink.open()?;
            super::super::tabular::write_tabular_prolog(spec, &mut writer)?;
            writer.flush()?;
        }
    }
    Ok(())
}

/// The results of the query batches searched so far, indexed by query (two contexts per
/// query).
struct BatchResults<'a> {
    hit_lists: Vec<Option<BlastnHitList>>,
    /// The warnings written before each query's report (`QueryWarnings`).
    warnings: &'a [Posts],
    /// The texts of the subjects that each query's report reads one by one
    /// (`BlastnReportInputs::query_reads`).
    reads: &'a [crate::report::query_warnings::QueryReads],
    lengths: &'a [usize],
    unsearched: Vec<bool>,
    query_masks: &'a [Vec<MaskedInterval>],
    /// With `-filtering_db` and the BLAST archive, each query's repeat hits that its batch's
    /// set-up found (`repeat_batch_masks`; 0-based, inclusive, in the searched letters), which
    /// the archive's masks join (`archive::CallBatches::repeat_hits`); empty otherwise.
    repeat_hits: &'a [Vec<(usize, usize)>],
    context_karlin: &'a [Option<ContextKarlin>],
    query_eff_searchsp: &'a [i64],
    /// The searched queries of each query batch of the call, in order, for the reports of
    /// every batch written at once (the BLAST archive writes one document per batch, and each
    /// batch's prefetch counts its own Seq-aligns, `common::database::BatchPrefetch`); empty
    /// when the call is one batch.
    batches: &'a [std::ops::Range<usize>],
}

/// Writes what NCBI has written when a query batch fails, before its error: the reports of
/// the batches before (`before`, the queries searched so far; without the epilog), or, for
/// the first batch, the outfmt 0 prolog of the sinks that write it with the reports (the
/// others wrote it before the search, `write_pairwise_prologs`).
#[allow(clippy::too_many_arguments)]
fn write_before_failed_batch(
    args: &BlastnArgs,
    before: x_query_stream::ReportQueries<'_>,
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
) -> Result<()> {
    if before.searched.is_empty() {
        write_pairwise_prologs(
            outputs,
            output_formats,
            blastn_intro(args),
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
struct ReportSubjects<'a> {
    source: SubjectSource<'a>,
    shown: &'a ShownSubjects,
    titles: &'a database::DatabaseTitles,
}

/// Writes the reports of `query_records` (the queries of the batches searched), with
/// NCBI's epilog (`PrintEpilog`) when `epilog`; `x_pool` as `write_reports`'.
#[allow(clippy::too_many_arguments)]
fn write_batch_reports(
    args: &BlastnArgs,
    query_records: x_query_stream::ReportQueries<'_>,
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
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

/// `write_batch_reports` with what earlier calls wrote (`written`: nothing, except for
/// LOSAT_X_RPTSTREAM's calls, `x_rpt_stream`).
///
/// EXPERIMENT (LOSAT_X_PAIRWISEPAR): `x_pool` is the search's pool, whose threads are idle between
/// the loop's searches; `None` keeps the serial writer.
#[allow(clippy::too_many_arguments)]
fn write_reports(
    args: &BlastnArgs,
    query_records: x_query_stream::ReportQueries<'_>,
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
    queries: x_query_stream::ReportQueries<'_>,
    subjects: ReportSubjects<'_>,
    ranges: &ReportRanges<'_>,
    outputs: &mut ReportOutputs<'_>,
    output_formats: &[BlastnOutputFormat],
    results: BatchResults<'_>,
    epilog: bool,
    mut written: x_rpt_stream::Written<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<()> {
    let query_records = queries.searched;
    // Outfmt 0 only, each such output's refusal (`refuse_shared_masks_in`): the query-anchored
    // sections (outfmt 1 to 4) draw these masks as NCBI does (`query_anchored::fill_loc_list`).
    if args.parse_deflines {
        crate::blastinput::fasta_reader::believed::refuse_shared_masks_in(
            outputs
                .formats
                .iter_mut()
                .zip(output_formats)
                .filter(|(_, output_format)| {
                    **output_format
                        == BlastnOutputFormat::Pairwise(
                            crate::report::pairwise::PairwiseView::Pairwise,
                        )
                })
                .map(|(format, _)| &mut format.refusal),
            query_records,
            |query| {
                results
                    .query_masks
                    .get(query)
                    .is_some_and(|m| !m.is_empty())
            },
            args.lcase_masking,
            "BLASTN",
        )?;
    }
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
    // titles, and their stored data where a report shows letters (outfmt 0, the hit records,
    // BLAST XML, the XML2 and JSON formats and the tabular fields of the alignment rows, such
    // as `sseq`; `database::reported_subjects`).
    let reported = match subjects.source {
        SubjectSource::Prepared(_) => None,
        SubjectSource::Database(database) => {
            let letters = outputs.hits.is_some()
                || output_formats.iter().any(|format| format.is_pairwise())
                || output_formats
                    .iter()
                    .any(|format| format.is_structured_xml())
                || output_formats
                    .iter()
                    .filter_map(|format| format.tabular_spec())
                    .any(|spec| spec.needs().rows);
            // `results.batches`: the batches of this call, or one, whose prefetch decides the
            // subjects that only the prefetch reads.
            Some(database::reported_subjects(
                database,
                subjects.titles,
                &mut hit_lists,
                letters,
                (
                    database::parsed_shown(output_formats),
                    database::prefetch_format(output_formats),
                ),
                (
                    database::after_throw(
                        output_formats,
                        outputs.other_outfmts,
                        args.html_links.is_some(),
                    ),
                    database::pairwise_reading(
                        output_formats,
                        args.formatting.report_sort(),
                        args.report_counts,
                    ),
                ),
                results.batches,
                written.queries,
                !args.ungapped,
                x_pool,
            )?)
        }
    };
    // The library's records name a database's sequences through this event.
    if let (Some(reported), SubjectSource::Database(database), true, Some(observer)) = (
        &reported,
        subjects.source,
        outputs.hits.is_some(),
        outputs.observer.as_deref_mut(),
    ) {
        observer.database_subjects(reported.library_subjects(database)?);
    }
    let database_subject_ids = match (&reported, subjects.source) {
        (Some(reported), SubjectSource::Database(database))
            if output_formats.iter().any(|format| format.is_tabular()) =>
        {
            Some(
                reported.tabular_subject_ids(
                    database.oid_map(),
                    output_formats
                        .iter()
                        .any(|format| format.tabular_reads_subject()),
                )?,
            )
        }
        _ => None,
    };
    let subject_letters = match (&reported, subjects.source) {
        (Some(reported), _) => SubjectLetters::new(&reported.records, reported.stores.as_deref()),
        (None, SubjectSource::Prepared(prepared)) => prepared.letters(),
        (None, SubjectSource::Database(_)) => SubjectLetters::new(&[], None),
    };
    // The archive's `LCaseMask` holds the DUST intervals of every query, valid or not, whether
    // or not its batch was searched; with `-window_masker_db` those of WindowMasker (empty for a
    // query of a batch whose masks failed: its batch writes no archive).
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:456-459
    // (`ExtractUserSpecifiedMasks`, which the archive's request reads, remote_blast.cpp:1358-1359:
    // the query vector's masks after `x_CalculateMasks`, both frames of every interval.)
    let archive_dust = if output_formats.contains(&BlastnOutputFormat::Archive) {
        match &args.window_masker {
            Some(window_masker) => {
                let minus = QueryStrand::from_option(&args.strand) == QueryStrand::Minus;
                query_records
                    .iter()
                    .map(|record| window_masker.masks(record.seq(), minus).unwrap_or_default())
                    .collect()
            }
            None => {
                let serial = crate::utils::threading::SearchPool::x_serial();
                super::super::coordination::apply_dust_masking(
                    args,
                    query_records,
                    x_pool.unwrap_or(&serial),
                )
            }
        }
    } else {
        Vec::new()
    };
    let records_parsed = if args.parse_deflines && reported.is_none() {
        crate::report::db_subject::parsed_records(
            subject_letters
                .records()
                .iter()
                .map(|record| (record.believed(), record.title.as_slice())),
            args.show_gis,
        )
    } else {
        Vec::new()
    };
    // LOSAT's refusals of these records in outfmt 0 to 4, each such format's (`FormatRefusal`).
    if !records_parsed.is_empty() {
        for (format, output_format) in outputs.formats.iter_mut().zip(output_formats) {
            if output_format.is_pairwise() {
                format
                    .refusal
                    .check(crate::report::db_subject::check_pairwise_records(
                        &records_parsed,
                        false,
                        "BLASTN",
                    ))?;
            }
        }
    }
    // The texts of the reads, made before the hit lists named the subjects by rank
    // (`database::reported_subjects`), by rank as the writers name them.
    let query_reads = match &reported {
        Some(reported) => {
            crate::report::query_warnings::reads_by_rank(results.reads, &reported.ordinals)
        }
        None => results.reads.to_vec(),
    };
    // The kept repeat hits are those of the reported queries, one for each (the archive reads
    // them by query index and would take a missing one for a query without hits).
    debug_assert!(
        results.repeat_hits.is_empty() || results.repeat_hits.len() == query_records.len(),
        "repeat hits for {} of {} queries",
        results.repeat_hits.len(),
        query_records.len()
    );
    // A subject set's subjects are all loaded (`database_loads`).
    let every_loaded = crate::report::db_subject::HitLoads::default();
    let report_inputs = BlastnReportInputs {
        query_reads: &query_reads,
        queries: query_records,
        query_records: crate::algorithm::blastn::pairwise::QueryRecords::new(
            query_records,
            ranges.queries,
            queries.input,
            queries.first,
        ),
        subjects: subject_letters,
        ranges,
        query_masks: results.query_masks,
        archive_dust: &archive_dust,
        // With `-filtering_db` the query vector's masks of a query with repeat hits are its masks
        // combined with the hits (`repeat_filter::RepeatFilter::add_masks`): those that the
        // batch's set-up found.
        archive_repeats: results.repeat_hits,
        lcase_masking: args.lcase_masking,
        database_mask: args.db_mask_resolved.as_ref().map(|mask| mask.id),
        query_eff_searchsp: results.query_eff_searchsp,
        first_context: QueryStrand::from_option(&args.strand).first_context(),
        megablast,
        indexed_megablast: args.megablast_index_loaded.is_some(),
        matrix_only: config.matrix_only,
        reward: config.reward,
        penalty: config.penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
        // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:2272
        zero_gap_extension_formula: matches!(args.task.as_str(), "megablast" | "blastn"),
        window_size: config.window_size,
        ungapped: args.ungapped,
        context_karlin: results.context_karlin,
        num_descriptions: args.report_counts.num_descriptions,
        num_alignments: args.report_counts.num_alignments,
        formatting: args.formatting,
        culling_limit: args.culling_limit_value(),
        best_hit: args.best_hit(),
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
        database_parsed: reported
            .as_ref()
            .map_or(records_parsed.as_slice(), |reported| {
                reported.parsed.as_slice()
            }),
        database_sources: reported
            .as_ref()
            .map_or(&[][..], |reported| reported.sources.as_slice()),
        database_loads: reported
            .as_ref()
            .map_or(&every_loaded, |reported| &reported.loads),
        finished_hit_lists: reported.is_some(),
        epilog,
        show_gis: args.show_gis,
        html: args.html_links.clone(),
        is_html: args.html,
        prolog: outputs
            .formats
            .iter()
            .zip(output_formats)
            .any(|(format, &output_format)| {
                output_format.is_pairwise() && !prolog_before_search(&format.sink)
            }),
        queries_before: written.queries,
        subject_cache: written.subjects,
        x_pairwise_pool: x_pool.and_then(crate::report::x_pairwise_par::pool),
        x_tab_pool: x_pool.and_then(crate::algorithm::blastn::x_tab_par::pool),
        out: args.out.as_deref(),
        expect: super::super::coordination::determine_evalue(args),
        dust: args.dust.params().map(|(level, window, linker)| {
            (level as i32, window as u32 as i32, linker as u32 as i32)
        }),
        repeat_filter: args.filtering_db.as_deref(),
        window_masker: args
            .window_masker
            .as_deref()
            .map(super::super::window_masker::WindowMasker::filter_option),
        soft_masking: args.soft_masking(),
        eff_db_length: if args.db_length() > 0 {
            args.db_length()
        } else {
            subjects.shown.total_letters
        },
        archive: args.archive.as_ref(),
        batch_spans: results.batches,
        both_strands: args.strand == "both",
        minus_strand: QueryStrand::from_option(&args.strand) == QueryStrand::Minus,
        tabular_taxdb_found: subjects.shown.tabular_taxdb_found,
        subject_taxdb: args.subject_taxdb.as_ref(),
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:482-496
    // The query ID of the tabular formats is the title's first word, or the local ID
    // (`FastaRecord::shown_id`).
    let query_ids: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.shown_id()))
        .collect();
    let query_titles: Vec<Arc<[u8]>> = query_records
        .iter()
        .map(|record| Arc::from(record.query_line().as_ref()))
        .collect();
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
/// `run_in_pool` gives them to the batches it searches (a split batch to its query chunks).
type BatchOutputs<'o, 'r> = Option<(&'o mut ReportOutputs<'r>, &'o [BlastnOutputFormat])>;

/// Why a search of the first query batch has its outputs (`BatchOutputs`).
const FIRST_BATCH_OUTPUTS: &str =
    "run_in_pool searches the first query batch itself, with the outputs";

/// `outputs` again, for one call.
fn reborrow_outputs<'o, 'r>(outputs: &'o mut BatchOutputs<'_, 'r>) -> BatchOutputs<'o, 'r> {
    outputs
        .as_mut()
        .map(|(outputs, formats)| (&mut **outputs, &**formats))
}

/// The length adjustment and the effective search space of a valid query context of `query_len`
/// letters against a database of `db_len` letters in `db_num_seqs` sequences
/// (`BLAST_CalcEffLengths`), with the context's gapped Karlin-Altschul block. The search space
/// is the context's `eff_searchsp` when no `-searchsp` replaces it. The reports compute the
/// length adjustment again, which the search does not keep.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:821-846
fn context_lengths(
    query_len: i64,
    db_len: i64,
    db_num_seqs: i32,
    gapped: &KarlinParams,
) -> Result<(i64, i64)> {
    let lengths =
        crate::stats::setup::calc_eff_lengths(query_len, db_len, db_num_seqs, gapped, 0, "BLASTN")?;
    Ok((i64::from(lengths.length_adjustment), lengths.eff_searchsp))
}

/// LOSAT's rejection of an error of the masks of `-window_masker_db` or `-window_masker_taxid`
/// that the batch loop does not report as NCBI does (`window_masker::rejection`).
fn window_masker_rejection(
    args: &BlastnArgs,
    error: &crate::utils::winmask::WinMaskError,
) -> anyhow::Error {
    super::super::window_masker::rejection(
        args.window_masker
            .as_deref()
            .map_or("", |masker| masker.name.as_str()),
        error,
    )
}

/// The masks of a query batch with `-filtering_db` (`coordination::query_masks`: DUST's, the
/// lower case and the hits of `filter`'s repeat search of the `searched` letters), and with
/// `-query_loc` and `-lcase_masking` (`ranged`: each query's record and the start of its
/// interval) the lower case next to the intervals that they join (`repeat_filter::MasksOutside`),
/// and the repeat hits of each query (`RepeatFilter::add_masks`), which the batch's BLAST archive
/// reads (`BatchResults::repeat_hits`).
/// NCBI makes them in the set-up of the batch's queries: an error of the repeat search that it
/// reports fails each query of the batch (`Ok(Err(..))`, `repeat_filter::batch_error`).
#[allow(clippy::type_complexity)]
fn repeat_batch_masks(
    args: &BlastnArgs,
    filter: &super::super::repeat_filter::RepeatFilter,
    searched: &[FastaRecord],
    ranged: Option<Vec<(&FastaRecord, usize)>>,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<
    std::result::Result<
        (
            Vec<Vec<MaskedInterval>>,
            Vec<MasksOutside>,
            Vec<Vec<(usize, usize)>>,
        ),
        crate::cli::NativeError,
    >,
> {
    let mut outside: Vec<MasksOutside> = ranged
        .unwrap_or_default()
        .iter()
        .zip(searched)
        .map(|(&(record, offset), part)| {
            MasksOutside::lower_case(record.seq(), offset, part.seq().len())
        })
        .collect();
    match super::super::coordination::query_masks(args, searched, &mut outside, pool) {
        Ok((masks, hits)) => Ok(Ok((masks, outside, hits))),
        Err(MaskError::Repeats(error)) => {
            match super::super::repeat_filter::batch_error(&error, searched.len()) {
                Some(failure) => Ok(Err(failure)),
                None => Err(super::super::repeat_filter::rejection(filter, error)),
            }
        }
        Err(error) => Err(mask_rejection(args, error)),
    }
}

/// LOSAT's rejection of an error of the masks (`coordination::MaskError`) that the batch loop does
/// not report as NCBI does: WindowMasker's (`window_masker_rejection`) and the repeat search's
/// (`repeat_filter::rejection`).
fn mask_rejection(args: &BlastnArgs, error: MaskError) -> anyhow::Error {
    match error {
        MaskError::WindowMasker(error) => window_masker_rejection(args, &error),
        MaskError::Repeats(error) => match args.repeat_filter.as_deref() {
            Some(filter) => super::super::repeat_filter::rejection(filter, error),
            None => anyhow::anyhow!("a repeat search without -filtering_db"),
        },
    }
}

/// `blastinput::app::check_diag_table_window` for a query batch (DI-6, `-window_size`), before
/// `search_query_batch`, whose code a check inside would change (its inlining; the default
/// search's code must stay the same). The concatenated length is that of the batch's contexts
/// (`search_query_batch`): each searched strand of a query with letters takes its letters and
/// a sentinel, and the searched block leaves out the last sentinel.
#[inline(never)]
fn check_batch_window(args: &BlastnArgs, records: &[FastaRecord]) -> Result<()> {
    let window = super::super::coordination::determine_window_size(args);
    let (plus_on, minus_on) = QueryStrand::from_option(&args.strand).searched();
    let strands = usize::from(plus_on) + usize::from(minus_on);
    let rooms: usize = records
        .iter()
        .map(|record| record.sequence.len())
        .filter(|&length| length > 0)
        .map(|length| strands * (length + 1))
        .sum();
    crate::blastinput::app::check_diag_table_window(
        i32::try_from(rooms.saturating_sub(1)).unwrap_or(i32::MAX),
        i32::try_from(window).unwrap_or(i32::MAX),
        "BLASTN",
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321 (cited at
// `SubjectScratch::db_masks`)
/// The ranges of database sequence `s_idx` (length `s_len_full`) of a search with `-db_soft_mask`
/// or `-db_hard_mask` (the algorithm `algo_id`, the kind `kind`): the complement of its stored
/// masks, read with the masks of its block (`BlastnDatabase::seq_ranges_into`) into the worker's
/// buffer, as its soft ranges (`-db_soft_mask`), or as its hard ranges in the worker's
/// `db_hard_split` (`HardSplitState::start`) with no soft ranges (`-db_hard_mask`,
/// blast_engine.c:162-176: NCBI's are the full range, which only its "no chunking" branch reads
/// and LOSAT's hard chunks never do, `next_hard_chunk`); the caller gives the soft ranges' buffer
/// back after the sequence's chunks and ends the hard ones (`HardSplitState::end`). A read error goes to `error_cell` (the first by ordinal id, which the batch
/// returns after the preliminary stage), and `None` ends the sequence's search. Out of line and
/// cold: the per-sequence path of a search without the options only tests them.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn read_db_seq_ranges(
    database: &database::BlastnDatabase,
    algo_id: i32,
    kind: crate::blastinput::db_mask::MaskKind,
    block: Option<&crate::seqdb::BlockRange>,
    s_idx: usize,
    s_len_full: usize,
    subject_scratch: &mut SubjectScratch,
    error_cell: &std::sync::Mutex<Option<(u32, anyhow::Error)>>,
) -> Option<Vec<(i32, i32)>> {
    let hard = kind == crate::blastinput::db_mask::MaskKind::Hard;
    let mut ranges = std::mem::take(if hard {
        &mut subject_scratch.db_hard_split.ranges
    } else {
        &mut subject_scratch.db_soft_ranges
    });
    let read = match block {
        Some(&block) => database.seq_ranges_into(
            block,
            s_idx as u32,
            algo_id,
            s_len_full,
            &mut subject_scratch.db_masks,
            &mut ranges,
        ),
        None => Err(anyhow::anyhow!(
            "LOSAT internal error: database sequence {s_idx} was searched without its block"
        )),
    };
    if let Err(error) = read {
        let mut first = error_cell.lock().expect("the error cell is not poisoned");
        if first
            .as_ref()
            .is_none_or(|(smallest, _)| (s_idx as u32) < *smallest)
        {
            *first = Some((s_idx as u32, error));
        }
        if hard {
            ranges.clear();
            subject_scratch.db_hard_split.ranges = ranges;
        } else {
            subject_scratch.db_soft_ranges = ranges;
        }
        return None;
    }
    if !hard {
        return Some(ranges);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:162-164
    // Every hard chunk carries its range inline, so the full range is not written.
    let mut none = std::mem::take(&mut subject_scratch.db_soft_ranges);
    none.clear();
    subject_scratch
        .db_hard_split
        .start(s_len_full, ranges, true);
    Some(none)
}

/// NCBI's error of the set-up of the records of a refused batch that are first in their scope
/// (`first_in_scope`, of the batch's searched records `batch`, numbered from `start` in the run),
/// which ends before the search (`BatchStage::Search::setup_only`), or `None` where the set-up
/// succeeds or NCBI splits the batch, whose chunks check the subjects (`ContextInputs::split`):
/// the refusal then comes (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20).
fn refused_batch_setup_error(
    args: &BlastnArgs,
    batch: &[FastaRecord],
    first_in_scope: impl Iterator<Item = usize>,
    start: usize,
    query_offsets: &[usize],
    window_masker_masks: Option<&[Vec<MaskedInterval>]>,
    lower_case_before: &[bool],
    masks_outside: &[MasksOutside],
    subjects: SubjectSource<'_>,
    outputs: BatchOutputs<'_, '_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: QueryBatching,
    first_batch: bool,
    posts: &mut Posts,
    query_keep: &mut Vec<u8>,
) -> Option<anyhow::Error> {
    let batch_letters: usize = batch.iter().map(|record| record.seq().len()).sum();
    if !args.ungapped && splits(batching.split, batch_letters) {
        return None;
    }
    let chosen: Vec<usize> = first_in_scope.collect();
    let records: Vec<FastaRecord> = chosen
        .iter()
        .map(|&query| batch[query - start].clone())
        .collect();
    let offsets: Vec<usize> = chosen.iter().map(|&query| query_offsets[query]).collect();
    let masks: Option<Vec<Vec<MaskedInterval>>> = window_masker_masks.map(|masks| {
        chosen
            .iter()
            .map(|&query| masks[query - start].clone())
            .collect()
    });
    let letter_before: Vec<bool> = if lower_case_before.is_empty() {
        Vec::new()
    } else {
        chosen
            .iter()
            .map(|&query| lower_case_before[query - start])
            .collect()
    };
    let outside: Vec<MasksOutside> = if masks_outside.is_empty() {
        Vec::new()
    } else {
        chosen
            .iter()
            .map(|&query| masks_outside[query - start])
            .collect()
    };
    let stage = BatchStage::Search {
        query_offsets: &offsets,
        query_masks: masks.as_deref(),
        lower_case_before: &letter_before,
        masks_outside: &outside,
        last: true,
        setup_only: true,
    };
    let queries = x_query_drop::BatchQueries::Shared(&records);
    let set_up = if args.ungapped {
        search_query_batch::<true>(
            args,
            queries,
            subjects,
            outputs,
            parallel_pool,
            batching,
            first_batch,
            stage,
            posts,
            None,
            query_keep,
        )
    } else {
        search_query_batch::<false>(
            args,
            queries,
            subjects,
            outputs,
            parallel_pool,
            batching,
            first_batch,
            stage,
            posts,
            None,
            query_keep,
        )
    };
    set_up
        .err()
        .map(|error| match error.downcast::<ErrorAfterInvalidBatches>() {
            Ok(ErrorAfterInvalidBatches(error)) => error,
            Err(error) => error,
        })
}

/// Searches one query batch against the subjects, with the query block, lookup table and
/// diagonal table of the batch, as NCBI's `CLocalBlast` does for each batch that
/// `run_in_pool` reads. `x_keep`: the batch loop's keep of LOSAT_X_MBLUTKEEP, which the batch's
/// megablast table is built on and goes back to (`lookup::x_mb_lut_keep`). `posts`: what NCBI has
/// posted for the batch so far (the batch loop's warnings before its first query), to which the
/// search adds what NCBI posts while it searches it (the masks of a split batch's chunks,
/// `search_query_chunks`), also when the search then fails.
fn search_query_batch<const UNGAPPED: bool>(
    args: &BlastnArgs,
    mut query_records: x_query_drop::BatchQueries<'_>,
    subjects: SubjectSource<'_>,
    mut outputs: BatchOutputs<'_, '_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: QueryBatching,
    first_batch: bool,
    stage: BatchStage<'_>,
    posts: &mut Posts,
    x_keep: Option<&Arc<x_mb_lut_keep::Keep>>,
    query_keep: &mut Vec<u8>,
) -> Result<QueryBatch> {
    // NCBI reference (598d8ae6): c++/include/algo/blast/blastinput/blast_args.hpp:1290-1297
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    let num_threads = parallel_pool.threads();

    // The batch's schedule (`batch_schedule::batch_schedule`): the `LOSAT_TIMING` clocks, the thread request,
    // the task's configuration and the speculative traceback's rules.
    let batch_schedule::BatchSchedule {
        timing_enabled,
        t_total,
        timing,
        requested_parallel,
        mut config,
        x_greedy_spec,
        speculative_traceback,
        x_ahead_window,
    } = batch_schedule::batch_schedule(args, num_threads);
    // Read sequences
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:485-651
    // The search reads the query records of its caller (`prepare_sequence_data`).
    let queries = query_records.records();
    // The records' letters for the readers below (`x_query_drop::Letters`).
    let letters = query_records.letters();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:1312-1325
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // Keep one fetched subject record per oid so metadata, lookup filtering,
    // and search reuse the same subject material instead of reparsing FASTA.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1427
    let subject_records: Option<&[FastaRecord]> =
        subjects.prepared().map(|prepared| prepared.records);
    let subject_metadata = subjects.metadata().clone();

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1475
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-500
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
    // The subjects as tasks of the pool; `use_parallel` after LOSAT_X_BNFEWSUBJ's rule
    // (`subject_env::subject_schedule`).
    let subjects_parallel = requested_parallel && subjects.num_subjects() > 1;

    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:836-847
    // The packed letters of each record for the lookup build, which reads them only for
    // `limit_lookup` (not a BLASTN option, always off); `None` with LOSAT_X_BNSTORE.
    let subject_packed_cache_ref = subjects
        .prepared()
        .filter(|prepared| prepared.stores.is_none())
        .map(|prepared| prepared.packed.as_slice());

    // Prepare sequence data (including DUST masking)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1401-1405
    // Diagnostic (LOSAT_TIMING): the batch's DUST for the common line (utils/xtiming.rs); a
    // chunk takes the batch's masks. No NCBI counterpart.
    let t_seg = crate::utils::xtiming::clock();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:445-446
    // The batch loop made the masks of `-window_masker_db` when it set the batch up (the
    // errors stop the batch there); the search takes them.
    let masks = match stage {
        BatchStage::Search {
            query_masks: Some(masks),
            ..
        } => masks.to_vec(),
        BatchStage::Search { .. } => query_masks(args, &letters, &mut [], parallel_pool)
            .map(|(masks, _)| masks)
            .map_err(|error| mask_rejection(args, error))?,
        BatchStage::ChunkPrelim { query_masks, .. } => query_masks.to_vec(),
    };
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Seg, t_seg);
    let seq_data = prepare_sequence_data(queries, masks, subject_metadata);

    // The subjects as the reports show them (the outfmt 0 prolog of an error below).
    let shown = ShownSubjects::new(args, subjects);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2589-2593
    let query_lengths = Arc::new(
        (0..letters.count())
            .map(|index| letters.length(index))
            .collect::<Vec<usize>>(),
    );

    let mut query_block = batch_setup::query_block(args, &letters, &seq_data);
    // The buffer is the one the caller keeps from batch to batch (`x_query_sentinels`).
    let encoded_query_buffer: &mut Vec<u8> = query_keep;
    let nomask_buffer = batch_setup::encode_query_block(
        args,
        &letters,
        parallel_pool,
        &mut query_block,
        encoded_query_buffer,
    );
    let batch_setup::QueryBlock {
        query_contexts,
        query_context_masks,
        query_concat_length,
        diag_container_query_length,
        query_context_offsets,
        query_context_index,
        ..
    } = query_block;
    let encoded_query_concat_blastna_with_sentinels: &[u8] = encoded_query_buffer;
    let encoded_query_concat_blastna: &[u8] = &encoded_query_buffer[1..query_concat_length + 1];
    let encoded_queries_blastna: Vec<&[u8]> = query_contexts
        .iter()
        .map(|ctx| {
            let start = ctx.query_offset.max(0) as usize;
            // The offset of an empty context can be past the block.
            if ctx.len == 0 {
                &[][..]
            } else {
                &encoded_query_concat_blastna[start..start + ctx.len]
            }
        })
        .collect();
    // The letters identities are counted on (NCBI's `sequence_nomask`); without masks written,
    // the searched contexts themselves.
    let nomask_contexts: Option<Vec<&[u8]>> = nomask_buffer.as_ref().map(|buffer| {
        encoded_queries_blastna
            .iter()
            .zip(&query_contexts)
            .map(|(searched, ctx)| {
                let start = ctx.query_offset.max(0) as usize + 1;
                // The offset of an empty context can be past the block.
                if searched.is_empty() {
                    &[][..]
                } else {
                    &buffer[start..start + searched.len()]
                }
            })
            .collect()
    });
    let encoded_queries_nomask: &[&[u8]] = nomask_contexts
        .as_deref()
        .unwrap_or(&encoded_queries_blastna);
    // The same letters concatenated as `encoded_query_concat_blastna` is (NCBI's
    // `sequence_nomask`, blast_filter.c:1382), which an `-ungapped` HSP's identities are counted
    // in: an HSP that NCBI keeps outside its context reads the letters outside it.
    let encoded_query_concat_nomask: &[u8] = nomask_buffer
        .as_deref()
        .map_or(encoded_query_concat_blastna, |buffer| {
            &buffer[1..query_concat_length + 1]
        });

    // Finalize configuration with query-dependent parameters (adaptive lookup table selection)
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:46-48
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:62-67
    let discontig_template = config.mb_template_length > 0;
    // The Karlin-Altschul blocks, alpha/beta and the gap-cost table check score with the engine
    // pair (`TaskConfig::score_reward`: 1/-3 under matrix-only scoring).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:89-96
    let scoring_spec = NuclScoringSpec {
        reward: config.score_reward,
        penalty: config.score_penalty,
        gap_open: config.gap_open,
        gap_extend: config.gap_extend,
    };
    // EXPERIMENT (LOSAT_X_QPREPPAR): with the switch on and a search pool of two or more threads
    // the contexts' blocks are computed by the pool, in context order; the blocks are the same
    // (`x_query_prep_par.rs`).
    let context_ungapped =
        x_query_prep_par::ungapped_blocks(&encoded_queries_blastna, &scoring_spec, parallel_pool);
    let (approx_table_entries, max_q_off) =
        compute_lookup_query_stats(&query_contexts, &query_context_masks);
    let max_query_length = max_q_off.saturating_add(1).max(1);
    finalize_task_config(
        &mut config,
        approx_table_entries,
        max_query_length,
        discontig_template,
    );

    if args.verbose {
        eprintln!(
            "[INFO] Adaptive lookup: approx_entries={}, max_q_off={}, lut_word_length={}, two_stage={}, scan_step={}",
            approx_table_entries,
            max_q_off,
            config.lut_word_length,
            config.use_two_stage,
            config.scan_step
        );
    }
    // The word extension of NCBI's small-query lookup table reads a compressed query in
    // which ambiguity codes are bases (`build_compressed_query`); the other tables compare
    // the query letters themselves.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_nalookup.c:237-240
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/na_ungapped.c:1806-1819
    let small_na_compressed_query =
        if config.small_na_lookup && config.effective_word_size > config.lut_word_length {
            build_compressed_query(encoded_query_concat_blastna)
        } else {
            Vec::new()
        };
    let small_na_aligned_one_byte = config.lut_word_length.is_multiple_of(COMPRESSION_RATIO)
        && config.scan_step.is_multiple_of(COMPRESSION_RATIO)
        && config.effective_word_size - config.lut_word_length <= 4;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1002-1003
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
    // `query->length` is the whole query block (both strands and the sentinel, the
    // offsets of the word finder), not one strand: a table sized from one strand made the
    // diagonals of the two strands share entries and lost HSPs (S07+ independent audit).
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:133-142
    // The two-hit window of the task or of `-window_size` (`configure_task`; a window whose sum
    // with the query would pass 2^30 is refused before this search, `check_batch_window`).
    let window_size = config.window_size;
    let diag_array_length = {
        let mut diag_array_length = 1usize;
        while diag_array_length < query_concat_length + window_size {
            diag_array_length <<= 1;
        }
        diag_array_length
    };
    let diag_mask = diag_array_length - 1;

    let intro = blastn_intro(args);
    // Whether NCBI searches the batch chunk by chunk when it has a valid context (a gapped
    // search, `search_query_chunks`): each chunk's set-up then checks the subjects
    // (`ContextInputs::split`).
    let split = matches!(stage, BatchStage::Search { .. })
        && !UNGAPPED
        && splits(
            batching.split,
            (0..letters.count())
                .map(|index| letters.length(index))
                .sum(),
        );
    // The Karlin-Altschul blocks of the batch's contexts, the checks of the set-up and the
    // effective search spaces (`context_statistics::context_statistics`).
    let context_statistics::ContextStatistics {
        context_karlin,
        round_down_evalue_score,
        context_valid,
        any_valid_context,
        search_karlin,
        query_eff_searchsp,
    } = context_statistics::context_statistics::<UNGAPPED>(
        &context_statistics::ContextInputs {
            args,
            context_ungapped: &context_ungapped,
            scoring_spec,
            query_contexts: &query_contexts,
            query_lengths: &query_lengths,
            batching,
            first_batch,
            stage,
            split,
            subjects,
            seq_data: &seq_data,
            intro,
            shown: &shown,
        },
        reborrow_outputs(&mut outputs),
    )?;
    if let BatchStage::Search {
        setup_only: true, ..
    } = stage
    {
        // The search's first read of the database after the set-up (`prelim_pass`): a sequence
        // file that cannot be opened fails the batch there (line WQA D10 round 6, review 4 M2).
        if let (true, SubjectSource::Database(database)) = (any_valid_context && !split, subjects) {
            if let Some(error) = database.sequence_file_read_error() {
                return Err(database_search_error(
                    error,
                    first_batch,
                    reborrow_outputs(&mut outputs),
                    intro,
                    &shown,
                ));
            }
        }
        return Err(SetUpOnly.into());
    }
    // The preliminary stage's option spaces and its own spaces where they differ from the setup's
    // (`batch_setup::batch_spaces`).
    let batch_setup::BatchSpaces {
        per_subject,
        option_spaces,
        prelim_eff_searchsp,
    } = batch_setup::batch_spaces(&batch_setup::SpaceInputs {
        args,
        stage,
        subjects,
        query_contexts: &query_contexts,
        context_valid: &context_valid,
        search_karlin: &search_karlin,
        query_eff_searchsp: &query_eff_searchsp,
    })?;

    // NCBI searches a batch of at least two query chunks chunk by chunk in its preliminary
    // stage, after it has set up the whole batch (the errors and warnings above are the
    // batch's), and runs the traceback on the batch. A batch without a valid context is not
    // searched.
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:232-236
    let split_prelim_lists = match stage {
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:60-61
        BatchStage::Search {
            query_offsets,
            lower_case_before,
            masks_outside,
            last,
            ..
        } if any_valid_context && !UNGAPPED => {
            let chunks = search_query_chunks(
                args,
                &letters,
                &seq_data.query_masks,
                query_offsets,
                lower_case_before,
                masks_outside,
                option_spaces.as_ref().unwrap_or(&query_eff_searchsp),
                &query_context_offsets,
                subjects,
                reborrow_outputs(&mut outputs),
                parallel_pool,
                batching,
                first_batch,
                posts,
                (x_keep, last),
            )?;
            // A batch with a valid context whose chunks have none (the chunks' own masks would
            // cover the letters that the batch's leave open) meets `BLAST_GapAlignSetUp` first
            // in the set-up of its traceback, after the chunks:
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:231-253
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1831-1837
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback_mt_priv.c:219-226
            // The greedy set-up refuses there as in a batch that is not split. The check of the
            // subjects before it (NCBI: `Traceback failed`, exit status 3) is not reproduced: no
            // input that reaches this was found (line RG D1a-6, review A F1).
            if let Some((_, false)) = &chunks {
                if let Some(error) = context_statistics::greedy_setup_refusal(args, false) {
                    return Err(context_statistics::setup_failure(
                        error,
                        first_batch,
                        reborrow_outputs(&mut outputs),
                        intro,
                        &shown,
                    ));
                }
            }
            chunks.map(|(lists, _)| lists)
        }
        _ => None,
    };
    // EXPERIMENT (LOSAT_X_BNQPACK): the encoding and the query chunks were the last readers of the
    // batch's letters; records held by the caller alone drop them here, before the lookup tables
    // are built (`x_query_drop.rs`). The search reads only the records' count from here on.
    let SequenceData {
        queries: _,
        query_masks,
        db_len_total,
        db_num_seqs,
    } = seq_data;
    query_records.after_encoding();
    let seq_data = SequenceData {
        queries: query_records.records(),
        query_masks,
        db_len_total,
        db_num_seqs,
    };

    // Whether the batch is searched with a MegaBLAST database index (`index_seeds` below; only a
    // database search keeps one, `mb_index::MegablastIndex::check`), and whether the index leaves
    // database volumes to the scan (a mixed search).
    let index_search = args.megablast_index_loaded.is_some() && split_prelim_lists.is_none();
    let index_mixed = args
        .megablast_index_loaded
        .as_deref()
        .is_some_and(crate::dbindex::search::LoadedIndex::mixed);
    // Build lookup tables
    let (lookup_tables, scan_step, offset_array_size) =
        batch_setup::batch_lookup_tables(&batch_setup::LookupInputs {
            config: &config,
            args,
            encoded_queries: &encoded_queries_blastna,
            query_context_masks: &query_context_masks,
            query_context_offsets: &query_context_offsets,
            subject_records,
            subject_packed_cache: subject_packed_cache_ref,
            approx_table_entries,
            x_keep,
            subjects,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:110-163
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:460-462
            // A search with a MegaBLAST database index has no lookup table, except a mixed
            // one (`eMixedMBLookupTable`: the table without an index). (NCBI also builds the
            // table for an old-style index found after the new style left a database volume
            // without an index, `partial`, and never reads it: every OID is the index's.)
            split: split_prelim_lists.is_some() || (index_search && !index_mixed),
        });
    // EXPERIMENT (LOSAT_X_CSRQSKEEP): the loop builds no table after these
    // (`BatchStage::x_last_tables`), so its keep frees the `qs` it still holds (these tables are
    // chains, which take none) before their scan, and takes none back from here on: a last CSR took
    // the kept `qs` and frees its `qs` when its tables drop, as before the switch. No NCBI
    // counterpart (memory only).
    if let (true, Some(keep)) = (stage.x_last_tables(), x_keep) {
        keep.close_qs();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:654 (cutoff_score_min computed per subject)
    // Cutoff scores are computed per subject in the main loop; no global map needed.

    if args.verbose {
        eprintln!("Searching...");
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-501
    let t_search_start = if timing_enabled {
        Some(std::time::Instant::now())
    } else {
        None
    };
    // The values of the batch's subject search (`subject_env::subject_search_env`).
    let subject_env::SubjectSearchEnv {
        progress_bar,
        hitlist_size,
        prelim_hitlist_size,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        verbose,
        debug_window,
        lcase_masking,
        effective_word_size,
        min_ungapped_score,
        use_dp,
        use_direct_lookup,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_gapped,
        x_drop_final,
        scan_range,
        min_diag_separation,
        // Read only by the chunk strips of `prelim_subject::collect_subject_prelim`, which a build
        // without the parallel feature (or serial Wasm) does not have.
        #[cfg_attr(
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
            expect(
                unused_variables,
                reason = "read only by the chunk strips, built with the parallel feature outside serial Wasm"
            )
        )]
        x_strips,
        #[cfg_attr(
            not(all(
                feature = "parallel",
                any(not(target_arch = "wasm32"), feature = "wasm-threads")
            )),
            expect(
                unused_variables,
                reason = "read only by the chunk strips, built with the parallel feature outside serial Wasm"
            )
        )]
        x_kind,
        use_parallel,
        evalue_threshold,
        subject_besthit,
        x_qgroup_batch,
        percent_identity,
        score_matrix,
        nucl_score_table,
        x_dropoff_init,
    } = subject_env::subject_search_env::<UNGAPPED>(&subject_env::SubjectEnvInputs {
        args,
        config: &config,
        subjects,
        subjects_parallel,
        parallel_pool,
        num_threads,
        two_stage_lookup_ref: lookup_tables.two_stage_lookup.as_ref(),
        diag_container_query_length,
        window_size,
        small_na_compressed_query: &small_na_compressed_query,
        context_karlin: &context_karlin,
        search_karlin: &search_karlin,
    })?;

    // The queries and their contexts for `PrelimBatch` and `TraceBatch` (the lookup tables are
    // the preliminary pass's, `SubjectPasses::prelim_pass`)
    let queries_ref = &seq_data.queries;
    let query_contexts_ref = &query_contexts;
    let query_eff_searchsp_ref = &query_eff_searchsp;

    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    let search_karlin_ref = &search_karlin;
    // The cutoffs of every context of the batch, those of the preliminary stage's own spaces where
    // they differ (`prelim_eff_searchsp`), and the statistics of one subject
    // (`batch_setup::batch_cutoff_set`, `batch_setup::subject_statistics`).
    let cutoff_inputs = batch_setup::CutoffInputs {
        subjects,
        db_len_total: seq_data.db_len_total,
        db_num_seqs: seq_data.db_num_seqs,
        query_contexts: &query_contexts,
        context_valid: &context_valid,
        query_eff_searchsp: &query_eff_searchsp,
        search_karlin: search_karlin_ref,
        x_dropoff_init: &x_dropoff_init,
        min_raw_gapped_score: config.min_raw_gapped_score,
        matrix_only: config.matrix_only,
        evalue_threshold,
        subject_length: None,
    };
    let batch_setup::BatchCutoffSet {
        batch: x_batch_cutoffs,
        prelim: x_prelim_cutoffs_owned,
    } = batch_setup::batch_cutoff_set::<UNGAPPED>(&cutoff_inputs, prelim_eff_searchsp.as_deref())?;
    let prelim_eff_searchsp_ref: &Vec<i64> =
        prelim_eff_searchsp.as_ref().unwrap_or(&query_eff_searchsp);
    let x_prelim_cutoffs = x_prelim_cutoffs_owned.as_ref().unwrap_or(&x_batch_cutoffs);
    let traceback_option_spaces = vec![args.options_searchsp(); query_contexts.len()];
    let subject_space_overflow = crate::stats::setup::SubjectSpaceOverflow::default();
    let subject_statistics =
        batch_setup::subject_statistics::<UNGAPPED>(cutoff_inputs, args.db_length());
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:478-536
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:301-310
    // The initial hits that the word finder saved in every subject chunk of the batch; they
    // size the next batch.
    let good_init_extends = std::sync::atomic::AtomicU64::new(0);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:318-321 (cited at
    // `SubjectScratch::db_masks`)
    // A database search with `-db_soft_mask` or `-db_hard_mask` (the formatter's `ResolvedMask`):
    // the database, the algorithm whose masks every database sequence's preliminary search reads
    // (`BlastnDatabase::seq_ranges_into`) and the kind, and the first error of those reads by
    // ordinal id, which the batch returns after the preliminary stage (`search_subject` returns
    // nothing).
    let db_mask: Option<(
        &database::BlastnDatabase,
        i32,
        crate::blastinput::db_mask::MaskKind,
    )> = match (subjects, args.db_mask_resolved.as_ref()) {
        (SubjectSource::Database(database), Some(mask)) => Some((database, mask.id, mask.kind)),
        _ => None,
    };
    // The algorithm of `-db_hard_mask`, whose masked letters the traceback and the ungapped
    // re-evaluation read as 'B' (`database::mask_hard_letters`; decision H2).
    let db_hard_mask: Option<i32> = db_mask
        .filter(|&(_, _, kind)| kind == crate::blastinput::db_mask::MaskKind::Hard)
        .map(|(_, algo_id, _)| algo_id);
    let db_mask_error: std::sync::Mutex<Option<(u32, anyhow::Error)>> = std::sync::Mutex::new(None);
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:302-306
    // The seeds of the MegaBLAST database index for the batch's query (`mb_index_seeds`), the
    // index volumes that the preliminary stage reaches searched here, one after another on the
    // whole pool (`IndexSeeds::search_volumes`); a batch split into query chunks has none (each
    // chunk searched its own above).
    let index_seeds = match (subjects, args.megablast_index_loaded.as_deref()) {
        (SubjectSource::Database(database), Some(index)) if index_search => {
            Some(mb_index_seeds::IndexSeeds::new(
                index,
                database.oid_map(),
                encoded_query_concat_blastna,
                &query_contexts,
                &context_valid,
                crate::dbindex::search::SearchOptions {
                    word_size: effective_word_size as u64,
                    two_hits: crate::algorithm::blastn::coordination::determine_window_size(args)
                        as u64,
                },
                db_hard_mask.map(|algo_id| (database, algo_id)),
            ))
        }
        _ => None,
    };
    if let Some(index) = &index_seeds {
        index.search_volumes(parallel_pool);
    }
    // What the preliminary search of a subject reads of the batch
    // (`prelim_subject::PrelimBatch`), built once: the preliminary stage's cutoffs and
    // spaces (the setup's, or its own where they differ, `prelim_eff_searchsp`). With the
    // statistics of each subject, `SubjectPasses::special_subject` copies it with the subject's.
    let prelim_batch = prelim_subject::PrelimBatch {
        query_context_index: &query_context_index,
        query_concat_length,
        diag_container_query_length,
        encoded_query_concat_blastna,
        encoded_query_concat_blastna_with_sentinels,
        encoded_queries_blastna: &encoded_queries_blastna,
        query_lengths: &query_lengths,
        queries: queries_ref,
        query_contexts: query_contexts_ref,
        query_eff_searchsp: prelim_eff_searchsp_ref,
        cutoff_scores: &x_prelim_cutoffs.0,
        hit_saving_cutoff_scores: &x_prelim_cutoffs.1,
        x_dropoff_scores: &x_prelim_cutoffs.2,
        reduced_cutoff_scores: &x_prelim_cutoffs.3,
        timing_ref: timing.as_deref(),
        effective_word_size,
        scan_step,
        scan_range,
        window_size,
        offset_array_size,
        use_direct_lookup,
        small_na_compressed_query: &small_na_compressed_query,
        small_na_aligned_one_byte,
        min_ungapped_score,
        diag_array_length,
        diag_mask,
        score_matrix: &score_matrix,
        nucl_score_table: &nucl_score_table,
        reward,
        penalty,
        gap_open,
        gap_extend,
        use_dp,
        x_drop_gapped,
        min_diag_separation,
        good_init_extends: &good_init_extends,
        requested_parallel,
        use_parallel,
        num_threads,
        parallel_pool,
        x_qgroup_batch,
        x_ahead_window,
        x_greedy_spec,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        x_strips,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        x_kind,
        subject_besthit,
        lcase_masking,
        db_mask,
        db_mask_error: &db_mask_error,
        mb_index: index_seeds.as_ref(),
        encoded_query_concat_nomask,
        search_karlin_ref,
        round_down_evalue_score,
        prelim_reap: PrelimReap::new(
            config.matrix_only,
            config.min_raw_gapped_score,
            evalue_threshold,
        ),
        percent_identity,
        matrix_only: config.matrix_only,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
        verbose,
        debug_window,
    };
    // What the traceback of a subject reads of the batch (`traceback_subject::TraceBatch`), built
    // once: the setup's cutoffs and spaces, or, in an ungapped search, which has no traceback, its
    // preliminary stage's (blast_traceback.c:1526, `perform_traceback`). With the statistics of
    // each subject, `SubjectPasses::special_subject` copies it with the subject's.
    let (trace_cutoffs, trace_eff_searchsp) = if UNGAPPED {
        (x_prelim_cutoffs, prelim_eff_searchsp_ref)
    } else {
        (&x_batch_cutoffs, query_eff_searchsp_ref)
    };
    let trace_batch = traceback_subject::TraceBatch {
        // A batch split into query chunks traces what its chunks' indexed searches kept. (With
        // an index of only some volumes this is 1 for every subject: the traceback sees only
        // subjects with kept HSPs, and a scanned subject shorter than a word has none.)
        shortest_subject: prelim_subject::shortest_subject(
            effective_word_size,
            args.megablast_index_loaded.is_some(),
        ),
        encoded_queries_blastna: &encoded_queries_blastna,
        encoded_queries_nomask,
        encoded_query_concat_blastna,
        encoded_query_concat_nomask,
        query_lengths: &query_lengths,
        query_concat_length,
        search_karlin_ref,
        round_down_evalue_score,
        queries: queries_ref,
        query_contexts: query_contexts_ref,
        query_eff_searchsp: trace_eff_searchsp,
        cutoff_scores: &trace_cutoffs.0,
        hit_saving_cutoff_scores: &trace_cutoffs.1,
        timing_ref: timing.as_deref(),
        score_matrix: &score_matrix,
        reward,
        penalty,
        gap_open,
        gap_extend,
        x_drop_final,
        use_dp,
        percent_identity,
        min_diag_separation,
        evalue_threshold,
        speculative_traceback,
        x_greedy_spec,
        x_ahead_window,
        num_threads,
        parallel_pool,
        x_qgroup_batch,
        progress_bar: &progress_bar,
        debug_mode,
        blastn_debug,
        blastn_trace_enabled,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1434-1443 and
    // c++/src/algo/blast/core/blast_traceback.c:1598-1628 (the statistics of each subject,
    // `subject_statistics`): with the statistics of each subject (`per_subject`, whose option
    // spaces both stages have), a database's subjects read their own; `SubjectPasses` then takes
    // the subject's call apart.
    let special_statistics =
        option_spaces
            .as_deref()
            .filter(|_| per_subject)
            .map(|option_spaces| subject_passes::SpecialStatistics {
                option_spaces,
                traceback_option_spaces: &traceback_option_spaces,
                subject_statistics: &subject_statistics,
                subject_space_overflow: &subject_space_overflow,
            });
    // The batch's first fetch that read a UTF8String header (`SubjectPasses::fetched_utf8`).
    let fetched_utf8 = crate::seqdb::Utf8Reads::default();
    // The passes over the batch's subjects (`subject_passes::SubjectPasses`), built once.
    let passes = subject_passes::SubjectPasses {
        args,
        subjects,
        prelim_batch: &prelim_batch,
        trace_batch: &trace_batch,
        any_valid_context,
        special: special_statistics.as_ref(),
        prelim_hitlist_size,
        hitlist_size,
        db_hard_mask,
        fetched_utf8: &fetched_utf8,
    };

    // The batch's results, read after its subjects are searched.
    let batch_result = |hit_lists: Vec<Option<BlastnHitList>>,
                        chunk_prelim_lists: PrelimHitLists| QueryBatch {
        hit_lists,
        chunk_prelim_lists,
        good_init_extends: good_init_extends.load(std::sync::atomic::Ordering::Relaxed),
        searched: any_valid_context,
        query_masks: seq_data.query_masks.clone(),
        context_karlin: context_karlin.clone(),
        query_eff_searchsp: query_eff_searchsp.clone(),
        timing: String::new(),
        fetched_utf8: fetched_utf8.first().map(|note| note.at),
    };
    // The preliminary pass over the subjects (`SubjectPasses::prelim_pass`); its error is the database's
    // read error as it is, reported before the masks' read errors below.
    let prelim_lists = passes
        .prelim_pass::<UNGAPPED>(split_prelim_lists, &lookup_tables)
        .map_err(|error| {
            database_search_error(
                error,
                first_batch,
                reborrow_outputs(&mut outputs),
                intro,
                &shown,
            )
        })?;
    // The first error of the database masks' reads (`db_mask_error`), as the preliminary stage's
    // read errors are reported (`database_search_error`).
    let db_mask_error = db_mask_error
        .lock()
        .expect("the error cell is not poisoned")
        .take();
    if let Some((_, error)) = db_mask_error {
        return Err(database_search_error(
            error,
            first_batch,
            reborrow_outputs(&mut outputs),
            intro,
            &shown,
        ));
    }
    // An index volume that did not read or search is LOSAT's refusal, after the preliminary
    // stage's read errors.
    if let Some(index_seeds) = &index_seeds {
        index_seeds.check()?;
    }
    // An overflow of a subject's search space (`subject_statistics`) is LOSAT's refusal, after the
    // preliminary stage's read errors (the searched blocks' above, the masks' here).
    subject_space_overflow.check()?;
    if let BatchStage::ChunkPrelim { .. } = stage {
        return Ok(batch_result(Vec::new(), prelim_lists));
    }
    // LOSAT_X_BNLUTDROP (graduated 2026-10-27, always on): the traceback reads the kept HSPs only,
    // so the tables are dropped here (NCBI's traceback stage, traceback_stage.cpp:208-211)
    // (`x_lut_drop.rs`).
    // The keep of LOSAT_X_MBLUTKEEP frees its table too when no later batch builds one on it.
    x_lut_drop::after_prelim(
        lookup_tables,
        x_keep.filter(|_| stage.x_last_tables()).map(|keep| &**keep),
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:111-114,1605
    // NCBI frees the initial (ungapped) hit list at the end of its preliminary search. No NCBI
    // counterpart (memory placement only): the idle scratch of the ungapped hits' sort (`X_SORT`)
    // serves the preliminary stage only, so it is freed here and the traceback runs without its
    // touched pages (Sakai x EDL933 `-subject` blastn: 11.4 MiB); the next batch's first sort maps
    // one again.
    X_SORT.x_clear_idle();
    // The traceback pass over the subjects (`SubjectPasses::traceback_pass`).
    let traced = passes.traceback_pass::<UNGAPPED>(prelim_lists);
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-279
    // No NCBI counterpart: NCBI writes nothing to its standard error while it searches a batch.
    // The bar moves only in the traceback pass; its end clears it
    // (`subject_env::batch_progress_bar`) whether the pass ends or fails, before anything else is
    // written (the batch's report, or the prologs of a database error), so nothing of it stays on
    // the terminal.
    if let Some(bar) = progress_bar.as_ref() {
        bar.finish_using_style();
    }
    let hit_lists = traced.map_err(|error| {
        database_search_error(
            error,
            first_batch,
            reborrow_outputs(&mut outputs),
            intro,
            &shown,
        )
    })?;
    subject_space_overflow.check()?;
    let mut batch = batch_result(hit_lists, Vec::new());
    if let Some(timing) = timing.as_ref() {
        batch.timing = blastn_timing_text(timing.as_ref(), t_search_start, t_total);
    }
    Ok(batch)
}

/// NCBI's preliminary search of a split query batch: each query chunk is searched as a batch
/// of its query parts (`BatchStage::ChunkPrelim`), and its preliminary hit lists are mapped
/// onto the batch's contexts and merged with those of the chunks before it
/// (`merge_query_chunk`). The merged hit lists of the batch's queries; `None` when NCBI does
/// not split the batch.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:237-289
/// A chunk has its own query block, masks, Karlin-Altschul blocks, lookup table and
/// diagnostics (a split batch gives the next batch size 0 initial hits), and NCBI drops its
/// warnings. A chunk without a valid context finds nothing, as NCBI ignores its error. With the
/// merged lists it returns whether a chunk had a valid context (its set-up checked the subjects
/// and the greedy set-up). `x_keep`: the batch loop's keep (LOSAT_X_MBLUTKEEP) and whether the
/// batch is the loop's last (`BatchStage::Search`'s `last`; LOSAT_X_CSRQSKEEP, the chunks'
/// `x_last`).
#[allow(clippy::too_many_arguments)]
fn search_query_chunks(
    args: &BlastnArgs,
    queries: &x_query_drop::Letters<'_>,
    query_masks: &[Vec<MaskedInterval>],
    query_offsets: &[usize],
    lower_case_before: &[bool],
    masks_outside: &[MasksOutside],
    batch_eff_searchsp: &[i64],
    context_offsets: &[i32],
    subjects: SubjectSource<'_>,
    mut outputs: BatchOutputs<'_, '_>,
    parallel_pool: &crate::utils::threading::SearchPool<'_>,
    batching: QueryBatching,
    first_batch: bool,
    posts: &mut Posts,
    (x_keep, last_batch): (Option<&Arc<x_mb_lut_keep::Keep>>, bool),
) -> Result<Option<(PrelimHitLists, bool)>> {
    let lengths: Vec<usize> = (0..queries.count())
        .map(|index| queries.length(index))
        .collect();
    // A negative CHUNK_SIZE splits a batch only with a negative OVERLAP_CHUNK_SIZE below it
    // (`calculate_num_chunks`); NCBI's `size_t` chunk ranges then wrap (split_query_cxx.cpp:
    // 145-171): NCBI stops with a CCoreException that names its build's files where a chunk
    // would be split again, and otherwise searches chunk ranges with gaps between them.
    if batching.split.negative_chunk_size() && splits(batching.split, lengths.iter().sum()) {
        anyhow::bail!(crate::algorithm::blastn::query_split::negative_split_refusal("BLASTN"));
    }
    let Some(chunks) = split_query_batch(
        &lengths,
        batching.split,
        QueryStrand::from_option(&args.strand),
    ) else {
        return Ok(None);
    };
    // With an overlap close to the chunk size NCBI's setup of a chunk would split it again
    // and stop with a CCoreException that names its build's files. LOSAT searches each chunk
    // once, as NCBI does for every overlap that does not split a chunk again (approved
    // exception 1 of PD-LOSAT-NCBI-DEFECTS).
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:190-201
    // (each chunk's setup builds its own query splitter; only a debug build asserts.)
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-208
    let mut merged: PrelimHitLists = Vec::with_capacity(queries.count());
    merged.resize_with(queries.count(), || None);
    // No NCBI counterpart (memory placement only, `x_query_sentinels`): the chunks' query buffer,
    // kept from chunk to chunk (the batch's own buffer is in use until its traceback).
    let mut query_keep: Vec<u8> = Vec::new();
    // Whether a chunk had a valid context, so that its set-up (`context_statistics`) checked the
    // subjects and the greedy set-up.
    let mut chunk_set_up = false;
    for (index, chunk) in chunks.iter().enumerate() {
        // NCBI fails to make the query factory of a chunk without a query ("Empty
        // CBlastQueryVector"); `split_query_batch` makes none.
        anyhow::ensure!(
            !chunk.queries.is_empty(),
            "a query chunk without a query is not supported by LOSAT's BLASTN"
        );
        let parts: Vec<FastaRecord> = chunk
            .queries
            .iter()
            .map(|part| queries.cut(part.query, part.from, part.to))
            .collect();
        // A query whose letter before its `-query_loc` interval is masked in the batch's masks
        // restricts that letter too: a lower-case letter (`lower_case_before`), or with
        // `-filtering_db` one that the batch's repeat masks joined (`masks_outside`, bit 0).
        let restricted = chunk
            .queries
            .iter()
            .map(|part| {
                let (masks, offset) = (&query_masks[part.query], query_offsets[part.query]);
                let letter_before = lower_case_before.get(part.query) == Some(&true)
                    || masks_outside
                        .get(part.query)
                        .is_some_and(|outside| outside.before & 1 != 0);
                if letter_before {
                    restrict_masks_with_letter_before(masks, part, offset)
                } else {
                    restrict_masks(masks, part, offset)
                }
            })
            .collect();
        // The chunk's query source loads the registry of `-window_masker_taxid` before it
        // makes the masks (`WindowMasker::load_posts`): after the posts of the batch's set-up and
        // of the chunks before, and before every error of the chunk's search.
        if let Some(window_masker) = &args.window_masker {
            posts.append(window_masker.load_posts(false));
        }
        // With `-filtering_db` the masks that the chunk keeps before a part's letters join the
        // part's masks when it combines the part's repeat hits with them
        // (`query_split::masked_before_part`).
        let mut parts_outside: Vec<MasksOutside> = if args.repeat_filter.is_some() {
            chunk
                .queries
                .iter()
                .map(|part| {
                    let before = masks_outside
                        .get(part.query)
                        .map_or(0, |outside| outside.before);
                    MasksOutside::before_part(masked_before_part(
                        &query_masks[part.query],
                        before,
                        part,
                        query_offsets[part.query],
                    ))
                })
                .collect()
        } else {
            Vec::new()
        };
        let part_masks =
            chunk_query_masks(args, &parts, restricted, &mut parts_outside, parallel_pool)
                .map_err(|error| mask_rejection(args, error))?;
        let batch = search_query_batch::<false>(
            args,
            x_query_drop::BatchQueries::Shared(&parts),
            subjects,
            reborrow_outputs(&mut outputs),
            parallel_pool,
            batching,
            first_batch,
            BatchStage::ChunkPrelim {
                query_masks: &part_masks,
                batch_eff_searchsp,
                x_last: last_batch && index + 1 == chunks.len(),
            },
            &mut *posts,
            x_keep,
            &mut query_keep,
        )?;
        chunk_set_up |= batch.searched;
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:174-178
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:519-522
        merge_query_chunk(
            &mut merged,
            batch.chunk_prelim_lists,
            chunk,
            context_offsets,
            batching.split.overlap as u32 as i32,
        )?;
    }
    Ok(Some((merged, chunk_set_up)))
}

/// NCBI's `BlastHSPStreamMerge` of a query chunk: the HSPs of each query part's hit list
/// move to the part's contexts in the batch, and the list merges with the query's hit list of
/// the chunks before (`merge_prelim_hit_list`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:476-522
/// The overlap is 100 and gaps are allowed (`CSplitQueryBlk` of a gapped search).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:527-534
fn merge_query_chunk(
    merged: &mut PrelimHitLists,
    chunk_lists: PrelimHitLists,
    chunk: &QueryChunk,
    context_offsets: &[i32],
    chunk_overlap_size: i32,
) -> Result<()> {
    for (local_query, hit_list) in chunk_lists.into_iter().enumerate() {
        let Some(mut hit_list) = hit_list else {
            continue;
        };
        let part = &chunk.queries[local_query];
        for list in &mut hit_list.hsplist_array {
            for hsp in &mut list.hsps {
                let local_context = hsp.context_idx as usize;
                // The chunk's HSPs are in the contexts that it searches.
                let context = chunk.contexts[local_context].ok_or_else(|| {
                    anyhow::anyhow!(
                        "internal error: an HSP of a query chunk is in a context that is not searched"
                    )
                })?;
                // The offsets are the starts of the parts (`query_split::context_offsets`). NCBI
                // adds them unchecked: an HSP that NCBI keeps outside its chunk's context (an
                // ungapped walk past two sentinels with `-xdrop_ungap 1e10`, the megablast greedy
                // alignment from it) keeps its place relative to the part, inside or outside the
                // query's context. The sums lie within the query block, far inside an Int4.
                let offset = chunk.context_offsets[local_context];
                let shift = |value: i32| {
                    value.checked_add(offset).ok_or_else(|| {
                        anyhow::anyhow!(
                            "internal error: an HSP offset of a query chunk leaves the 32-bit range"
                        )
                    })
                };
                hsp.prelim_qs = shift(hsp.prelim_qs)?;
                hsp.prelim_qe = shift(hsp.prelim_qe)?;
                hsp.seed_qs = shift(hsp.seed_qs)?;
                hsp.context_idx = context as u32;
                hsp.query_idx = part.query as u32;
                hsp.query_frame = if context % 2 == 0 { 1 } else { -1 };
                hsp.query_context_offset = context_offsets[context];
            }
        }
        // A context that the chunk does not have (a strand that is not searched) keeps -1.
        let split_point = |strand: usize| {
            let position = 2 * local_query + strand;
            if chunk.contexts[position].is_some() {
                chunk.context_offsets[position]
            } else {
                -1
            }
        };
        let split_points = [split_point(0), split_point(1)];
        merge_prelim_hit_list(
            hit_list,
            &mut merged[part.query],
            split_points,
            chunk_overlap_size,
        );
    }
    for hit_list in merged.iter_mut().flatten() {
        for list in &mut hit_list.hsplist_array {
            sort_prelim_hits_by_score(&mut list.hsps);
        }
    }
    Ok(())
}

/// NCBI's `Blast_HitListMerge`: a query chunk's hit list of a query merges with the query's
/// hit list of the chunks before into a new hit list of the same size, subject by subject in
/// subject order (`Blast_HitListUpdate` keeps the best `prelim_hitlist_size` subjects by
/// their preliminary e-values). The HSP lists of a subject in both merge with
/// `Blast_HSPListsMerge` when a split point is positive, else `Blast_HSPListAppend`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2132-2217
fn merge_prelim_hit_list(
    mut hitlist1: HitList<PrelimHspList>,
    combined: &mut Option<HitList<PrelimHspList>>,
    split_offsets: [i32; 2],
    chunk_overlap_size: i32,
) {
    let Some(mut hitlist2) = combined.take() else {
        *combined = Some(hitlist1);
        return;
    };
    let mut new_hitlist = HitList::new(hitlist1.hsplist_max);
    // The subjects of a hit list are distinct, so the sort order is total.
    hitlist1.hsplist_array.sort_by_key(|list| list.oid);
    hitlist2.hsplist_array.sort_by_key(|list| list.oid);
    let query_is_split = split_offsets.iter().any(|&offset| offset > 0);
    let mut lists1 = std::mem::take(&mut hitlist1.hsplist_array)
        .into_iter()
        .peekable();
    let mut lists2 = std::mem::take(&mut hitlist2.hsplist_array)
        .into_iter()
        .peekable();
    while let (Some(list1), Some(list2)) = (lists1.peek(), lists2.peek()) {
        if list1.oid < list2.oid {
            new_hitlist.update(lists1.next().unwrap());
        } else if list1.oid > list2.oid {
            new_hitlist.update(lists2.next().unwrap());
        } else {
            let list1 = lists1.next().unwrap();
            let mut list2 = lists2.next().unwrap();
            if query_is_split {
                merge_prelim_hit_lists(
                    &mut list2.hsps,
                    list1.hsps,
                    HspListSplit::Query {
                        offsets: split_offsets,
                    },
                    chunk_overlap_size as usize,
                    true,
                );
            } else {
                append_prelim_hit_list(&mut list2.hsps, list1.hsps);
            }
            new_hitlist.update(list2);
        }
    }
    for list in lists1.chain(lists2) {
        new_hitlist.update(list);
    }
    *combined = Some(new_hitlist);
}

#[cfg(test)]
#[path = "run_tests.rs"]
mod tests;
