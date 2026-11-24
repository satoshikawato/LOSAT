//! BLASTP search coordination and execution
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c

use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::{bail, Context, Result};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::prelude::*;
use std::io::Write;
use std::path::Path;
use std::sync::{
    atomic::{AtomicU64, Ordering as AtomicOrdering},
    Arc,
};
use std::time::Instant;

use crate::algorithm::blastn::blast_engine::QueryReading;
use crate::algorithm::blastn::interval_tree::{BlastIntervalTree, IndexMethod, TreeHsp};
use crate::algorithm::common::db_mask_search::SearchMask;
use crate::algorithm::tblastx::blast_aascan::{
    s_blast_aa_scan_subject_one_range, BlastOffsetPair as OffsetPair,
};
use crate::algorithm::tblastx::blast_extend::DiagStruct;
use crate::algorithm::tblastx::lookup::compressed::x_clut_direct::{
    x_bp_build_compressed_lookup, x_bp_build_compressed_lookup_of, x_bpclut_on, XCompressedLookup,
    XFlatCompressedLookup, X_BLASTP, X_BLASTP_MATRIX,
};
use crate::algorithm::tblastx::lookup::compressed::BlastCompressedAaLookupTable;
use crate::algorithm::tblastx::lookup::{
    x_bp_build_ncbi_lookup, x_bp_build_ncbi_lookup_as, BackboneCell16, BlastAaLookupTable,
    QueryContext, XBpLookup,
};
use crate::algorithm::tblastx::ncbi_cutoffs::{gap_trigger_raw_score, x_drop_raw_score};
use crate::algorithm::tblastx::translation::QueryFrame;
use crate::api::local_blast::{FormatObserver, FormatProbe, HspIndex, OutputSink, ReportOutputs};
use crate::api::QueryEnd;
use crate::blastinput::fasta_reader::FastaRecord;
use crate::blastinput::seq_range;
use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::core::composition_adjustment::adjust_scores::BlastCompositionWorkspace;
use crate::core::composition_adjustment::redo_alignment::BlastCompoAdjustMode;
use crate::report::diag_log::Posts;
use crate::report::row_letters::RowSide;
use crate::report::{
    write_blastp_pairwise_report, BlastpPairwiseReport, OutputConfig, OutputFormat, PairwiseConfig,
    PairwiseHit, ReportContext,
};
use crate::stats::compute_blosum62_ideal_karlin_params;
use crate::stats::search_space::SearchSpace;
use crate::stats::{
    blast_spouge_etos, blast_spouge_stoe, lookup_protein_gumbel_params,
    lookup_protein_params_ungapped, BlastGumbelBlk, KarlinParams,
};
use crate::utils::matrix::standard_score;

// `-best_hit_overhang` / `-best_hit_score_edge` (series NP, O6a-2): the best hit writer and pipe.
mod besthit;
mod culling;
mod database;
// `-db_soft_mask` and `-db_hard_mask` (series NP, M4b): the kernels of a masked database search,
// chosen once per call in place of the kernel (`prelim_run::run_preliminary_with`).
mod db_masked;
// What the replay of the parallel redo's windows (`x_redo_par.rs`, built where they run) reads of
// the one-thread loop's Smith-Waterman redo: the marks of the lists that the redo's fetch excludes
// and the loop's reader of a query's marks, the errors of a batch and of a match (with the refusal
// or the warning of a header read before them), and the notes they read (the status of a query
// batch is core's `smith_waterman_match::SwBatchStatus`, which both import).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) use {
    database::{FetchExcludedMatch, FetchNotes},
    redo_run::{query_fetch_excluded_reached, sw_batch_failed_after, sw_endless_after},
};
// The errors that the windows' tests expect (`x_redo_par_tests.rs`).
#[cfg(all(
    test,
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
pub(super) use redo_run::{sw_batch_failed, sw_endless};
// The kernel of one subject (`BlastpPrelimKernel::search_subject`): its scan and ungapped
// extensions, then its gapped half, from the subject's letters to its HSP lists.
mod prelim_kernel;
// The preliminary stage of a call after its kernel (`run_call`): the choice of the pool and the
// call's branch (RPS-BLAST, PSSM, query blocks, a database's blocks, a subject set); the close of
// its stream and the sort of its preliminary lists.
mod prelim_run;
// EXPERIMENT (LOSAT_X_BPIHPOOL): the arrays of initial hits past NCBI's first size in one pool
// shared by the threads.
mod x_ih_pool;
// EXPERIMENT (LOSAT_X_BPHITSFLAT): a subject's preliminary HSPs grouped by query through one array
// of their positions.
mod x_hits_flat;
// EXPERIMENT (LOSAT_X_BPSUBJRUN): the parallel `-subject` loop in runs of subjects, one diagonal
// table fill per run.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
mod x_subject_runs;
// EXPERIMENT (LOSAT_X_BPQBLOCK): the subjects' tabular IDs and titles of the block search made
// once for all blocks (held by `x_query_blocks::XBlockSubjects`).
pub(super) mod x_subject_meta;
// `-comp_based_stats 0`: NCBI's ordinary traceback in place of the composition-based redo; the
// steps after any traceback on the final hit lists.
mod traceback;
// The composition-based redo: its parameters of each query, its plan and its paths (`run_call`).
mod redo_run;
// psiblast `-in_pssm` (series PS, item PS4a): the stages of a PSSM search, reached through the
// `psi_search` hooks of `run_call` and its preliminary stage (`blast_engine/prelim_run.rs`).
mod pssm_hooks;
// rpsblast (series PR, item PS11a): the RPS engine's search, reached through the `rps_search`
// hooks of `search_cli`, `run_call` and its preliminary stage (`blast_engine/prelim_run.rs`).
mod rps_hooks;
// psiblast's rounds of FASTA queries (series PS, item PS5a): the round driver, reached from
// `search`, and the round's outcome hook of `run_call`.
mod psi_round;
// The rendered rows of the final hits a window of queries at a time, as the reports read them.
pub(super) mod query_rows;
// The setup of a call before its search: its requested reports, query frames, subjects, scoring
// values, lookup table, search spaces, cutoffs, hit lists and scan's sizes (`run_call`).
mod call_setup;
// The warning lines of a call's queries, the checks and the prefetch of its reported subjects and
// its refusals before any format is written (`run_call`).
mod report_prep;
// The reports of a call's final hits, one function per arm of the loop over the formats
// (`run_call`).
mod report_write;
// An ungapped HSP of length 0 in the final hit lists: outfmt 1 to 4's refusal and `-subject`'s
// exception (`run_call`).
mod empty_alignments;

use super::alignment::build_alignment_view_with_matrix;
use super::args::{
    BlastpApp, BlastpArgs, BlastpCompositionMode, BlastpLookupTableType, ResolvedBlastpArgs,
};
use super::encoding::{
    encode_protein_query_frame_lcase, encode_protein_query_frame_soft_masked,
    encode_protein_query_frame_with_seg, encode_protein_query_frame_with_seg_timed,
    encode_protein_sequence, ncbistdaa_to_ascii, EncodedProtein,
};
use super::extension::BlastpUngappedData;
#[cfg(test)]
use super::extension::{extend_one_hit, extend_two_hit};
use super::gapalign::{
    blastp_get_start_for_gapped_alignment, blastp_score_only_gapped_alignment_with_scratch,
    BlastpGappedAlignmentMode, BlastpPreliminaryHsp, GapAlignScratch,
};
use super::hsp::{
    collect_hits_and_methods_from_hit_lists, fill_results_from_compo_heaps,
    hit_list_to_seqalign_order, reap_hsplist_by_evalue, BlastCompoHeap, BlastpHitList, BlastpHsp,
    BlastpHspList,
};
use super::kappa::{
    blastp_kappa_range_counters_env_enabled, print_blastp_kappa_range_counters,
    reset_blastp_kappa_range_counters, BlastpTraceHspTarget,
};
use super::protein_query_split::protein_query_split_sizes;
use super::x_prelim_compact::{x_expanded_lists, x_finish_compact_lists, XPrelimHsp};
// No NCBI counterpart: matches per Newton lanes prefetch in the BLASTP redo loops; it does not
// change any value NCBI computes.
const X_NEWTON_LANES_MATCHES: usize = 16;
// No NCBI counterpart: the size limit of the preliminary interval tree's grid index under
// LOSAT_X_BPTREECELLS (the tree's own limit is 1 << 16 cells); it does not change any value NCBI
// computes.
const X_BP_TREE_MAX_CELLS: usize = 1 << 14;
// No NCBI counterpart: reads the LOSAT_X_BPTREECELLS switch once; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPTREECELLS): strict set (on unless `0`). The grid cells of the
/// preliminary interval tree at most, per subject.
fn x_bp_tree_max_cells() -> usize {
    static CELLS: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *CELLS.get_or_init(|| {
        if crate::utils::xswitch::on("LOSAT_X_BPTREECELLS") {
            X_BP_TREE_MAX_CELLS
        } else {
            1 << 16
        }
    })
}
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_extend.h:142-155
#[derive(Clone, Copy, Debug)]
struct InitHSP {
    // NCBI BlastInitHSP::offsets.qs_offsets.{q_off,s_off}.
    q_seed_absolute: i32,
    s_seed: i32,
    // NCBI BlastInitHSP::ungapped_data payload.
    ungapped_data: BlastpUngappedData,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3905
impl InitHSP {
    #[inline]
    fn query_end_absolute(self) -> i32 {
        self.ungapped_data.q_start + self.ungapped_data.length
    }

    #[inline]
    fn subject_end(self) -> i32 {
        self.ungapped_data.s_start + self.ungapped_data.length
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1498
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1688
struct BlastpSubjectResult {
    subject_index: usize,
    // The lists with the compact element (`x_prelim_compact.rs`).
    x_compact_hsp_lists: Vec<(usize, BlastpHspList<XPrelimHsp>)>,
}

impl BlastpSubjectResult {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
    fn new(subject_index: usize) -> Self {
        Self {
            subject_index,
            x_compact_hsp_lists: Vec::new(),
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1554
    fn reset(&mut self, subject_index: usize) {
        self.subject_index = subject_index;
        self.x_compact_hsp_lists.clear();
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:991-1041
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:109-146
/// The statistics of one subject that NCBI's preliminary search computes for a database length of
/// 0 (`BLAST_OneSubjectUpdateParameters`, `call_setup::subject_parameters`): each context's search
/// space (the preliminary e-values without the Gumbel block and the linking) and the cutoffs at
/// the subject's length (`blastp_prelim_cutoffs`: hit-saving and word cutoffs, and with psiblast's
/// `-sum_stats true` the subject's sum statistics, the lowered hit-saving cutoffs and the
/// linking's small-gap cutoff).
struct BlastpSubjectParameters {
    search_spaces: Vec<SearchSpace>,
    cutoffs: BlastpPrelimCutoffs,
}

impl BlastpSubjectParameters {
    /// The values the kernel reads of these statistics; the preliminary e-value of a subject's own
    /// statistics is the options' (`-evalue`, `run_call`'s `prelim_evalue` without the
    /// composition-based statistics' stretch).
    fn values(&self, evalue: f64) -> BlastpPrelimValues<'_> {
        BlastpPrelimValues {
            word_cutoff_scores: &self.cutoffs.word_cutoff_scores,
            hit_cutoff_scores: &self.cutoffs.hit_cutoff_scores,
            search_spaces: &self.search_spaces,
            prelim_evalue: evalue,
            sum_stats: self.cutoffs.sum_stats.as_ref(),
        }
    }
}

/// What `run_call`'s kernel (`BlastpPrelimKernel::search_subject`, `blast_engine/prelim_kernel.rs`)
/// reads of the preliminary statistics: the word and
/// hit-saving cutoffs of each context, the search spaces, the preliminary e-value and psiblast's
/// sum statistics; the setup's (made once per call) or a subject's own
/// (`BlastpSubjectParameters::values`). One pointer per subject, so that the default path reads
/// the setup's as it did before the statistics of each subject.
struct BlastpPrelimValues<'a> {
    word_cutoff_scores: &'a [i32],
    hit_cutoff_scores: &'a [i32],
    search_spaces: &'a [SearchSpace],
    prelim_evalue: f64,
    sum_stats: Option<&'a BlastpSumStats>,
}

struct BlastpSubjectScratch {
    offset_pairs: Vec<OffsetPair>,
    diag_array: Vec<DiagStruct>,
    diag_offset: i32,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:253-267
    init_hsps: Vec<InitHSP>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3524-3527
    chaining_nodes: Vec<BlastpInitHspNode>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3628-3633
    chaining_keep: Vec<bool>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3971-4008
    preliminary_hsp_list: Vec<BlastpPreliminaryHsp>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3997
    //
    // Rust tracks how many accepted preliminary HSPs belong to each query so
    // the redo path can tell when NCBI's remove/rebuild loop would leave the
    // interval tree unchanged.
    preliminary_hsp_counts_by_query: Vec<usize>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:117-149
    preliminary_hits_by_query: Vec<Vec<BlastpPreliminaryHsp>>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:112
    // EXPERIMENT (LOSAT_X_BPHITSFLAT): in place of `preliminary_hits_by_query`, the positions of
    // the HSPs in query order and each query's start (`x_hits_flat.rs`).
    x_hits_order: Vec<u32>,
    x_hits_starts: Vec<u32>,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3748-3786
    restricted_align_found: Vec<bool>,
    restricted_align_array: Vec<bool>,
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:69-80
    preliminary_gap_scratch: GapAlignScratch,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
    preliminary_interval_tree: BlastIntervalTree,
    // `-db_soft_mask` and `-db_hard_mask`: the masked kernels' reader, plan and buffers
    // (`db_masked.rs`), made at a worker's first masked sequence; a search without them keeps
    // `None` (one word in the scratch).
    db_mask: Option<Box<db_masked::BlastpMaskScratch>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1681
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1554
struct BlastpTiming {
    query_encoding_ns: AtomicU64,
    subject_encoding_ns: AtomicU64,
    lookup_ns: AtomicU64,
    scan_ungapped_ns: AtomicU64,
    scan_ungapped_calls: AtomicU64,
    preliminary_gapped_ns: AtomicU64,
    preliminary_gapped_calls: AtomicU64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3955-4009
    preliminary_exact_retry_count: AtomicU64,
    preliminary_retry_removed_hsps: AtomicU64,
    preliminary_interval_tree_rebuilds: AtomicU64,
    preliminary_interval_tree_rebuild_nodes: AtomicU64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2525
    preliminary_endpoint_duplicates_purged: AtomicU64,
    prelim_merge_sort_ns: AtomicU64,
    prelim_merge_sort_calls: AtomicU64,
    kappa_redo_ns: AtomicU64,
    kappa_redo_calls: AtomicU64,
    final_traceback_ns: AtomicU64,
    final_traceback_calls: AtomicU64,
    output_format_ns: AtomicU64,
    output_format_calls: AtomicU64,
}

impl BlastpTiming {
    fn new() -> Self {
        Self {
            query_encoding_ns: AtomicU64::new(0),
            subject_encoding_ns: AtomicU64::new(0),
            lookup_ns: AtomicU64::new(0),
            scan_ungapped_ns: AtomicU64::new(0),
            scan_ungapped_calls: AtomicU64::new(0),
            preliminary_gapped_ns: AtomicU64::new(0),
            preliminary_gapped_calls: AtomicU64::new(0),
            preliminary_exact_retry_count: AtomicU64::new(0),
            preliminary_retry_removed_hsps: AtomicU64::new(0),
            preliminary_interval_tree_rebuilds: AtomicU64::new(0),
            preliminary_interval_tree_rebuild_nodes: AtomicU64::new(0),
            preliminary_endpoint_duplicates_purged: AtomicU64::new(0),
            prelim_merge_sort_ns: AtomicU64::new(0),
            prelim_merge_sort_calls: AtomicU64::new(0),
            kappa_redo_ns: AtomicU64::new(0),
            kappa_redo_calls: AtomicU64::new(0),
            final_traceback_ns: AtomicU64::new(0),
            final_traceback_calls: AtomicU64::new(0),
            output_format_ns: AtomicU64::new(0),
            output_format_calls: AtomicU64::new(0),
        }
    }

    #[inline]
    fn record_duration(counter: &AtomicU64, start: Option<Instant>) {
        if let Some(start) = start {
            let elapsed_ns = start.elapsed().as_nanos().min(u64::MAX as u128) as u64;
            counter.fetch_add(elapsed_ns, AtomicOrdering::Relaxed);
        }
    }

    #[inline]
    fn record_call(counter: &AtomicU64, calls: &AtomicU64, start: Option<Instant>) {
        Self::record_duration(counter, start);
        calls.fetch_add(1, AtomicOrdering::Relaxed);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1705
    #[inline]
    fn record_count(counter: &AtomicU64, value: usize) {
        counter.fetch_add(value.min(u64::MAX as usize) as u64, AtomicOrdering::Relaxed);
    }

    #[inline]
    fn seconds(counter: &AtomicU64) -> f64 {
        counter.load(AtomicOrdering::Relaxed) as f64 / 1e9
    }

    #[inline]
    fn calls(counter: &AtomicU64) -> u64 {
        counter.load(AtomicOrdering::Relaxed)
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1470-1475
#[inline]
fn blastp_timing_start(enabled: bool) -> Option<Instant> {
    if enabled {
        Some(Instant::now())
    } else {
        None
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1674-1681
#[inline]
fn blastp_timing_env_enabled() -> bool {
    #[cfg(target_arch = "wasm32")]
    {
        false
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        std::env::var_os("LOSAT_TIMING").is_some()
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1705
fn print_blastp_timing(timing: &BlastpTiming, total_start: Option<Instant>) {
    eprintln!(
        "[TIMING] blastp_query_encoding_seg: {:.3}s",
        BlastpTiming::seconds(&timing.query_encoding_ns)
    );
    eprintln!(
        "[TIMING] blastp_subject_encoding: {:.3}s",
        BlastpTiming::seconds(&timing.subject_encoding_ns)
    );
    eprintln!(
        "[TIMING] blastp_build_lookup: {:.3}s",
        BlastpTiming::seconds(&timing.lookup_ns)
    );
    eprintln!(
        "[TIMING] blastp_scan_ungapped: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.scan_ungapped_ns),
        BlastpTiming::calls(&timing.scan_ungapped_calls)
    );
    eprintln!(
        "[TIMING] blastp_prelim_gapped: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.preliminary_gapped_ns),
        BlastpTiming::calls(&timing.preliminary_gapped_calls)
    );
    eprintln!(
        "[TIMING] blastp_prelim_exact_retries: {} (removed_hsps={}, interval_tree_rebuilds={}, rebuild_nodes={}, endpoint_duplicates_purged={})",
        BlastpTiming::calls(&timing.preliminary_exact_retry_count),
        BlastpTiming::calls(&timing.preliminary_retry_removed_hsps),
        BlastpTiming::calls(&timing.preliminary_interval_tree_rebuilds),
        BlastpTiming::calls(&timing.preliminary_interval_tree_rebuild_nodes),
        BlastpTiming::calls(&timing.preliminary_endpoint_duplicates_purged)
    );
    eprintln!(
        "[TIMING] blastp_prelim_merge_sort: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.prelim_merge_sort_ns),
        BlastpTiming::calls(&timing.prelim_merge_sort_calls)
    );
    eprintln!(
        "[TIMING] blastp_kappa_redo: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.kappa_redo_ns),
        BlastpTiming::calls(&timing.kappa_redo_calls)
    );
    eprintln!(
        "[TIMING] blastp_final_traceback_stats: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.final_traceback_ns),
        BlastpTiming::calls(&timing.final_traceback_calls)
    );
    eprintln!(
        "[TIMING] blastp_output_format: {:.3}s (calls={})",
        BlastpTiming::seconds(&timing.output_format_ns),
        BlastpTiming::calls(&timing.output_format_calls)
    );
    // EXPERIMENT (LOSAT_X_BPDIAG): the cells the fused loop tested in batches (`x_seed_fused.rs`).
    super::x_seed_fused::x_bp_diag_print_timing();
    if let Some(total_start) = total_start {
        eprintln!(
            "[TIMING] blastp_total: {:.3}s",
            total_start.elapsed().as_secs_f64()
        );
    }
    // The redone lists kept with fitted room (`x_redo_fit.rs`).
    super::x_redo_fit::x_bp_redo_fit_print_timing();
}

impl BlastpSubjectScratch {
    fn new(offset_array_size: i32, diag_array_size: i32, window: i32, query_span: i32) -> Self {
        Self {
            offset_pairs: vec![
                OffsetPair::default();
                usize::try_from(offset_array_size)
                    .expect("NCBI BLAST offset pair array size must fit in usize")
            ],
            diag_array: vec![
                DiagStruct::default();
                usize::try_from(diag_array_size)
                    .expect("NCBI BLAST diagonal array size must fit in usize")
            ],
            diag_offset: window,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:253-267
            init_hsps: x_ih_pool::own_array(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3524-3527
            chaining_nodes: Vec::new(),
            chaining_keep: Vec::new(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3971-3975
            preliminary_hsp_list: Vec::new(),
            preliminary_hsp_counts_by_query: Vec::new(),
            preliminary_hits_by_query: Vec::new(),
            x_hits_order: Vec::new(),
            x_hits_starts: Vec::new(),
            restricted_align_found: Vec::new(),
            restricted_align_array: Vec::new(),
            preliminary_gap_scratch: GapAlignScratch::new(),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3830-3831
            preliminary_interval_tree: BlastIntervalTree::new(0, query_span, 0, 1),
            db_mask: None,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:109-146
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn prepare_independent_subject(&mut self, diag_offset: i32) {
        self.diag_offset = diag_offset;
        self.diag_array.fill(DiagStruct::default());
        self.reset_subject_state();
    }

    // No NCBI counterpart as a function: the per-subject resets of `prepare_independent_subject`
    // other than the diagonal table's, split out of it unchanged (each with its NCBI reference
    // below); EXPERIMENT (LOSAT_X_BPSUBJRUN) runs them alone before each later subject of a run
    // (`blast_engine/x_subject_runs.rs`). It changes no value NCBI computes.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    fn reset_subject_state(&mut self) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:229-235
        self.init_hsps.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3524-3578
        self.chaining_nodes.clear();
        self.chaining_keep.clear();
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3971-3975
        self.preliminary_hsp_list.clear();
        self.preliminary_hsp_counts_by_query.clear();
        for query_hits in &mut self.preliminary_hits_by_query {
            query_hits.clear();
        }
        self.restricted_align_found.clear();
        self.restricted_align_array.clear();
        self.preliminary_interval_tree.reset();
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:152-172
    fn finish_subject(&mut self, subject_length: usize, window: i32) {
        if self.diag_offset >= i32::MAX / 4 {
            self.diag_offset = window;
            for diag in &mut self.diag_array {
                *diag = DiagStruct::clear(window);
            }
        } else {
            self.diag_offset += subject_length as i32 + window;
        }
    }
}

#[inline]
fn path_label(path: &Path) -> String {
    path.file_name()
        .map(|name| name.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.display().to_string())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
struct BlastpPreparedSubject {
    encoded: EncodedProtein,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1659-1668
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn prepare_blastp_subjects_preserving_order(
    subject_records: &[FastaRecord],
    requested_threads: usize,
    effective_threads: usize,
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<Vec<BlastpPreparedSubject>> {
    let prepare_subject = |record: &FastaRecord| BlastpPreparedSubject {
        encoded: encode_protein_sequence(record.seq()),
    };

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
    if blastp_should_parallelize_subject_prep(subject_records, requested_threads, effective_threads)
    {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1688
        let pool = blastp_parallel_pool;
        return Ok(pool.install(|| subject_records.par_iter().map(prepare_subject).collect()));
    }

    Ok(subject_records.iter().map(prepare_subject).collect())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
#[cfg(any(
    not(feature = "parallel"),
    all(target_arch = "wasm32", not(feature = "wasm-threads"))
))]
fn prepare_blastp_subjects_preserving_order(
    subject_records: &[FastaRecord],
    _num_threads: usize,
) -> Result<Vec<BlastpPreparedSubject>> {
    let prepare_subject = |record: &FastaRecord| BlastpPreparedSubject {
        encoded: encode_protein_sequence(record.seq()),
    };

    Ok(subject_records.iter().map(prepare_subject).collect())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3698-3699
const BLASTP_RESTRICTED_MULT: f64 = 0.68;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:672-676
const BLASTP_KAPPA_SCALING_FACTOR: f64 = 32.0;
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2397-2399,2402-2404
const BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION: f64 = 1.74;
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_options.h:163
/// The inclusion threshold of the options of every program (`PSIBlastOptionsNew`), which
/// psiblast's `-inclusion_ethresh` changes (`ResolvedBlastpArgs::inclusion_ethresh`).
pub(super) const PSI_INCLUSION_ETHRESH: f64 = 0.002;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2380-2382
// `final_bits` is `-xdrop_gap_final` or BLASTP's default (`ResolvedBlastpArgs`); the conversion is
// the protein programs' (`core::x_drop::redo_x_dropoff`).
fn blastp_redo_x_dropoff(
    final_bits: f64,
    scaled_gapped_lambda: f64,
    x_drop_gapped_final: i32,
) -> i32 {
    crate::core::x_drop::redo_x_dropoff(final_bits, scaled_gapped_lambda, x_drop_gapped_final)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2438-2448
fn blastp_redo_near_identical_cutoff(scaled_gapped_lambda: f64) -> f64 {
    (BLASTP_NEAR_IDENTICAL_BITS_PER_POSITION * std::f64::consts::LN_2) / scaled_gapped_lambda
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2450-2462
fn blastp_redo_cutoff_score(
    do_link_hsps: bool,
    cutoff_score: i32,
    local_scaling_factor: f64,
) -> i32 {
    if do_link_hsps {
        ((cutoff_score as f64) * local_scaling_factor) as i32
    } else {
        1
    }
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:187-193
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:844-855
//
// The program is blastp's only: psiblast's first round of a FASTA query runs this engine, but
// its program is `eBlastTypePsiBlast` and its gapped alignments are exact.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_options.cpp:50
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:739-741
fn blastp_restricted_align_enabled(app: BlastpApp, expect_value: f64) -> bool {
    app == BlastpApp::Blastp && expect_value <= 10.0
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-941
fn blastp_prelim_evalue(expect_value: f64, comp_based_stats: BlastpCompositionMode) -> f64 {
    let cbs_stretch = match comp_based_stats {
        BlastpCompositionMode::CompositionMatrixAdjust
        | BlastpCompositionMode::ForceFullMatrixAdjust => 5.0,
        _ => 1.0,
    };
    cbs_stretch * expect_value
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:932-941
fn blastp_cutoff_score_max(
    expect_value: f64,
    query_length: i32,
    subject_length_for_hit_cutoff: i32,
    comp_based_stats: BlastpCompositionMode,
    gapped_params: &crate::stats::KarlinParams,
    gumbel_params: Option<&crate::stats::BlastGumbelBlk>,
    search_space: i64,
) -> i32 {
    match gumbel_params {
        Some(gumbel_params) => blast_spouge_etos(
            blastp_prelim_evalue(expect_value, comp_based_stats),
            gapped_params,
            gumbel_params,
            query_length,
            subject_length_for_hit_cutoff,
        ),
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:942-944
        None => crate::algorithm::tblastx::ncbi_cutoffs::cutoff_score_from_evalue(
            expect_value,
            search_space,
            gapped_params,
        )
        .max(1),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:825-886
fn blastp_compo_adjust_mode(mode: BlastpCompositionMode) -> BlastCompoAdjustMode {
    match mode {
        BlastpCompositionMode::NoCompositionBasedStats => {
            BlastCompoAdjustMode::NoCompositionBasedStats
        }
        BlastpCompositionMode::CompositionBasedStats => BlastCompoAdjustMode::CompositionBasedStats,
        BlastpCompositionMode::CompositionMatrixAdjust => {
            BlastCompoAdjustMode::CompositionMatrixAdjust
        }
        BlastpCompositionMode::ForceFullMatrixAdjust => BlastCompoAdjustMode::ForceFullMatrixAdjust,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3099-3107
// No NCBI counterpart: code placement only; it changes no value NCBI computes.
/// Runs `f` in a function of its own: each closure gives an instance that is not inlined into the
/// caller (the per-subject kernel's branches of another matrix than BLOSUM62, item SC2b).
#[inline(never)]
fn x_outlined<R>(f: impl FnOnce() -> R) -> R {
    f()
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
/// The composition-based redo's workspace of another matrix than BLOSUM62 (item SC2b), put in
/// place of `workspace` (BLOSUM62's, which its callers make as before): cold and in place, so that
/// the callers' frames do not grow by a second workspace (SC1's `use_matrix_workspace` of
/// TBLASTN). Every matrix whose search runs the redo has NCBI's data; IDENTITY keeps BLOSUM62's
/// workspace unused (`blastp_redo_matrix`).
#[cold]
#[inline(never)]
pub(super) fn blastp_use_matrix_workspace(
    workspace: &mut BlastCompositionWorkspace,
    matrix: ScoringMatrix,
) {
    if blastp_redo_matrix(matrix) == ScoringMatrix::Blosum62 {
        return;
    }
    *workspace = BlastCompositionWorkspace::new(matrix)
        .expect("BLASTP searches only the matrices with NCBI's composition data");
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1499
/// The matrix whose composition data the redo's workspace and initial matrix take: the search's,
/// but BLOSUM62 for IDENTITY, which has none. IDENTITY is searched in mode 0 only
/// (`ResolvedBlastpArgs::identity_cbs_reset`), whose redo has no input (the preliminary lists are
/// taken before it, item SC4a), so those tables stay unused, as NCBI, which does not run the redo
/// below mode 1, makes none.
fn blastp_redo_matrix(matrix: ScoringMatrix) -> ScoringMatrix {
    if matrix == ScoringMatrix::Identity {
        ScoringMatrix::Blosum62
    } else {
        matrix
    }
}

fn blastp_local_scaling_factor(mode: BlastpCompositionMode, matrix: ScoringMatrix) -> f64 {
    if mode == BlastpCompositionMode::NoCompositionBasedStats {
        return 1.0;
    }
    let _ = matrix;
    BLASTP_KAPPA_SCALING_FACTOR
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3119-3128
fn build_redone_match_heaps(
    num_queries: usize,
    heap: impl Fn() -> BlastCompoHeap,
) -> Vec<BlastCompoHeap> {
    // `heap`: one query's heap as BlastCompo_HeapInitialize leaves it
    // (`ResolvedBlastpArgs::redone_match_heap`), built anew for each query.
    (0..num_queries).map(|_| heap()).collect()
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1475
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn blastp_should_parallelize_stage(
    stage: &str,
    _requested_threads: usize,
    effective_threads: usize,
    _subject_records: usize,
    worker_jobs: usize,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1475
    // Iterate subject OIDs with `BlastSeqSrcIteratorNext` until `BLAST_SEQSRC_EOF`.
    // This stage selects existing independent jobs within the verified pool.
    let parallel = effective_threads > 1 && worker_jobs > 1;
    crate::utils::threading::report_stage("blastp", stage, worker_jobs, parallel);
    parallel
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn blastp_should_parallelize_subject_prep(
    subject_records: &[FastaRecord],
    requested_threads: usize,
    effective_threads: usize,
) -> bool {
    blastp_should_parallelize_stage(
        "subject-prep",
        requested_threads,
        effective_threads,
        subject_records.len(),
        subject_records.len(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:152-172
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn blastp_subject_diag_offsets(subjects: &[EncodedProtein], window: i32) -> Vec<i32> {
    let mut offsets = Vec::with_capacity(subjects.len());
    let mut diag_offset = window;
    for subject in subjects {
        offsets.push(diag_offset);
        // A subject without letters is not searched (blast_engine.c:1429-1432,
        // `process_subject`), so the offset does not move for it.
        if subject.aa_len == 0 {
            continue;
        }
        if diag_offset >= i32::MAX / 4 {
            diag_offset = window;
        } else {
            diag_offset += subject.aa_len as i32 + window;
        }
    }
    offsets
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3429-3449
// The matches of query `q_idx` (`hsplist_count` of its preliminary hit list), from the compact
// list while it is not expanded yet (`x_prelim_compact.rs`).
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn x_prelim_hsplist_count(
    preliminary_hit_lists: &[Option<BlastpHitList>],
    x_compact_hit_lists: &[Option<BlastpHitList<XPrelimHsp>>],
    q_idx: usize,
) -> usize {
    match (&preliminary_hit_lists[q_idx], &x_compact_hit_lists[q_idx]) {
        (Some(hit_list), _) => hit_list.hsplist_count,
        (None, Some(hit_list)) => hit_list.hsplist_count,
        (None, None) => 0,
    }
}

// No NCBI counterpart: NCBI's `BlastHSPStreamWrite` writes into the HSP stream's writer, a
// collector, a culling writer or a best hit writer (`CreateHspWriter`). The database search takes
// the writer as a type, so the default search is the collector's instance of
// `database::search_blocks` and a filtering option's writer (`PrelimFilter`) is a second instance
// (the choice is made once per call, outside the loops).
/// Where a subject's preliminary HSP lists go: the order of the calls is the subject (OID) order.
trait PrelimSink {
    /// `BlastHSPStreamWrite` of one subject's lists; the subject result is emptied.
    fn add(&mut self, subject_result: &mut BlastpSubjectResult);
}

/// The default writer: `hspfilter_collector.c`, a heap of `prelim_hitlist_size` lists per query.
struct PrelimCollector<'a> {
    x_compact_hit_lists: &'a mut [Option<BlastpHitList<XPrelimHsp>>],
    prelim_hitlist_size: usize,
}

impl PrelimSink for PrelimCollector<'_> {
    #[inline(always)]
    fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        update_blastp_preliminary_hit_lists(
            self.x_compact_hit_lists,
            subject_result,
            self.prelim_hitlist_size,
        );
    }
}

/// The writer of a filtering option in place of the collector (`CreateHspWriter`,
/// setup_factory.cpp:305-342): `-culling_limit`'s culling writer or the best hit writer of
/// `-best_hit_overhang` / `-best_hit_score_edge`, which exclude each other.
enum PrelimFilter {
    Culling(culling::PrelimCulling),
    BestHit(besthit::PrelimBestHit),
}

impl PrelimFilter {
    /// The writer of the call's options, if one of them is given.
    fn of_options(args: &ResolvedBlastpArgs, prelim_hitlist_size: usize) -> Option<Self> {
        if args.culling_limit > 0 {
            Some(Self::Culling(culling::PrelimCulling::new(
                prelim_hitlist_size,
                args.culling_limit,
            )))
        } else {
            args.best_hit.map(|params| {
                Self::BestHit(besthit::PrelimBestHit::new(prelim_hitlist_size, params))
            })
        }
    }

    /// `BlastHSPStreamWrite` of one subject's lists; the subject result is emptied.
    fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        match self {
            Self::Culling(culling) => culling.add(subject_result),
            Self::BestHit(best_hit) => best_hit.add(subject_result),
        }
    }

    /// The writer's `Final`: the hit lists of the queries go into `compact_lists`.
    fn finish(self, compact_lists: &mut [Option<BlastpHitList<XPrelimHsp>>]) {
        match self {
            Self::Culling(culling) => culling.finish(compact_lists),
            Self::BestHit(best_hit) => best_hit.finish(compact_lists),
        }
    }
}

impl PrelimSink for PrelimFilter {
    fn add(&mut self, subject_result: &mut BlastpSubjectResult) {
        PrelimFilter::add(self, subject_result);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:139-161
fn update_blastp_preliminary_hit_lists(
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    subject_result: &mut BlastpSubjectResult,
    prelim_hitlist_size: usize,
) {
    let _subject_index = subject_result.subject_index;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:139-159,161
    for (q_idx, provisional_hsp_list) in subject_result.x_compact_hsp_lists.drain(..) {
        x_compact_hit_lists[q_idx]
            .get_or_insert_with(|| BlastpHitList::new(prelim_hitlist_size))
            .update(provisional_hsp_list);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:271-327
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn merge_blastp_subject_results_in_order(
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    subject_results: Vec<Option<BlastpSubjectResult>>,
    prelim_hitlist_size: usize,
) {
    for (subject_index, subject_result) in subject_results.into_iter().enumerate() {
        let Some(mut subject_result) = subject_result else {
            continue;
        };
        debug_assert_eq!(subject_result.subject_index, subject_index);
        update_blastp_preliminary_hit_lists(
            x_compact_hit_lists,
            &mut subject_result,
            prelim_hitlist_size,
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
/// `merge_blastp_subject_results_in_order` with a filtering option's writer in place of the
/// collector: the same order, the subjects by index.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
fn merge_blastp_subject_results_filtered(
    filter: &mut PrelimFilter,
    subject_results: Vec<Option<BlastpSubjectResult>>,
) {
    for (subject_index, subject_result) in subject_results.into_iter().enumerate() {
        let Some(mut subject_result) = subject_result else {
            continue;
        };
        debug_assert_eq!(subject_result.subject_index, subject_index);
        filter.add(&mut subject_result);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:1098-1101
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:500-503
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
fn build_blastp_concatenated_query(contexts: &[QueryContext]) -> Vec<u8> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
    //
    // `ctx.aa_seq[1..]` is the per-context `query->sequence` slice
    // (residues plus that context's trailing NULLB). Concatenating those slices
    // already produces the full `query->sequence` buffer; the extra calloc'd
    // byte that exists after `query->length` in NCBI is not part of
    // `query->sequence` itself and must not contribute to Rust slice length.
    let total_len = contexts
        .iter()
        .map(|ctx| ctx.aa_seq.len() - 1)
        .sum::<usize>();
    let mut concatenated = Vec::with_capacity(total_len);
    for ctx in contexts {
        concatenated.extend_from_slice(&ctx.aa_seq[1..]);
    }
    concatenated
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:246-250
//
// LOSAT's `query_length` mirrors `BLAST_SequenceBlk->length`, which excludes
// the final trailing NULLB but includes shared inter-context boundary sentinels.
fn blastp_concatenated_query_length(contexts: &[QueryContext]) -> i32 {
    contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:204-221
// `bits` is `-xdrop_ungap` or BLASTP's default (`ResolvedBlastpArgs`).
fn blastp_x_dropoff_init_for_context(ctx: &QueryContext, bits: f64) -> i32 {
    if !ctx.is_valid {
        return 0;
    }
    x_drop_raw_score(bits, &ctx.karlin_params, 1.0)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:904-940
fn blastp_hit_cutoff_score_for_context(
    ctx: &QueryContext,
    evalue: f64,
    subject_length: i32,
    mode: BlastpCompositionMode,
    gapped_params: &KarlinParams,
    gapped_gumbel: Option<&BlastGumbelBlk>,
    search_space: i64,
) -> i32 {
    if !ctx.is_valid {
        return i32::MAX;
    }
    blastp_cutoff_score_max(
        evalue,
        ctx.aa_len as i32,
        subject_length,
        mode,
        gapped_params,
        gapped_gumbel,
        search_space,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:324-345
// `bits` is psiblast's `-gap_trigger` or the protein default (`ResolvedBlastpArgs`).
fn blastp_word_cutoff_score_for_context(
    hit_cutoff_score: i32,
    ctx: &QueryContext,
    bits: f64,
) -> i32 {
    if !ctx.is_valid {
        return i32::MAX;
    }
    let gap_trigger = gap_trigger_raw_score(bits, &ctx.karlin_params);
    gap_trigger.min(hit_cutoff_score)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:346-374
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_parameters.h:77
/// The word cutoff of a valid context in an ungapped search: `BLAST_Cutoffs` for the e-value
/// `CUTOFF_E_BLASTP` (1e-300) in the search space of the shorter of the subject and the query
/// times the subject's length, with the context's own block, at most the gap trigger and the hit-saving cutoff
/// `hit_cutoff_score` (`cutoff_score_max`): `word_cutoff_score_for_update`, the function that
/// ports these lines for TBLASTX too, with BLASTP's E-value. `word_subject_length` is the sequence
/// source's average length (`BlastSeqSrcGetAvgSeqLen`), or the subject's length with the
/// statistics of each subject (`BlastInitialWordParametersUpdate`). There is no link parameter
/// set in BLASTP, so the gap decay rate is 0 (the hit parameters have no link parameters) and the
/// scale factor 1.
fn blastp_ungapped_word_cutoff_score_for_context(
    hit_cutoff_score: i32,
    ctx: &QueryContext,
    bits: f64,
    word_subject_length: i32,
) -> i32 {
    if !ctx.is_valid {
        return i32::MAX;
    }
    let gap_trigger = gap_trigger_raw_score(bits, &ctx.karlin_params);
    crate::algorithm::tblastx::ncbi_cutoffs::word_cutoff_score_for_update(
        crate::algorithm::tblastx::ncbi_cutoffs::CUTOFF_E_BLASTP,
        ctx.aa_len as i64,
        i64::from(word_subject_length),
        gap_trigger,
        hit_cutoff_score,
        0.0,
        &ctx.karlin_params,
        1.0,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:949-976
/// The sum statistics' cutoff score of a context (`BLAST_Cutoffs` with the gap decay, from an
/// e-value of 1 in the search space of the average query and subject lengths); `BLAST_Cutoffs`
/// keeps the `new_cutoff` of 1 when its score is lower (blast_stat.c:4108-4129).
pub(super) fn sum_stats_hit_cutoff(
    contexts: &[QueryContext],
    avg_subject_length: i32,
    gapped_params: &KarlinParams,
) -> i32 {
    let concat_qlen = contexts
        .last()
        .map(|ctx| ctx.frame_base + ctx.aa_len as i32)
        .unwrap_or(0);
    let avg_qlen = concat_qlen / contexts.len().max(1) as i32;
    crate::algorithm::tblastx::ncbi_cutoffs::cutoff_score_sum_stats(
        avg_qlen,
        avg_subject_length,
        crate::core::link_hsps::BLAST_GAP_DECAY_RATE_GAPPED,
        gapped_params,
    )
}

/// psiblast's `-sum_stats true` in a search with the matrix (round 1): the lowered cutoff of the
/// HSPs the search keeps and the linking parameters of the preliminary stage.
pub(super) struct BlastpSumStats {
    /// The sum statistics' cutoff (`sum_stats_hit_cutoff`).
    hit_cutoff: i32,
    /// The preliminary stage's linking parameters (`cutoff_small_gap` = the words' smallest
    /// cutoff, `word_params->cutoff_score_min`).
    pub(super) preliminary: crate::core::link_hsps::EvenGapLinkParams,
}

impl BlastpSumStats {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:396-416
    pub(super) fn new(
        contexts: &[QueryContext],
        avg_subject_length: i32,
        gapped_params: &KarlinParams,
        word_cutoff_scores: &[i32],
    ) -> Self {
        let word_cutoff_score_min = contexts
            .iter()
            .zip(word_cutoff_scores)
            .filter(|(ctx, _)| ctx.is_valid)
            .map(|(_, &cutoff)| cutoff)
            .min()
            .unwrap_or(i32::MAX);
        Self {
            hit_cutoff: sum_stats_hit_cutoff(contexts, avg_subject_length, gapped_params),
            preliminary: crate::core::link_hsps::EvenGapLinkParams::preliminary(
                word_cutoff_score_min,
            ),
        }
    }

    /// Sets `cutoff_score` to the smaller of `new_cutoff` and `cutoff_score` for each valid context.
    pub(super) fn lower_hit_cutoffs(
        &self,
        contexts: &[QueryContext],
        mut cutoffs: Vec<i32>,
    ) -> Vec<i32> {
        for (cutoff, ctx) in cutoffs.iter_mut().zip(contexts) {
            if ctx.is_valid {
                *cutoff = (*cutoff).min(self.hit_cutoff);
            }
        }
        cutoffs
    }
}

/// The cutoffs of a BLASTP preliminary stage with the matrix at one subject length
/// (`blastp_prelim_cutoffs`), for each context.
struct BlastpPrelimCutoffs {
    /// The hit-saving cutoffs (`hit_params->cutoffs[context].cutoff_score`), lowered by psiblast's
    /// `-sum_stats true`.
    hit_cutoff_scores: Vec<i32>,
    /// The initial words' cutoffs (`word_params->cutoffs[context].cutoff_score`).
    word_cutoff_scores: Vec<i32>,
    /// psiblast's `-sum_stats true`: the sum statistics and the preliminary linking.
    sum_stats: Option<BlastpSumStats>,
}

/// The cutoffs of a BLASTP preliminary stage with the matrix at one subject length: the setup's
/// (`cutoff_subject_length`, `-comp_based_stats`, the preliminary stage's search spaces) or, with
/// the statistics of each subject, the subject's (`BLAST_OneSubjectUpdateParameters`: its length,
/// a composition-based statistics of 0, its search spaces). The hit-saving cutoffs
/// (`cutoff_score_max`) come first, the words' cutoffs read them before the sum statistics lower
/// them (`BlastInitialWordParametersUpdate` reads `cutoff_score_max`), and with psiblast's
/// `-sum_stats true` the sum statistics at that length take the words' smallest cutoff as the
/// linking's small-gap cutoff (`BlastLinkHSPParametersUpdate`) and lower each hit-saving cutoff to
/// their own where that is smaller (`cutoff_score`). `exact_spaces` are the contexts' effective
/// search spaces (0 where none is exact). The PSSM rounds' counterpart is
/// `pssm_hooks::psi_prelim_cutoffs`.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1015-1025 (cited at
/// `call_setup::subject_parameters`)
fn blastp_prelim_cutoffs(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    subject_length: i32,
    word_subject_length: i32,
    mode: BlastpCompositionMode,
    exact_spaces: &[Option<i64>],
    gapped_params: &KarlinParams,
    gapped_gumbel: Option<&BlastGumbelBlk>,
) -> BlastpPrelimCutoffs {
    if args.ungapped {
        blastp_prelim_cutoffs_of::<true>(
            args,
            contexts,
            subject_length,
            word_subject_length,
            mode,
            exact_spaces,
            gapped_params,
            gapped_gumbel,
        )
    } else {
        blastp_prelim_cutoffs_of::<false>(
            args,
            contexts,
            subject_length,
            word_subject_length,
            mode,
            exact_spaces,
            gapped_params,
            gapped_gumbel,
        )
    }
}

/// `blastp_prelim_cutoffs` of a search whose mode is `UNGAPPED` (`args.ungapped`), for the
/// callers that decide the mode once outside their loops (`call_setup::subject_parameters`): the
/// instance of a gapped search reads the call's gapped block and takes the gapped word cutoff for
/// every context, with no test of the mode in the loops over the contexts. An ungapped search
/// reads each context's block (`call_setup::block_of`) and takes the word cutoff of
/// `BlastInitialWordParametersUpdate` for its programs without gapped extension
/// (`blastp_ungapped_word_cutoff_score_for_context`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:1015-1025 (cited at
/// `call_setup::subject_parameters`)
#[inline(always)]
fn blastp_prelim_cutoffs_of<const UNGAPPED: bool>(
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    subject_length: i32,
    word_subject_length: i32,
    mode: BlastpCompositionMode,
    exact_spaces: &[Option<i64>],
    gapped_params: &KarlinParams,
    gapped_gumbel: Option<&BlastGumbelBlk>,
) -> BlastpPrelimCutoffs {
    let hit_cutoff_scores: Vec<i32> = contexts
        .iter()
        .zip(exact_spaces)
        .map(|(ctx, space)| {
            blastp_hit_cutoff_score_for_context(
                ctx,
                args.evalue,
                subject_length,
                mode,
                call_setup::block_of::<UNGAPPED>(gapped_params, ctx),
                gapped_gumbel,
                space.unwrap_or(0),
            )
        })
        .collect();
    let word_cutoff_scores: Vec<i32> = hit_cutoff_scores
        .iter()
        .zip(contexts)
        .map(|(&hit_cutoff_score, ctx)| {
            if UNGAPPED {
                blastp_ungapped_word_cutoff_score_for_context(
                    hit_cutoff_score,
                    ctx,
                    args.gap_trigger_bits,
                    word_subject_length,
                )
            } else {
                blastp_word_cutoff_score_for_context(hit_cutoff_score, ctx, args.gap_trigger_bits)
            }
        })
        .collect();
    let sum_stats = args
        .sum_stats
        .then(|| BlastpSumStats::new(contexts, subject_length, gapped_params, &word_cutoff_scores));
    let hit_cutoff_scores = match &sum_stats {
        Some(sum_stats) => sum_stats.lower_hit_cutoffs(contexts, hit_cutoff_scores),
        None => hit_cutoff_scores,
    };
    BlastpPrelimCutoffs {
        hit_cutoff_scores,
        word_cutoff_scores,
        sum_stats,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:969-1017
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1369-1439
//
// For multi-sequence subject sets NCBI keeps the hit-saving cutoff baseline
// fixed from setup-time sequence-source statistics; it does not recompute
// blastp protein cutoffs for each subject unless `db_length == 0`.
// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:133-135
// NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_multiseq.cpp:226-240
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_seqsrc.h:205
// The local `-subject` sequences are a multi-sequence source: the shortest record (an
// empty one counts as 0) is raised to 10.
#[inline]
fn blastp_preliminary_cutoff_subject_length(subjects: &[EncodedProtein]) -> i32 {
    const BLAST_SEQSRC_MINLENGTH: i32 = 10;
    subjects
        .iter()
        .map(|subject| subject.aa_len as i32)
        .min()
        .unwrap_or(i32::MAX)
        .max(BLAST_SEQSRC_MINLENGTH)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:234-242
// The preliminary HSP list of a subject for a query with the compact element, with the E-values
// `Blast_HSPListGetEvalues` gives (`x_prelim_compact.rs`).
#[allow(clippy::too_many_arguments)]
fn x_build_compact_preliminary_hsp_list<'h>(
    preliminary_hits: impl ExactSizeIterator<Item = &'h BlastpPreliminaryHsp>,
    oid: u32,
    query_index: u32,
    query_length: usize,
    subject_length: usize,
    karlin_params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    search_space: &crate::stats::search_space::SearchSpace,
    link: Option<&crate::core::link_hsps::EvenGapLinkParams>,
) -> BlastpHspList<XPrelimHsp> {
    if let Some(link) = link {
        return x_link_compact_preliminary_hsp_list(
            preliminary_hits,
            oid,
            query_index,
            query_length,
            subject_length,
            karlin_params,
            search_space,
            link,
        );
    }
    let search_space = search_space.effective_space;
    let hsps: Vec<XPrelimHsp> = match gumbel_params {
        Some(gumbel_params) => preliminary_hits
            .map(|hit| XPrelimHsp {
                hsp: *hit,
                e_value: blast_spouge_stoe(
                    hit.raw_score,
                    karlin_params,
                    gumbel_params,
                    query_length as i32,
                    subject_length as i32,
                ),
            })
            .collect(),
        None => preliminary_hits
            .map(|hit| XPrelimHsp {
                hsp: *hit,
                e_value: search_space
                    * (-karlin_params.lambda * hit.raw_score as f64 + karlin_params.k.glibc_ln())
                        .glibc_exp(),
            })
            .collect(),
    };
    let mut hsp_list = BlastpHspList {
        oid,
        query_index,
        hsps,
        best_evalue: f64::MAX,
    };
    super::hsp::update_best_evalue(&mut hsp_list);
    hsp_list
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:870-873
// With sum statistics the preliminary list of a subject takes the e-values of the linked sets
// (`core::link_hsps`) in the place of `Blast_HSPListGetEvalues`, and is sorted by score
// (`BLAST_LinkHsps` ends with `Blast_HSPListSortByScore`). NCBI links the subject's list of all
// queries; the linking never joins HSPs of different contexts and the sort is a total order, so
// each query's list is the same.
#[allow(clippy::too_many_arguments)]
fn x_link_compact_preliminary_hsp_list<'h>(
    preliminary_hits: impl ExactSizeIterator<Item = &'h BlastpPreliminaryHsp>,
    oid: u32,
    query_index: u32,
    query_length: usize,
    subject_length: usize,
    karlin_params: &KarlinParams,
    search_space: &crate::stats::search_space::SearchSpace,
    link: &crate::core::link_hsps::EvenGapLinkParams,
) -> BlastpHspList<XPrelimHsp> {
    let mut hsp_list = BlastpHspList {
        oid,
        query_index,
        hsps: preliminary_hits
            .map(|hit| XPrelimHsp {
                hsp: *hit,
                e_value: f64::MAX,
            })
            .collect(),
        best_evalue: f64::MAX,
    };
    super::hsp::link_hsp_list(
        &mut hsp_list,
        blastp_link_context(karlin_params, query_length, search_space),
        subject_length as i32,
        link,
    );
    hsp_list
}

/// The context of a query as the linking reads it: the gapped Karlin block (`kbp_gap`), the
/// query's length and its length adjustment and effective search space (`query_info->contexts`).
pub(super) fn blastp_link_context(
    karlin_params: &KarlinParams,
    query_length: usize,
    search_space: &crate::stats::search_space::SearchSpace,
) -> crate::core::link_hsps::LinkContext {
    crate::core::link_hsps::LinkContext {
        lambda: karlin_params.lambda,
        log_k: karlin_params.k.glibc_ln(),
        query_length: query_length as i32,
        length_adjustment: search_space.length_adjustment as i32,
        eff_searchsp: search_space.effective_space as i64,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:151-178
pub(super) fn blastp_hsp_from_preliminary(hit: &BlastpPreliminaryHsp, e_value: f64) -> BlastpHsp {
    BlastpHsp {
        num: 1,
        identity: 0.0,
        length: usize::try_from(hit.query_end - hit.query_start)
            .expect("preliminary blastp HSP query span must be non-negative"),
        mismatch: 0,
        gapopen: 0,
        q_start: usize::try_from(hit.query_start + 1)
            .expect("preliminary blastp query start must be non-negative"),
        q_end: usize::try_from(hit.query_end)
            .expect("preliminary blastp query end must be non-negative"),
        s_start: usize::try_from(hit.subject_start + 1)
            .expect("preliminary blastp subject start must be non-negative"),
        s_end: usize::try_from(hit.subject_end)
            .expect("preliminary blastp subject end must be non-negative"),
        gapped_q_start: hit.gapped_query_start,
        gapped_s_start: hit.gapped_subject_start,
        e_value,
        bit_score: 0.0,
        num_ident: 0,
        query_context: hit.query_context,
        query_frame: hit.query_frame,
        subject_frame: hit.subject_frame,
        query_length: hit.query_length,
        q_idx: hit.q_idx,
        s_idx: hit.s_idx,
        raw_score: hit.raw_score,
        gap_info: None,
        num_positives: 0,
        // Blast_HSPInit callocs the HSP: no composition adjustment yet.
        comp_adjustment_method: 0,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:516-518
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:776-779
#[inline]
fn blastp_diag_coord(offset_delta: i32, diag_mask: i32) -> usize {
    usize::try_from(offset_delta & diag_mask)
        .expect("NCBI BLAST masked diagonal index must be non-negative")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:358-372
fn blastp_save_init_hsp(
    init_hsps: &mut Vec<InitHSP>,
    q_start_absolute: i32,
    s_start: i32,
    q_off_absolute: i32,
    s_off: i32,
    len: i32,
    score: i32,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:336-341
    if init_hsps.len() == init_hsps.capacity() {
        x_ih_pool::before_full_push(init_hsps);
    }
    init_hsps.push(InitHSP {
        q_seed_absolute: q_off_absolute,
        s_seed: s_off,
        ungapped_data: BlastpUngappedData {
            q_start: q_start_absolute,
            s_start,
            length: len,
            score,
        },
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2370-2374
//
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_query_info.h:85-88
#[derive(Clone, Copy, Debug, Default)]
struct BlastpContextLookupBounds {
    min_length: i32,
    max_length: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:178-181
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:227-228
#[inline]
fn blastp_context_lookup_bounds(contexts: &[QueryContext]) -> BlastpContextLookupBounds {
    let mut min_length = i32::MAX;
    let mut max_length = 0i32;
    for ctx in contexts {
        let length =
            i32::try_from(ctx.aa_len).expect("NCBI BLAST protein query length must fit Int4");
        min_length = min_length.min(length);
        max_length = max_length.max(length);
    }
    if min_length == i32::MAX {
        min_length = 0;
    }
    BlastpContextLookupBounds {
        min_length,
        max_length,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-213
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:192-201
// The compressed lookup tables that NCBI builds for the queries of this search, as ranges of
// `locations` (the lookup segments of all queries, in query order): one for each query batch
// (`query_batch_size`: GetQueryBatchSize's 10000 letters, or the environment variable
// BATCH_SIZE) that NCBI does not split, none for a split batch, and one for each query chunk
// of a split batch, which is a search of its own here (`search_split_query_batches`,
// `query_chunk`) whose queries make that one table.
fn ncbi_compressed_lookup_tables(
    query_batches: &[seq_range::RangedBatch],
    query_lengths: &[usize],
    contexts: &[QueryContext],
    locations: &[(i32, i32)],
    query_chunk: bool,
    splits_queries: bool,
) -> Result<Vec<std::ops::Range<usize>>> {
    if query_chunk {
        return Ok(vec![0..locations.len()]);
    }
    let sizes = protein_query_split_sizes(10000)?;
    // The first lookup segment of query `q_idx` or of a query after it.
    let first_segment = |q_idx: usize| {
        let context = contexts.partition_point(|context| (context.q_idx as usize) < q_idx);
        contexts.get(context).map_or(locations.len(), |context| {
            locations.partition_point(|&(left, _)| left < context.frame_base)
        })
    };
    Ok(query_batches
        .iter()
        .map(|batch| batch.searched.clone())
        .filter(|searched| !searched.is_empty())
        .filter(|searched| {
            // An ungapped search splits no query batch (`search_split_query_batches`).
            !splits_queries
                || crate::algorithm::blastp::protein_query_split::split_protein_batch(
                    &query_lengths[searched.clone()],
                    sizes,
                )
                .is_none()
        })
        .map(|searched| first_segment(searched.start)..first_segment(searched.end))
        .collect())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-242
#[cfg(test)]
#[inline]
fn blastp_context_idx_for_query_offset(contexts: &[QueryContext], concat_off: i32) -> usize {
    blastp_context_idx_for_query_offset_with_bounds(
        contexts,
        blastp_context_lookup_bounds(contexts),
        concat_off,
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:225-236
#[inline]
fn blastp_context_idx_for_query_offset_with_bounds(
    contexts: &[QueryContext],
    bounds: BlastpContextLookupBounds,
    concat_off: i32,
) -> usize {
    let size = contexts.len();
    if size == 0 {
        return 0;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_query_info.c:219-242
    if size == 1 {
        return 0;
    }

    let mut lo = 0usize;
    let mut hi = size;
    if bounds.min_length > 0 && bounds.max_length > 0 {
        let n = usize::try_from(concat_off)
            .expect("NCBI BLAST query offset for context lookup must be non-negative");
        let max_length =
            usize::try_from(bounds.max_length).expect("NCBI BLAST max query length must fit usize");
        let min_length =
            usize::try_from(bounds.min_length).expect("NCBI BLAST min query length must fit usize");
        lo = (n / (max_length + 1)).min(size - 1);
        hi = (n / (min_length + 1) + 1).min(size);
    }

    while lo < hi.saturating_sub(1) {
        let mid = (lo + hi) / 2;
        if concat_off < contexts[mid].frame_base {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    lo
}

#[inline]
fn blastp_get_ungapped_hsp_context_with_bounds(
    contexts: &[QueryContext],
    bounds: BlastpContextLookupBounds,
    init_hsp: &InitHSP,
) -> usize {
    blastp_context_idx_for_query_offset_with_bounds(contexts, bounds, init_hsp.q_seed_absolute)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3870-3879
#[inline]
fn blastp_advance_preliminary_hit_index(
    hit_index: usize,
    redo_index: Option<usize>,
    next_same_query_index: &[Option<usize>],
) -> usize {
    if let Some(redo_index) = redo_index {
        if hit_index < redo_index && !next_same_query_index.is_empty() {
            return next_same_query_index
                .get(hit_index)
                .and_then(|index| *index)
                .map(|index| index.min(redo_index))
                .unwrap_or(redo_index);
        }
    }
    hit_index + 1
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2374-2389
#[inline]
fn blastp_adjust_initial_hsp_offsets(init_hsp: &mut InitHSP, query_start: i32) {
    init_hsp.q_seed_absolute -= query_start;
    init_hsp.ungapped_data.q_start -= query_start;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2399-2466
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2445-2467
// Inlined into both callers (`prelim_kernel`'s `gapped_subject` and `ungapped_subject`): run once
// per initial HSP, a call of its own cost gapped blastp 0.019 % more instructions than dev's one
// inlined caller (callgrind, `O4b` A/B).
#[inline(always)]
fn blastp_adjust_hsp_offsets_and_get_query_data<'a>(
    concatenated_query: &'a [u8],
    contexts: &'a [QueryContext],
    context_lookup_bounds: BlastpContextLookupBounds,
    init_hsp: &mut InitHSP,
) -> (usize, &'a QueryContext, &'a [u8]) {
    let context_idx =
        blastp_get_ungapped_hsp_context_with_bounds(contexts, context_lookup_bounds, init_hsp);
    let ctx = &contexts[context_idx];
    let query_start = ctx.frame_base;

    assert!(
        init_hsp.ungapped_data.q_start >= query_start,
        "NCBI BLAST requires local-query adjustment to keep q_start non-negative"
    );
    blastp_adjust_initial_hsp_offsets(init_hsp, query_start);
    assert!(
        init_hsp.ungapped_data.q_start >= 0,
        "NCBI BLAST requires adjusted ungapped q_start >= 0"
    );
    let query_start_usize = usize::try_from(query_start)
        .expect("NCBI BLAST concatenated blastp query offsets must be non-negative");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:2426-2466
    //
    // `query_tmp.sequence` points at the current query context only; it does
    // not run into later concatenated queries. Keep the trailing NULLB so the
    // local slice matches the C buffer layout.
    let query_end_usize = query_start_usize
        .checked_add(ctx.aa_len + 1)
        .expect("NCBI BLAST local blastp query span must fit in usize");
    let query_sequence = &concatenated_query[query_start_usize..query_end_usize];
    (context_idx, ctx, query_sequence)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
const BLASTP_OFFSET_ARRAY_SIZE: i32 = 4096;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
#[inline]
fn blastp_offset_array_size(longest_chain: i32) -> i32 {
    BLASTP_OFFSET_ARRAY_SIZE + longest_chain.max(0)
}

enum BlastpRuntimeLookup {
    Aa(BlastAaLookupTable),
    // EXPERIMENT (LOSAT_X_BPLK16): the exact table of a block whose lookup query is at most 65,535
    // letters, with 16-bit hits and 8-byte cells (lookup/backbone.rs `BackboneCell16`); the same
    // hits in the same order. Only for a search with a window (the fused two-hit loop reads it).
    Aa16(BlastAaLookupTable<BackboneCell16>),
    Compressed(BlastCompressedAaLookupTable),
    // EXPERIMENT (LOSAT_X_BPCLUT): the compressed table built flat (lookup/x_clut_direct.rs); the
    // same cells, PV array and parameters.
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:104-107
    CompressedFlat(XFlatCompressedLookup),
}

impl BlastpRuntimeLookup {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:481-496
    #[inline]
    fn word_length(&self) -> i32 {
        match self {
            Self::Aa(lookup) => lookup.word_length,
            Self::Aa16(lookup) => lookup.word_length,
            Self::Compressed(lookup) => lookup.word_length,
            Self::CompressedFlat(lookup) => lookup.word_length,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/lookup_wrap.c:255-288
    #[inline]
    fn longest_chain(&self) -> i32 {
        match self {
            Self::Aa(lookup) => lookup.longest_chain,
            Self::Aa16(lookup) => lookup.longest_chain,
            Self::Compressed(lookup) => lookup.longest_chain,
            Self::CompressedFlat(lookup) => lookup.longest_chain,
        }
    }
}

#[derive(Clone, Copy)]
enum BlastpOffsetPairScanMode {
    OneHit,
    // The two-hit search scans with the fused loop; the offset-pair two-hit scan stays as the unit tests' reference.
    #[cfg(test)]
    TwoHit,
}

impl BlastpOffsetPairScanMode {
    /// aa_ungapped.c clamps the scan end for the two-hit search only.
    fn clamps_scan_end(self) -> bool {
        match self {
            Self::OneHit => false,
            #[cfg(test)]
            Self::TwoHit => true,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:86-93
#[inline(always)]
fn blastp_visit_scanned_offset_pairs<F>(offset_pairs: &[OffsetPair], hits: i32, visit: &mut F)
where
    F: FnMut(OffsetPair),
{
    let offset_pairs_ptr = offset_pairs.as_ptr();
    for i in 0..hits as usize {
        // SAFETY: the selected NCBI-shaped scan callback returns at most
        // `array_size`, which is set to `offset_pairs.len()` for this buffer.
        let pair = unsafe { *offset_pairs_ptr.add(i) };
        visit(pair);
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:762-769
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_util.c:112-127
fn blastp_for_each_offset_pair<F>(
    lookup: &BlastpRuntimeLookup,
    subject_sequence: &[u8],
    subject_length: usize,
    offset_pairs: &mut [OffsetPair],
    scan_mode: BlastpOffsetPairScanMode,
    mut visit: F,
) where
    F: FnMut(OffsetPair),
{
    let word_size = lookup.word_length();
    let seq_ranges = [(0, subject_length as i32)];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_size];

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:498-500
    if scan_mode.clamps_scan_end() && scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:492-505,763-769
    let array_size =
        i32::try_from(offset_pairs.len()).expect("NCBI BLAST offset pair array size must fit Int4");

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aascan.c:552-568
    match lookup {
        BlastpRuntimeLookup::Aa(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = s_blast_aa_scan_subject_one_range(
                    lookup,
                    subject_sequence,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
        // EXPERIMENT (LOSAT_X_BPLK16): the table of 16-bit hits is built only for a search with a
        // window, which scans with the fused two-hit loop; this offset-pair scan (a search without a
        // window) never gets one.
        BlastpRuntimeLookup::Aa16(_) => {
            unreachable!("a table of 16-bit hits is built only for a search with a window")
        }
        BlastpRuntimeLookup::Compressed(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = lookup.scan_subject(
                    subject_sequence,
                    &seq_ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
        // EXPERIMENT (LOSAT_X_BPCLUT): the same scan on the flat table (`scan_subject`'s statements, a
        // cell's hits copied from its slice; lookup/x_clut_direct.rs).
        BlastpRuntimeLookup::CompressedFlat(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = lookup.scan_subject::<X_BLASTP>(
                    subject_sequence,
                    &seq_ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
    }
}
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:763-769 (see
// `blastp_for_each_offset_pair`)
/// [`blastp_for_each_offset_pair`] for the tables of another matrix than BLOSUM62 (item SC2b), with
/// the same statements: the scan of the flat compressed table is a copy of its own
/// (`X_BLASTP_MATRIX`), so that BLOSUM62's keeps its one caller.
fn blastp_for_each_offset_pair_of<F>(
    lookup: &BlastpRuntimeLookup,
    subject_sequence: &[u8],
    subject_length: usize,
    offset_pairs: &mut [OffsetPair],
    scan_mode: BlastpOffsetPairScanMode,
    mut visit: F,
) where
    F: FnMut(OffsetPair),
{
    let word_size = lookup.word_length();
    let seq_ranges = [(0, subject_length as i32)];
    let mut scan_range = [0, seq_ranges[0].0, seq_ranges[0].1 - word_size];
    if scan_mode.clamps_scan_end() && scan_range[2] < scan_range[1] {
        scan_range[2] = scan_range[1];
    }
    let array_size =
        i32::try_from(offset_pairs.len()).expect("NCBI BLAST offset pair array size must fit Int4");
    match lookup {
        BlastpRuntimeLookup::Aa(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = s_blast_aa_scan_subject_one_range(
                    lookup,
                    subject_sequence,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
        BlastpRuntimeLookup::Aa16(_) => {
            unreachable!("a table of 16-bit hits is built only for a search with a window")
        }
        BlastpRuntimeLookup::Compressed(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = lookup.scan_subject(
                    subject_sequence,
                    &seq_ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
        BlastpRuntimeLookup::CompressedFlat(lookup) => {
            while scan_range[1] <= scan_range[2] {
                let hits = lookup.scan_subject::<X_BLASTP_MATRIX>(
                    subject_sequence,
                    &seq_ranges,
                    offset_pairs,
                    array_size,
                    &mut scan_range,
                );
                blastp_visit_scanned_offset_pairs(offset_pairs, hits, &mut visit);
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1163-1169
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3036-3062
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:4209-4319
fn validate_requested_blastp_support(args: &ResolvedBlastpArgs) -> Result<()> {
    let program = args.app.name();
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:917-918
    // LOSAT's limits, checked where NCBI starts the search: the options that NCBI runs and
    // LOSAT's blastp engine does not implement. `-ungapped` is searched (series NP, O4b: the
    // ungapped search, `prelim_kernel::BlastpPrelimKernel::ungapped_subject`).
    // An infinite -evalue, or DBL_MAX with the composition-based redo (decision D12;
    // `protein_options_check::check_evalue_limit`). A smaller -evalue such as 1e308 runs.
    crate::blastinput::protein_options_check::check_evalue_limit(
        args.evalue,
        args.comp_based_stats.is_enabled(),
        program,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:1486-1487
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:678-686
    // `-use_sw_tback` (item X4): BLASTP runs NCBI's Smith-Waterman redo with composition-based
    // statistics 1 to 3 (`kappa::smith_waterman`). Without them (`-comp_based_stats 0`, and
    // IDENTITY's reset, `identity_cbs_reset`) NCBI runs the redo too,
    // without an adjustment, on the preliminary lists in subject order (the stream is not sorted
    // by score), which LOSAT's redo loop does not read yet (the next item, X4-c0). psiblast's redo
    // of a PSSM keeps its refusal (the position-based Smith-Waterman redo).
    if args.use_sw_tback {
        if args.app != BlastpApp::Blastp || args.psi_search.is_some() {
            bail!("{}", crate::cli::run_refusal(program, "use_sw_tback"));
        }
        // `-ungapped` (which only `-comp_based_stats 0` allows): NCBI crashes at the redo's set-up
        // in the traceback of each query batch it searches, before any match is redone, and LOSAT
        // refuses there (`redo_run::sw_ungapped_check`, in `run_call`), after the reports of the
        // batches before, which have no valid context. A run without such a batch ends as NCBI's
        // does.
        if !args.ungapped && !args.comp_based_stats.is_enabled() {
            bail!(
                "-use_sw_tback without composition-based statistics (-comp_based_stats 0, and -matrix IDENTITY, which NCBI searches without them) is not supported by LOSAT's {program} yet"
            );
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_options.c:1389-1392
    // The compressed lookup table runs at each of its word sizes, 5 to 7
    // (`validate_protein_options` has refused the others), and so does the lookup table at 2
    // (blastp-short, item SC3), 3 and 4 (item SC8b); this refuses no word size that the options
    // handler lets through.
    if !matches!(args.word_size, 2..=4)
        && args.lookup_table_type != BlastpLookupTableType::CompressedAaLookupTable
        && args.psi_search.is_none()
    {
        bail!(
            "-word_size {} (the lookup table) is not supported by LOSAT's {program}",
            args.word_size
        );
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:258-271
    // BLASTP reads the gap costs of the search at run time (the gapped alignments, the chaining
    // and the composition-based redo) and their Karlin and Gumbel values from NCBI's table
    // (`lookup_protein_params`, `lookup_protein_gumbel_params`), and the options handler has
    // checked them against the table (`validate_protein_options`). It reads the matrix at run
    // time too: the lookup table's words and the contexts' Karlin blocks (`blastp_aa_lookup`,
    // `blastp_compressed_lookup`), the ideal Karlin block (`ideal_karlin_params_of`), the
    // ungapped extensions (`x_candidate_of!`), the gapped ones (`BlastpScoreMatrix`), the
    // chaining (the typed spelling) and the composition-based redo (its matrix information and
    // workspace). Every matrix of NCBI's tables runs with any of its rows; IDENTITY, which has no
    // composition data, in mode 0 (the options handler resets the others, `identity_cbs_reset`;
    // item SC4b). psiblast's first round is BLASTP's search, and its later rounds search the PSSM
    // of the previous one with the matrix and the gap costs of the search (its background and
    // frequency ratios, `score_block`; the contexts' Karlin blocks, `pssm_hooks::psi_contexts`;
    // the positives, `pssm_search::traceback0::final_hsp`; and the table's gapped block): every
    // matrix of NCBI's tables but IDENTITY, which psiblast's options handler refuses as NCBI's
    // does (`std_matrix_only`). A PSSM query (`-in_pssm`, `psi_search` before the search) keeps
    // BLOSUM62 11/1 (line PS, PS7b).
    if args.psi_search.is_some()
        && (args.scoring.matrix != ScoringMatrix::Blosum62
            || args.scoring.gap_open != 11
            || args.scoring.gap_extend != 1)
    {
        bail!(
            "the matrix {} with gap costs {}/{} is not supported by LOSAT's {program}",
            args.matrix_name,
            args.scoring.gap_open,
            args.scoring.gap_extend
        );
    }
    // Unified P-values (a `u` after the mode, BLASTP only) run in the composition-based redo
    // (`kappa::evalue_adjust_contexts`: each query's E-values take the lengths and the search
    // space that NCBI 2.17.0 reads, the next query's of the batch; item X5 of line RG). A search
    // that IDENTITY resets to mode 0 has none (`identity_cbs_reset`).
    // Every composition mode runs (item SC5 of line NP). Mode 0 runs BLASTP's ordinary traceback
    // (`traceback.rs`), psiblast's first round of a FASTA query too (with psiblast's hit list
    // size, `traceback_hitlist_size`), and the PSSM's of psiblast `-in_pssm` and of its later
    // rounds (`pssm_search::traceback0`). Modes 1 to 3 run the composition-based redo, whose
    // matrix adjustment rule the mode picks (`blast_adjust_scores`: 1 scales the matrix by the
    // ratio of the Lambdas, 3 always adjusts it to the target relative entropy), psiblast's first
    // round of a FASTA query too; psiblast's PSSM searches run modes 2 and 3 as 1
    // (`PsiBlastSetupScoreBlock`, the options), and mode 1 as the position-based redo
    // (`pssm_search::redo1`). rpsblast's search reads its modes in `rps_search` (item PS11b of
    // line PS).
    Ok(())
}

/// `protein_options_check::check_identity_compressed_lookup` of a BLASTP or psiblast search.
fn check_identity_compressed_lookup_of(args: &ResolvedBlastpArgs) -> Result<()> {
    crate::blastinput::protein_options_check::check_identity_compressed_lookup(
        &args.matrix_name,
        args.lookup_table_type == BlastpLookupTableType::CompressedAaLookupTable,
        args.word_size as i32,
        args.app.ncbi_name(),
        args.app.name(),
    )
}

/// The report of one requested BLASTP format: the pairwise report (outfmt 0, and the
/// query-anchored outfmt 1 to 4: `PairwiseView`), the default outfmt 6 and 7, the tabular reports with a field list or a delimiter and outfmt 10
/// and 20 (`tabular`), the Seq-annot of Seq-aligns of each query (outfmt 8, 9, 12;
/// `seqalign`), BLAST XML (outfmt 5; `blast_xml`), or the XML2 / JSON documents (outfmt 15 and
/// 16; `blast_xml2`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BlastpReport {
    /// Outfmt 0 to 4 (`Text` holds 6 and 7 only).
    Pairwise(crate::report::pairwise::PairwiseView),
    Text(OutputFormat),
    /// NCBI's `CBlastTabularInfo` with a field list or a delimiter
    /// (`tabular::custom_tabular_spec`, `TabularSpec::of_app`), through the shared writer
    /// (`report::tabular`, `tabular::write_blastp_tabular`).
    CustomTabular(&'static crate::report::tabular::TabularSpec),
    Seqalign(crate::report::serial::SerialFormat),
    Xml,
    Xml2(crate::report::blast_xml2::Xml2Kind),
    /// The BLAST archive (outfmt 11, `archive`): one document per NCBI query batch.
    Archive,
}

impl BlastpReport {
    /// Whether the report's writer writes the warnings that come after the last query's report
    /// before its epilog (`QueryWarnings::after_last_query`): psiblast's reports, outfmt 0, 6, 7,
    /// 10 and 20.
    fn writes_after_last_query(self) -> bool {
        matches!(
            self,
            BlastpReport::Pairwise(_)
                | BlastpReport::Text(OutputFormat::Tabular | OutputFormat::TabularWithComments)
                | BlastpReport::CustomTabular(_)
        )
    }

    /// Whether NCBI prefetches the subjects that each query batch reports before the batch's
    /// reports and the report shows no taxonomy of them (outfmt 0, 5, 6 and 7; with
    /// `PRE_FETCH_SEQS_LIMIT=2147483647` also 8, 9 and 12, `database::prefetches_subjects`): as
    /// the call's first report, the taxonomy database's messages of that prefetch go with the
    /// query warnings.
    fn prefetches_subjects(self) -> bool {
        crate::algorithm::common::database::prefetches_subjects(match self {
            BlastpReport::Pairwise(view) => Some(view.outfmt()),
            BlastpReport::Text(outfmt) => Some(outfmt as i32),
            // NCBI's `s_IsPrefetchFormat` prefetches for every tabular format (6, 7, 10 and 20):
            // they read the subjects as 6 and 7 do; the rows' own lookups of the taxonomy names
            // come with the rows (`report::tabular::TabularNames`).
            BlastpReport::CustomTabular(spec) => Some(match spec.kind {
                crate::report::tabular::TabularKind::Commented => 7,
                _ => 6,
            }),
            BlastpReport::Xml => Some(5),
            BlastpReport::Seqalign(serial) => Some(serial.outfmt()),
            // The archive fetches nothing; XML2 and JSON write the messages at each `Hit`
            // (`TaxdbPrinter`).
            BlastpReport::Archive | BlastpReport::Xml2(_) => None,
        })
    }

    /// Whether the report reads the sequence that NCBI's formatters find by a subject's Seq-align
    /// id: every one but the Seq-annots, the archive and a table of fields that read nothing of
    /// the subject (`report::tabular::reads_subject`).
    fn reads_subjects(self) -> bool {
        match self {
            BlastpReport::Seqalign(_) | BlastpReport::Archive => false,
            BlastpReport::CustomTabular(spec) => {
                spec.reads_subject(crate::blastinput::app::Program::Blastp)
            }
            BlastpReport::Pairwise(_)
            | BlastpReport::Text(_)
            | BlastpReport::Xml
            | BlastpReport::Xml2(_) => true,
        }
    }

    /// Whether the report is tabular (outfmt 6, 7, 10 and 20), which reads the tabular subject
    /// IDs (`report::tabular::writes_rows`: not the organism report).
    fn is_tabular(self) -> bool {
        match self {
            BlastpReport::Text(OutputFormat::Tabular | OutputFormat::TabularWithComments) => true,
            BlastpReport::CustomTabular(spec) => spec.writes_rows(),
            _ => false,
        }
    }

    /// The field list of a custom tabular report (`None`: the default fields, or no table).
    fn tabular_spec(self) -> Option<&'static crate::report::tabular::TabularSpec> {
        match self {
            BlastpReport::CustomTabular(spec) => Some(spec),
            _ => None,
        }
    }

    /// The subject ids that the report prints (`PrintedIds`): none but in a table.
    fn printed_ids(self) -> crate::report::tabular::PrintedIds {
        match self {
            BlastpReport::CustomTabular(spec) => crate::report::tabular::PrintedIds::of(Some(spec)),
            _ if self.is_tabular() => crate::report::tabular::PrintedIds::of(None),
            _ => crate::report::tabular::PrintedIds::default(),
        }
    }
}

/// The Seq-annot format of a `-outfmt` value (outfmt 8, 12), read as NCBI reads the value
/// (`app::parse_formatting_string`, which drops a custom specification of these formats);
/// `None` for the other formats.
fn blastp_seqalign_format(outfmt: &str) -> Option<crate::report::serial::SerialFormat> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    crate::report::serial::SerialFormat::of(crate::blastinput::app::ReportFormat::from_number(
        choice.number,
    )?)
}

/// Whether NCBI's result of each query has no alignment set: a query of a batch that was not
/// searched (`query_batch_skipped`), but the merge of several volumes gives each query a set, empty
/// for such a batch, so outfmt 7 counts `# 0 hits found` there (review A of PS12a, S2's rpsblast
/// half), XML2 and JSON write `hits` and the plain message, the organism report its header
/// (PS15a-2 review B1).
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/rpsblast_local.cpp:155-180
fn no_alignment_sets(
    query_batch_skipped: &[bool],
    rps_multi_volume: bool,
) -> std::borrow::Cow<'_, [bool]> {
    if rps_multi_volume {
        std::borrow::Cow::Owned(vec![false; query_batch_skipped.len()])
    } else {
        std::borrow::Cow::Borrowed(query_batch_skipped)
    }
}

/// Whether rpsblast recounts the identities of its HSPs for the reports of the call: none of them is
/// a custom table, whose spec reads them its own way, or the BLAST archive or a Seq-annot (outfmt 8,
/// 9, 12), whose Seq-aligns show the traceback's count (`PrelimRun::rps_recount_identities`;
/// `ReportFormat::writes_seqaligns`, item PS15a-2).
fn recounts_rps_identities(reports: &[BlastpReport]) -> bool {
    reports.iter().all(|report| {
        !matches!(
            report,
            BlastpReport::CustomTabular(_) | BlastpReport::Archive | BlastpReport::Seqalign(_)
        )
    })
}

/// The pairwise view (`PairwiseView::of`) of the format number `number` where `app` writes it:
/// outfmt 0, and the query-anchored 1 to 4 of BLASTP and RPSBLAST (`Program::writes`; psiblast
/// refuses those, item PS15b). NCBI's application writes 0 to 4 through one frame
/// (`PrintProlog`, the query batches, `PrintOneResultSet`), and so does LOSAT; `None` for the
/// other formats.
fn blastp_pairwise_view_of(
    app: BlastpApp,
    number: i32,
) -> Option<crate::report::pairwise::PairwiseView> {
    crate::blastinput::app::ReportFormat::from_number(number)
        .and_then(crate::report::pairwise::PairwiseView::of)
        .filter(|view| {
            *view == crate::report::pairwise::PairwiseView::Pairwise || app != BlastpApp::Psiblast
        })
}

/// The format of a `-outfmt` value that is outfmt 0, 6 or 7 to LOSAT's own writers (not a custom
/// table, `tabular::custom_tabular_spec`), its number read as NCBI reads it
/// (`app::parse_formatting_string`: `06` and `+6` are 6), on the command line (its value
/// `FormatChoice::normalized`) and in the library (the caller's value as given) alike; an error
/// for the other numbers.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2829
fn blastp_text_format(outfmt: &str) -> std::result::Result<OutputFormat, String> {
    match crate::blastinput::app::parse_formatting_string(outfmt).map(|choice| choice.number) {
        Ok(0) => Ok(OutputFormat::Pairwise),
        Ok(6) => Ok(OutputFormat::Tabular),
        Ok(7) => Ok(OutputFormat::TabularWithComments),
        _ => Err(format!(
            "Unsupported output format: {}. Supported: 0, 6, 7",
            OutputFormat::split_fields(outfmt).0
        )),
    }
}

/// The query-anchored view of a `-outfmt` value of 1 to 4 that `app` writes
/// (`blastp_pairwise_view_of`), read as NCBI reads the value; `None` for the other formats. The
/// words after the number (the library's value as given) are dropped, as NCBI drops them for
/// every format but the tables (`parse_formatting_string`) and the command line's value does
/// (`FormatChoice::normalized`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2845-2851
fn blastp_query_anchored_view(
    app: BlastpApp,
    outfmt: &str,
) -> Option<crate::report::pairwise::PairwiseView> {
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    blastp_pairwise_view_of(app, choice.number)
        .filter(|view| *view != crate::report::pairwise::PairwiseView::Pairwise)
}

/// Whether a `-outfmt` value is the BLAST archive (outfmt 11), read as NCBI reads it (a
/// specification after the number is ignored).
fn blastp_is_archive(outfmt: &str) -> bool {
    crate::blastinput::app::parse_formatting_string(outfmt).is_ok_and(|choice| {
        crate::blastinput::app::ReportFormat::from_number(choice.number)
            == Some(crate::blastinput::app::ReportFormat::ArchiveFormat)
    })
}

/// Whether a `-outfmt` value is BLAST XML (outfmt 5), read as NCBI reads it.
fn blastp_is_xml(outfmt: &str) -> bool {
    crate::blastinput::app::parse_formatting_string(outfmt).is_ok_and(|choice| {
        crate::blastinput::app::ReportFormat::from_number(choice.number)
            == Some(crate::blastinput::app::ReportFormat::Xml)
    })
}

/// The XML2 / JSON document of a `-outfmt` value (outfmt 13 to 16), read as NCBI reads it;
/// `None` for the other formats. Outfmt 13 and 14 write files next to `-out`: only the command
/// line reaches them (`resolve_for_formats` stops the other callers at NCBI's XInclude check).
fn blastp_xml2_kind(outfmt: &str) -> Option<crate::report::blast_xml2::Xml2Kind> {
    use crate::blastinput::app::ReportFormat;
    let choice = crate::blastinput::app::parse_formatting_string(outfmt).ok()?;
    crate::report::blast_xml2::Xml2Kind::of(ReportFormat::from_number(choice.number)?)
}

/// Whether the block search takes this XML2 / JSON document: all four (outfmt 13 to 16). Outfmt 13
/// and 14 (`Files13`, `Files14`) write one side file per report and a master file at `-out`; the
/// side files are made by the same `Xml2Document` as the single-stream forms, whose only state
/// that outlives a report (`written`, which numbers the side file and the master's list, `started`
/// and the held report of a `-subject` run) is in `Xml2Carry`, and `FsSideFiles` keeps nothing
/// between calls (`x_query_blocks.rs`, module doc). No NCBI counterpart: chooses LOSAT's path.
fn x_xml2_kind_takes_blocks(kind: crate::report::blast_xml2::Xml2Kind) -> bool {
    use crate::report::blast_xml2::Xml2Kind;
    matches!(
        kind,
        Xml2Kind::Single15 | Xml2Kind::Single16 | Xml2Kind::Files13 | Xml2Kind::Files14
    )
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1264-1283
fn write_blastp_outfmt7_header<W: Write>(
    writer: &mut W,
    context: &ReportContext,
    query_title: &[u8],
    fields: &[crate::report::tabular::TabularField],
    num_hits: usize,
    batch_skipped: bool,
    iteration: Option<u32>,
) -> std::io::Result<()> {
    let version_str = crate::wording::REPORT_VERSION;
    writeln!(
        writer,
        "# {} {}",
        context.program.to_uppercase(),
        version_str
    )?;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1297-1298
    // psiblast's tabular reports carry their round (`x_PrintTabularReport` of `results` and `itr_num`,
    // blast_format.cpp:1458): 1, or the round of the round driver (`psi_round`).
    if let Some(iteration) = iteration {
        writeln!(writer, "# Iteration: {iteration}")?;
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1305-1308
    // The title's bytes (`write_outfmt7_query_line`), as NCBI's reader keeps them
    // (`FastaRecord::title`).
    crate::report::outfmt6::write_outfmt7_query_line(writer, query_title)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:790-803
    if let Some(db) = context
        .database_name
        .as_ref()
        .or(context.subject_name.as_ref())
    {
        writeln!(writer, "# Database: {}", db)?;
    }
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1277-1278
    // A query of an unsearched batch has no Seq-align set, so no hit count is written.
    if batch_skipped {
        return Ok(());
    }
    if num_hits > 0 {
        let field_list = fields
            .iter()
            .map(|field| field.header_name())
            .collect::<Vec<_>>()
            .join(", ");
        writeln!(writer, "# Fields: {}", field_list)?;
    }
    writeln!(writer, "# {} hits found", num_hits)?;
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/showdefline.cpp:225-231
/// The subject ID of the tabular formats of every subject, by subject index
/// (`report::defline::tabular_subject_id`: the title's first word, or the local ID, and
/// for `lcl|Subject_...` the first word of `GenerateDefline`, which is `unnamed` for a
/// protein without a title), as the other programs'; `None` where `GenerateDefline` throws
/// (`common::tabular::SubjectId`), whose rows NCBI writes with its catch's values
/// (`report::tabular::SubjectLookup::Failed`) or stops at (`report::tabular::check_failed_lookups`).
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:863-870
fn tabular_subject_ids(
    subjects: &[FastaRecord],
) -> Vec<crate::algorithm::common::tabular::SubjectId> {
    subjects
        .iter()
        .map(|subject| {
            subject
                .subject_tabular_id(true)
                .ok()
                .flatten()
                .map(Arc::from)
        })
        .collect()
}

/// The tabular subject ids of a subject set under `-parse_deflines`: `tabular_subject_ids`, with
/// `sseqid`, `sacc` and `saccver` of the records in `parsed` (`ParsedSubject::tabular_ids`).
#[allow(clippy::type_complexity)]
fn record_tabular_ids(
    subjects: &[FastaRecord],
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    program: &str,
) -> Result<(
    Vec<crate::algorithm::common::tabular::SubjectId>,
    Vec<Option<[Arc<[u8]>; 2]>>,
)> {
    let mut ids = tabular_subject_ids(subjects);
    let mut accessions = vec![None; ids.len()];
    for (index, subject) in parsed.iter().enumerate() {
        let Some(subject) = subject else {
            continue;
        };
        // No ids where the lookup throws (`ParsedSubject::tabular_ids`).
        let found = subject.tabular_ids(true).map_err(|problem| {
            anyhow::anyhow!(
                "the subject record {} has {problem}, which is not supported by LOSAT's {program}",
                String::from_utf8_lossy(&subject.alignment_fasta())
            )
        })?;
        ids[index] = found.as_ref().map(|[id, _, _]| Arc::from(id.as_slice()));
        accessions[index] =
            found.map(|[_, accession, version]| [Arc::from(accession), Arc::from(version)]);
    }
    Ok((ids, accessions))
}

/// The tabular subject id field of one subject that the default fields print: `saccver`. It is
/// the same id as `sseqid` and `sacc` for a `-subject` record and for a database sequence without
/// parsed Seq-ids; a database sequence with parsed Seq-ids has its own three
/// (`ParsedSubject::tabular_ids`).
#[derive(Clone, Copy)]
struct SubjectIdFields<'a> {
    accession_version: &'a [u8],
}

impl<'a> SubjectIdFields<'a> {
    /// The fields of one id.
    fn same(id: &'a [u8]) -> Self {
        SubjectIdFields {
            accession_version: id,
        }
    }

    /// The fields of subject `s_idx`: its id among `ids` (empty for none; a failed lookup's
    /// rows are refused before, `check_failed_lookups`), and its `saccver` among `accessions`
    /// where it has one.
    fn of(
        ids: &'a [crate::algorithm::common::tabular::SubjectId],
        accessions: &'a [Option<[Arc<[u8]>; 2]>],
        s_idx: usize,
    ) -> Self {
        let id = crate::algorithm::common::tabular::printed_id(ids, s_idx);
        match accessions.get(s_idx).and_then(Option::as_ref) {
            Some([_, version]) => SubjectIdFields {
                accession_version: version,
            },
            None => SubjectIdFields::same(id),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:271-327
/// The Seq-aligns NCBI makes of the results `hit_lists`: one for each HSP of a `gapped` search,
/// one for each subject of an ungapped one (`-ungapped`, `shared::subject_seq_aligns`), which
/// the prefetch counts against `PRE_FETCH_SEQS_LIMIT` (`shared::prefetch_off`). Outfmt 7's
/// count takes those of the formatter, one for each HSP in both modes (`gapped` true).
fn blastp_seq_aligns(hit_lists: &[Option<BlastpHitList>], gapped: bool) -> usize {
    hit_lists
        .iter()
        .filter_map(|hit_list| hit_list.as_ref())
        .map(|hit_list| {
            hit_list
                .hsplist_array
                .iter()
                .take(hit_list.hsplist_count)
                .map(|hsp_list| {
                    crate::algorithm::common::database::subject_seq_aligns(
                        hsp_list.hsps.len(),
                        gapped,
                    )
                })
                .sum::<usize>()
        })
        .sum()
}

/// The rows of subjects whose lookups fail where NCBI has no result
/// (`report::tabular::check_failed_lookups`), before any format is written: each tabular
/// format's refusal is that format's (`FormatRefusal::check`). `parsed`: the subjects' parsed
/// Seq-ids (a database's listed titles'); `loads`: the reports' loads (`HitLoads::of`), whose
/// throwing subjects are refused where the first row of a report stops NCBI
/// (`TabularSpec::may_stop_at_failed_lookup`).
#[allow(clippy::too_many_arguments)]
fn check_tabular_lookups(
    formats: &mut [crate::api::local_blast::FormatOutput<'_>],
    format_plans: &[BlastpReport],
    hit_lists: &[Option<BlastpHitList>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>],
    loads: &crate::report::db_subject::HitLoads,
    listed_titles: Option<&database::ListedTitles<'_>>,
    subject_records: &[FastaRecord],
    application: &str,
) -> Result<()> {
    for (format, outfmt) in formats
        .iter_mut()
        .zip(format_plans)
        .filter(|(_, outfmt)| outfmt.is_tabular())
    {
        let spec = outfmt.tabular_spec();
        if !crate::algorithm::common::tabular::checks_failed_lookups(
            spec,
            subject_ids,
            parsed,
            loads,
        ) {
            continue;
        }
        let checked = crate::report::tabular::check_failed_lookups(
            spec,
            hit_lists.iter().flatten().map(|hit_list| {
                let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
                hit_list.hsplist_array[..count]
                    .iter()
                    .filter_map(|hsp_list| hsp_list.hsps.first())
                    .map(|hsp| hsp.s_idx as usize)
            }),
            |s_idx| {
                crate::algorithm::common::tabular::subject_lookup_failure(
                    spec,
                    subject_ids,
                    parsed,
                    loads,
                    listed_titles.is_some(),
                    true,
                    s_idx,
                )
            },
            |s_idx| match listed_titles {
                Some(titles) => crate::report::tabular::failed_database_subject(
                    titles.oid(s_idx),
                    loads.thrown(s_idx),
                ),
                None => crate::report::tabular::failed_record_title(
                    s_idx,
                    subject_records[s_idx].shown_id(),
                ),
            },
            application,
        );
        format.refusal.check(checked)?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
/// The default tabular reports (outfmt 6 and 7 with the `std` fields and a tab), from the final
/// hit lists; the others are written by `tabular::write_blastp_tabular`.
fn write_blastp_tabular_hit_lists(
    hit_lists: &[Option<BlastpHitList>],
    outfmt: OutputFormat,
    writer: &mut impl Write,
    query_ids: &[Arc<[u8]>],
    query_titles: &[Arc<[u8]>],
    subject_ids: &[crate::algorithm::common::tabular::SubjectId],
    subject_accessions: &[Option<[Arc<[u8]>; 2]>],
    loads: &crate::report::db_subject::HitLoads,
    query_batch_skipped: &[bool],
    hits_found: Option<&[usize]>,
    context: &ReportContext,
    mut probe: Option<&mut FormatProbe<'_>>,
    mut warnings: Option<&mut crate::report::query_warnings::QueryWarnings<'_>>,
    epilog: bool,
    queries_before: usize,
    iteration: Option<u32>,
    query_other_coding: Option<i32>,
) -> Result<()> {
    let config = OutputConfig::ncbi_compat();
    // Index of the next HSP in the final hit list (collect_hits_from_hit_lists order).
    let mut hsp_index = 0;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1264-1283,1322-1325
    for (q_idx, hit_list_opt) in hit_lists.iter().enumerate() {
        // The query's warnings come before its lines (`QueryWarnings`).
        if let Some(warnings) = warnings.as_deref_mut() {
            warnings.before_query(q_idx, &mut *writer)?;
        }
        if outfmt == OutputFormat::TabularWithComments {
            write_blastp_outfmt7_header(
                writer,
                context,
                &query_titles[q_idx],
                &crate::report::tabular::STD_FIELDS,
                // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:763-764
                // The formatter splits an ungapped search's Seq-align of a subject into one for
                // each HSP before the header counts them.
                hits_found.map_or_else(
                    || blastp_seq_aligns(std::slice::from_ref(hit_list_opt), true),
                    |found| found[q_idx],
                ),
                query_batch_skipped[q_idx],
                iteration,
            )?;
        }
        let Some(hit_list) = hit_list_opt else {
            continue;
        };
        // The query's formatter (`report::tabular::TabularRows`): the subject that a row whose
        // load throws keeps (`TabularRows::write_thrown_default`).
        let mut rows = crate::report::tabular::TabularRows::default();
        for hsp_list in hit_list.hsplist_array.iter().take(hit_list.hsplist_count) {
            // The subject's load throws where the query's batch has its prefetch off
            // (`report::db_subject::HitLoads`): each row is written as the catch leaves it.
            let thrown = hsp_list
                .hsps
                .first()
                .is_some_and(|first| loads.thrown(first.s_idx as usize).is_some());
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:837-838
            // The subject's first row (the default fields read it) loads it where the query's
            // batch did not prefetch it (`QueryWarnings::subject_read`).
            if let (Some(first), Some(warnings), false) =
                (hsp_list.hsps.first(), warnings.as_deref_mut(), thrown)
            {
                warnings.subject_read(q_idx, first.s_idx, &mut *writer)?;
            }
            for hsp in &hsp_list.hsps {
                // The default fields read the letters of the query (`pident`, `mismatch`): a query
                // whose molecule is `other` ends the run at the first row.
                if let Some(coding) = query_other_coding {
                    return Err(crate::report::tabular::query_sequence_end(
                        writer,
                        &[],
                        coding,
                    ));
                }
                // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1100-1108
                // One row per HSP; the probe marks the row without changing its bytes.
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.begin(hsp_index);
                }
                // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
                let query_id = query_ids
                    .get(hsp.q_idx as usize)
                    .map_or(&b""[..], |id| id.as_ref());
                if thrown {
                    rows.write_thrown_default(
                        writer,
                        query_id,
                        |kept| {
                            SubjectIdFields::of(subject_ids, subject_accessions, kept)
                                .accession_version
                        },
                        hsp.e_value,
                        hsp.bit_score,
                        &config,
                    )?;
                    if let Some(probe) = probe.as_mut() {
                        writer.flush()?;
                        probe.end(hsp_index);
                    }
                    hsp_index += 1;
                    continue;
                }
                rows.found(hsp.s_idx as usize);
                let subject_id =
                    SubjectIdFields::of(subject_ids, subject_accessions, hsp.s_idx as usize);
                crate::report::outfmt6::write_hit_fields_bytes(
                    writer,
                    query_id,
                    // The default fields' subject is `saccver`.
                    subject_id.accession_version,
                    hsp.identity,
                    hsp.num_ident,
                    hsp.length,
                    hsp.mismatch,
                    hsp.gapopen,
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.e_value,
                    hsp.bit_score,
                    &config,
                )?;
                if let Some(probe) = probe.as_mut() {
                    writer.flush()?;
                    probe.end(hsp_index);
                }
                hsp_index += 1;
            }
        }
    }

    // The warnings after the last query's report (`QueryWarnings::after_last_query`).
    if let Some(warnings) = warnings {
        warnings.after_last_query(hit_lists.len(), &mut *writer)?;
    }
    // A search that stops at a query batch (`Empty CBlastQueryVector` with `-query_loc`)
    // writes no epilog (blastp_app.cpp:294-295, `PrintEpilog`). The count is every query
    // formatted (NCBI's `m_QueriesFormatted`): EXPERIMENT (LOSAT_X_BPQBLOCK) the last block of
    // the block search adds the queries of the blocks before it (`queries_before`).
    if outfmt == OutputFormat::TabularWithComments && epilog {
        writeln!(
            writer,
            "# BLAST processed {} queries",
            queries_before + query_titles.len()
        )?;
    }
    Ok(())
}

fn decode_alignment(seq: &[u8], start: usize, end: usize) -> Option<String> {
    if start > end || end > seq.len() {
        return None;
    }
    Some(
        seq[start..end]
            .iter()
            .map(|&aa| ncbistdaa_to_ascii(aa))
            .collect(),
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:518-523
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:649-653
#[inline]
fn query_nomask_sequence(ctx: &QueryContext) -> &[u8] {
    let query_nomask_full = ctx.aa_seq_nomask.as_deref().unwrap_or(&ctx.aa_seq);
    &query_nomask_full[1..query_nomask_full.len() - 1]
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1554
// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1547-1556
// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:936-951
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:4305-4318
// NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:2518-2521
/// The outfmt 0 rows `hits` (the report's copies; the library API's records keep their letters)
/// with their masks in lower case (`HspMasks::lower_row`): a query row's SEG masks of
/// `query_masks` (`query_record_masks`), a subject row's runs of each list of `subject_masks`
/// (the record's lower-case runs with `-lcase_masking`, `x_soft_subject::x_subject_masks`; a
/// masked database search's masks, `db_mask_search::shown_db_subject_masks`), at the columns of
/// `x_FillLocList` (`query_anchored::fill_row_loc_list`: a mask position outside the aligned range
/// is the row's first or last residue), which `x_OutputSeq` puts in lower case when it writes the
/// row. The lists are not kept (`query_anchored::output_seq`): one `HspMasks` serves every row,
/// except under `-html`, where each HSP's go to `html_masks` in the order of `hits`
/// (`HspMasks::copy_rows`, `report::pairwise::HtmlRowMasks`: their spans).
fn set_row_masks(
    hits: &mut [PairwiseHit],
    query_masks: Option<&[Vec<(i32, i32)>]>,
    subject_masks: [&std::collections::HashMap<(u32, u32), Vec<(i32, i32)>>; 2],
    mut html_masks: Option<&mut Vec<crate::report::query_anchored::HspMasks>>,
) {
    let mut masks = crate::report::query_anchored::HspMasks::default();
    for hit in hits {
        if html_masks.is_some() {
            masks.clear();
        }
        let query = query_masks.map_or(&[][..], |masks| &masks[hit.hit.q_idx as usize][..]);
        // A row without letters (`-ungapped`'s HSP of length 0, whose range ends before it
        // starts) has no mask.
        if let Some(shown) = hit
            .query_seq
            .as_mut()
            .filter(|row| !query.is_empty() && !row.is_empty())
        {
            let mut letters = std::mem::take(shown).into_bytes();
            let row = &letters;
            masks.set_row(
                0,
                0,
                row,
                // A `-` is a residue at a column of `query_stored_gaps` (a PSSM query's gap
                // residue).
                |column| {
                    row[column] != b'-'
                        || hit
                            .query_stored_gaps
                            .binary_search(&(column as u32))
                            .is_ok()
                },
                (hit.hit.q_start as u32 - 1, hit.hit.q_end as u32 - 1),
                true,
                1,
                query.iter().map(|&(from, to)| (from as u32, to as u32, 0)),
            );
            masks.lower_row(0, &mut letters);
            *shown = String::from_utf8(letters).expect("ASCII alignment row");
        }
        let key = (hit.hit.q_idx, hit.hit.s_idx);
        let mut runs = subject_masks
            .iter()
            .filter_map(|masks| masks.get(&key))
            .flatten()
            .peekable();
        if let Some(shown) = hit
            .subject_seq
            .as_mut()
            .filter(|row| !row.is_empty() && runs.peek().is_some())
        {
            let mut letters = std::mem::take(shown).into_bytes();
            let row = &letters;
            // A stored `-` is a residue of the subject (`PairwiseHit::subject_stored_gaps`).
            let stored = &hit.subject_stored_gaps;
            masks.set_row(
                1,
                0,
                row,
                |column| row[column] != b'-' || stored.contains(&(column as u32)),
                (hit.hit.s_start as u32 - 1, hit.hit.s_end as u32 - 1),
                true,
                1,
                runs.map(|&(from, to)| (from as u32, to as u32, 0)),
            );
            masks.lower_row(1, &mut letters);
            *shown = String::from_utf8(letters).expect("ASCII alignment row");
        }
        if let Some(table) = html_masks.as_deref_mut() {
            table.push(masks.copy_rows());
        }
    }
}

/// The query masks of the search's results, which every report of the call reads (outfmt 0's
/// lower case, BLAST XML's masked letters, the archive's masks, outfmt 0's refusal of a query
/// that is also a subject): `query_record_masks`, or none for each query of a PSSM search
/// (`report::record::carries_query_masks` of `app`'s program type: psiblast in every round).
fn result_query_masks(
    app: crate::algorithm::blastp::args::BlastpApp,
    query_frames: &[Vec<QueryFrame>],
    offsets: &seq_range::Placements,
) -> Vec<Vec<(i32, i32)>> {
    if !crate::report::record::carries_query_masks(app.program_type()) {
        return vec![Vec::new(); query_frames.len()];
    }
    query_record_masks(query_frames, offsets)
}

/// Outfmt 0's refusal of a query that is also a subject and has masks, each outfmt 0 output's
/// (`believed::refuse_shared_masks_in`; the query-anchored sections draw these masks as NCBI
/// does, `query_anchored::fill_loc_list`). A PSSM search's results carry no query masks
/// (`report::record::carries_query_masks`), so neither arm of the refusal applies to psiblast:
/// its Sbjct rows show the subject's masks on the subject's rows only (`SetSubjectMasks`, and only
/// query masks are drawn by id).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1583-1584
fn refuse_shared_record_masks(
    args: &ResolvedBlastpArgs,
    format_plans: &[BlastpReport],
    outputs: &mut ReportOutputs<'_>,
    query_records: &[FastaRecord],
    query_frames: &[Vec<QueryFrame>],
    offsets: &seq_range::Placements,
) -> Result<()> {
    let pairwise = |outfmt: &BlastpReport| {
        matches!(
            outfmt,
            BlastpReport::Pairwise(crate::report::pairwise::PairwiseView::Pairwise)
        )
    };
    if !args.parse_deflines
        || !crate::report::record::carries_query_masks(args.app.program_type())
        || !format_plans.iter().any(pairwise)
    {
        return Ok(());
    }
    let masks = query_record_masks(query_frames, offsets);
    crate::blastinput::fasta_reader::believed::refuse_shared_masks_in(
        outputs
            .formats
            .iter_mut()
            .zip(format_plans)
            .filter(|(_, outfmt)| pairwise(outfmt))
            .map(|(format, _)| &mut format.refusal),
        query_records,
        |query| masks.get(query).is_some_and(|m| !m.is_empty()),
        args.lcase_masking,
        args.app.name(),
    )
}

/// The masked regions of each query as NCBI's `TMaskedQueryRegions` holds them (`Map`,
/// api/blast_aux.cpp:825-842): the SEG masks merged, a mask that covers the whole query
/// dropped, and the rest moved into the query's record by the start of its interval
/// (`offsets`, with `-query_loc`, as the hits' coordinates); inclusive 0-based positions.
fn query_record_masks(
    query_frames: &[Vec<QueryFrame>],
    offsets: &seq_range::Placements,
) -> Vec<Vec<(i32, i32)>> {
    query_frames
        .iter()
        .map(|frames| {
            let Some(frame) = frames.first() else {
                return Vec::new();
            };
            let ranges = frame
                .seg_masks
                .iter()
                .map(|&(start, end)| (start as i32, end as i32 - 1))
                .collect();
            crate::algorithm::common::masks::combine_masks(ranges)
                .into_iter()
                .filter(|&(from, to)| from <= to && !(from == 0 && to == frame.aa_len as i32 - 1))
                .collect()
        })
        .enumerate()
        .map(|(q_idx, masks): (usize, Vec<(i32, i32)>)| {
            let offset = offsets.offset(q_idx) as i32;
            masks
                .into_iter()
                .map(|(from, to)| (from + offset, to + offset))
                .collect()
        })
        .collect()
}

/// Moves a hit's reported coordinates from the searched letters to the records, by the
/// start of the query's and the subject's interval, and sets its query length to the
/// `-query_loc` range as given (the tabular `qlen`). Without ranges nothing changes.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1630-1637
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:884-885
fn shift_to_records(hit: &mut Hit, ranges: &BlastpRanges<'_>) {
    let q_shift = ranges.query_input.placements.offset(hit.q_idx as usize);
    let s_shift = ranges.subject_placements.offset(hit.s_idx as usize);
    hit.q_start += q_shift;
    hit.q_end += q_shift;
    hit.s_start += s_shift;
    hit.s_end += s_shift;
    if let Some(length) = ranges.query_input.range_length {
        hit.query_length = length;
    }
}

/// `shift_to_records` for an HSP of the hit lists (the tabular rows that are written
/// without the alignment).
fn shift_hsp_to_records(hsp: &mut BlastpHsp, ranges: &BlastpRanges<'_>) {
    let q_shift = ranges.query_input.placements.offset(hsp.q_idx as usize);
    let s_shift = ranges.subject_placements.offset(hsp.s_idx as usize);
    hsp.q_start += q_shift;
    hsp.q_end += q_shift;
    hsp.s_start += s_shift;
    hsp.s_end += s_shift;
    if let Some(length) = ranges.query_input.range_length {
        hsp.query_length = length;
    }
}

/// psiblast's PSSMs after its round (`-save_pssm_after_last_round` with `-out_pssm` and/or
/// `-out_ascii_pssm`): one for each query whose results have alignments, in the order of the
/// queries, from the query's whole record and its Seq-aligns before the formatter's cut
/// (`algorithm::psiblast::pssm`). Their files are made ready here, before the reports, and
/// written after them (`PssmFiles`: the ASN.1 bytes staged beside each target, which is not
/// changed before `commit`). A query whose PSSM fails or whose ASN.1 file cannot be written stops
/// the run after its report (`PssmStop`, handled by `search`), with the PSSMs of the queries
/// before it.
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:147-162
// The query is its record's Bioseq, whole also with `-query_loc`; its Seq-align id is
// `Query_<n>` (`n` its position in the input, as the reports name it); the ancillary data are
// the query's (`CBlastAncillaryData`, blast_results.cpp:71-114: the ungapped block of its
// context and the gapped one, and in the first round the PSI blocks are copies of them).
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn psiblast_pssms(
    args: &ResolvedBlastpArgs,
    request: &crate::algorithm::psiblast::pssm::PssmRequest,
    hit_lists: &[Option<BlastpHitList>],
    ranges: &BlastpRanges<'_>,
    shown_subjects: ShownSubjects<'_>,
    contexts: &[QueryContext],
    gapped: &KarlinParams,
    keys: psi_round::SubjectKeys<'_>,
) -> Result<crate::algorithm::psiblast::pssm::PssmFiles> {
    use crate::algorithm::psiblast::pssm::{
        query_pssm, PssmAlignment, PssmAncillary, PssmFiles, PssmQueryRecord, PssmScoring, PssmStop,
    };
    use crate::core::psi::KarlinBlk;
    let subject = |s_idx: u32| keys.letters(shown_subjects, s_idx);
    let scoring = PssmScoring {
        matrix: args.scoring.matrix,
        matrix_name: &args.matrix_name,
        gap_open: args.scoring.gap_open,
        gap_extend: args.scoring.gap_extend,
    };
    let gapped = KarlinBlk {
        lambda: gapped.lambda,
        k: gapped.k,
        h: gapped.h,
    };
    // The rows of each query's PSSM are the handles of the Seq-align ids that it takes
    // (`SubjectKeys`): keys that NCBI's handles can take for one, and giim ids, are refused before
    // any PSSM is made (`rounds::HandleIndex`).
    for hit_list in hit_lists.iter().flatten() {
        let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
        let rows: std::collections::BTreeSet<crate::algorithm::psiblast::rounds::SetKey> = hit_list
            .hsplist_array[..count]
            .iter()
            .flat_map(|hsp_list| &hsp_list.hsps)
            .filter(|hsp| psi_round::included(hsp.e_value, args.inclusion_ethresh))
            .map(|hsp| keys.key(hsp.s_idx))
            .collect();
        crate::algorithm::psiblast::rounds::SetKey::check_handles(&rows)?;
    }
    let mut files = PssmFiles::new(request);
    for (query, hit_list) in hit_lists.iter().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
        let q_shift = ranges.query_input.placements.offset(query);
        let alignments: Vec<PssmAlignment<'_>> = hit_list.hsplist_array[..count]
            .iter()
            .flat_map(|hsp_list| &hsp_list.hsps)
            .map(|hsp| PssmAlignment {
                subject: hsp.s_idx,
                e_value: crate::report::serial::seqalign::seqalign_evalue(hsp.e_value),
                q_start: hsp.q_start - 1 + q_shift,
                s_start: hsp.s_start - 1,
                q_length: hsp.q_end + 1 - hsp.q_start,
                ops: hsp.gap_info.as_deref(),
            })
            .collect();
        if alignments.is_empty() {
            continue;
        }
        let ordinal = ranges.query_input.ordinal(query);
        let record = &ranges.input_records[ordinal];
        let id = crate::algorithm::psiblast::pssm::fasta_query_id(
            record,
            args.first_query_ordinal + ranges.x_ordinal_base + ordinal + 1,
        );
        let ungapped = &contexts[query].karlin_params;
        let ungapped = KarlinBlk {
            lambda: ungapped.lambda,
            k: ungapped.k,
            h: ungapped.h,
        };
        let ancillary = PssmAncillary {
            ungapped,
            gapped,
            psi_ungapped: ungapped,
            psi_gapped: gapped,
        };
        let mut unread = Vec::new();
        let made = query_pssm(
            request,
            scoring,
            gapped,
            args.inclusion_ethresh,
            PssmQueryRecord {
                letters: record.seq(),
                id: &id,
                title: record.title(),
            },
            &alignments,
            &subject,
            &|s_idx| keys.key(s_idx),
            ancillary,
            &mut unread,
            // The PSSM after the last round, which no archive holds.
            false,
        );
        let mut warnings = crate::algorithm::psiblast::pssm::PssmReadMessages::default();
        keys.unread_warnings(args.app.ncbi_name(), &unread, &mut warnings);
        files.add_read_warnings(query, warnings);
        // NCBI's posts of the bytes replaced in the file come as it is written, after the
        // PSSM's read messages (`PssmFiles::take_read_warnings`), and before a write fault of the
        // file, which NCBI meets at the first flush of its 4 KiB buffer (`PssmFiles::add`), as the
        // rounds write them (`psi_round::write_round_pssm`).
        let made = made.and_then(|pssm| {
            files.add_file_posts(
                query,
                crate::algorithm::psiblast::pssm::file_posts(args.app.ncbi_name(), &pssm.asn_posts),
            );
            files.add(pssm.asn.as_deref(), pssm.ascii)
        });
        if let Err(error) = made {
            return Err(PssmStop {
                query,
                error,
                files,
            }
            .into());
        }
    }
    Ok(files)
}

/// The stored residue of pyrrolysine (O): the only code that `ncbistdaa_to_ascii` shows as O
/// (`NCBISTDAA_TO_AMINOACID`, blast_encoding.c:115-121, cited above it).
const STORED_O: u8 = 26;

/// The stored gap character (`-`): the code that `ncbistdaa_to_ascii` shows as `-`, the letter
/// of a gap in the rows of an alignment.
const STORED_GAP: u8 = 0;

// No NCBI counterpart: how LOSAT's reports reach the letters they read; NCBI reads each
// sequence from its Bioseq (alnvec.cpp:116-121, cited above `restore_displayed_o`).
/// The letters of one sequence as the reports read them (ASCII, as read).
#[derive(Clone, Copy)]
enum ShownLetters<'a> {
    /// A record's letters.
    Ascii(&'a [u8]),
    /// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): a database subject's stored residues (without the
    /// sentinels), each read as `ncbistdaa_to_ascii` gives it, as
    /// `database::ListedSubjects::read_titles` writes the letters of its records.
    Stored(&'a [u8]),
    /// The letters of a PSSM query's Bioseq where they are not its record's
    /// (`PsiSearch::bioseq_letters`: iupacaa letters with one in lower case), which the reports
    /// show and compare as they are (`report::row_letters::QueryCase::Bioseq`).
    Bioseq(&'a [u8]),
}

impl<'a> ShownLetters<'a> {
    /// A query's letters: those of its Bioseq where they are not its record's (`bioseq`,
    /// `PsiSearch::bioseq_letters`), its record's (`record`) otherwise.
    fn query(record: &'a [u8], bioseq: Option<&'a [u8]>) -> Self {
        bioseq.map_or(Self::Ascii(record), Self::Bioseq)
    }

    /// Whether a letter is a pyrrolysine (O).
    fn has_o(self) -> bool {
        match self {
            // Two byte searches (memchr) instead of one case-folding loop: the same answer.
            Self::Ascii(letters) | Self::Bioseq(letters) => {
                letters.contains(&b'O') || letters.contains(&b'o')
            }
            Self::Stored(residues) => residues.contains(&STORED_O),
        }
    }

    /// Whether the reports show a letter that the search did not read: an O, which it read as X
    /// (`displayed_letter`), or a Bioseq's letter (`Bioseq`).
    fn shows_other_letters(self) -> bool {
        matches!(self, Self::Bioseq(_)) || self.has_o()
    }

    /// Whether the recount (`recount_shown_identities`) counts the HSPs of these letters again:
    /// the reports show a letter that the search did not read (`shows_other_letters`), or a gap
    /// character, a residue that shows as the letter of a gap (`-`: makeblastdb stores it in a
    /// protein database, a PSSM query holds it as the gap residue; a `-` of a FASTA input is
    /// dropped as NCBI drops it). Each test is a byte search of the letters (`contains`).
    fn recounted(self) -> bool {
        match self {
            Self::Ascii(letters) => letters.contains(&b'-') || self.has_o(),
            Self::Stored(residues) => {
                residues.contains(&STORED_GAP) || residues.contains(&STORED_O)
            }
            Self::Bioseq(_) => true,
        }
    }

    /// The letter at `index` (0-based).
    fn get(self, index: usize) -> Option<u8> {
        match self {
            Self::Ascii(letters) | Self::Bioseq(letters) => letters.get(index).copied(),
            Self::Stored(residues) => residues
                .get(index)
                .map(|&residue| ncbistdaa_to_ascii(residue) as u8),
        }
    }

    /// The letter that the reports show at `index` (0-based) where the search read `letter`
    /// (`ncbistdaa_to_ascii`): a Bioseq's letter as it is, otherwise `displayed_letter`.
    fn shown(self, letter: char, index: usize) -> char {
        match self {
            Self::Bioseq(letters) => letters
                .get(index)
                .map_or(letter, |&shown| char::from(shown)),
            _ => displayed_letter(letter, self.get(index)),
        }
    }

    /// How the reports compare a row of these letters (a query's).
    fn case(self) -> crate::report::row_letters::QueryCase {
        match self {
            Self::Bioseq(_) => crate::report::row_letters::QueryCase::Bioseq,
            Self::Ascii(_) | Self::Stored(_) => crate::report::row_letters::QueryCase::Masked,
        }
    }
}

/// The 0-based position of the first residue of `side` in the protein alignment `hit`.
fn row_start(hit: &Hit, side: RowSide) -> usize {
    match side {
        RowSide::Query => hit.q_start.saturating_sub(1),
        RowSide::Subject => hit.s_start.min(hit.s_end).saturating_sub(1),
    }
}

// No NCBI counterpart: how LOSAT's reports reach the queries' letters.
/// The queries whose letters the reports read, by query index: their records, and the letters of
/// a PSSM query's Bioseq where they are not its record's (`bioseq`, `PsiSearch::bioseq_letters`;
/// a run of such a query has no other query).
#[derive(Clone, Copy)]
struct ShownQueries<'a> {
    records: &'a [FastaRecord],
    bioseq: Option<&'a [u8]>,
}

impl<'a> ShownQueries<'a> {
    /// The queries `records` of the run `args`.
    fn of(records: &'a [FastaRecord], args: &'a ResolvedBlastpArgs) -> Self {
        let bioseq = args
            .psi_search
            .as_deref()
            .and_then(|psi| psi.bioseq_letters.as_deref());
        Self { records, bioseq }
    }

    /// The record letters of the query `q_idx` (none past the queries).
    fn record(self, q_idx: usize) -> &'a [u8] {
        self.records
            .get(q_idx)
            .map_or(&[][..], |record| record.seq())
    }

    /// The letters of the query `q_idx` as the reports read them (`ShownLetters::query`).
    fn letters(self, q_idx: usize) -> Option<ShownLetters<'a>> {
        self.records
            .get(q_idx)
            .map(|record| ShownLetters::query(record.seq(), self.bioseq))
    }
}

// No NCBI counterpart: how LOSAT's reports reach the subjects' letters and titles.
/// The subjects whose letters and titles the reports read, by subject index.
#[derive(Clone, Copy)]
enum ShownSubjects<'a> {
    /// The records of a subject set, or of a database's listed subjects
    /// (`database::ListedTitles::records`).
    Records(&'a [FastaRecord]),
    /// EXPERIMENT (LOSAT_X_BPLISTEDLEAN): a database's listed subjects, by rank: their stored
    /// residues (`database::ListedSubjects::encoded`) and the titles of the reported ones
    /// (`database::ListedTitles::titles`).
    Listed {
        subjects: &'a [EncodedProtein],
        titles: &'a [Arc<[u8]>],
    },
}

impl<'a> ShownSubjects<'a> {
    /// The letters of the subject `s_idx`.
    fn letters(self, s_idx: usize) -> Option<ShownLetters<'a>> {
        match self {
            Self::Records(records) => records
                .get(s_idx)
                .map(|record| ShownLetters::Ascii(record.seq())),
            // A subject without residues has no O and no HSP with
            // a query that has one; no report reads its letters (`ListedSubjects::x_read_shown`).
            Self::Listed { subjects, .. } => subjects.get(s_idx).map(|subject| {
                ShownLetters::Stored(subject.aa_seq.get(1..=subject.aa_len).unwrap_or_default())
            }),
        }
    }

    /// The title of the subject `s_idx`.
    fn title(self, s_idx: usize) -> Option<&'a [u8]> {
        match self {
            Self::Records(records) => records.get(s_idx).map(|record| &record.title[..]),
            Self::Listed { titles, .. } => titles.get(s_idx).map(|title| &title[..]),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/alnmgr/alnvec.cpp:116-121
/// The displayed row `row` of the `side` of the alignment `hit`, whose letters are `sequence`:
/// the reports read the input sequence, so a pyrrolysine (O) that the search read as X is shown
/// as O, and a PSSM query's Bioseq letters as they are (`ShownLetters::shown`). The columns of
/// the residues are those of the alignment's operations (`RowSide::reads`), not the row's
/// letters other than `-`: a gap character is a residue too (`recounted`).
fn restore_displayed_o(
    row: Option<String>,
    sequence: ShownLetters<'_>,
    hit: &Hit,
    side: RowSide,
) -> Option<String> {
    let row = row?;
    if !sequence.shows_other_letters() {
        return Some(row);
    }
    let implicit = [GapEditOp::Sub(hit.length as u32)];
    let ops = hit.gap_info.as_deref().unwrap_or(&implicit);
    let mut position = row_start(hit, side);
    let mut columns = row.chars();
    let mut shown = String::with_capacity(row.len());
    for op in ops {
        let reads = side.reads(op);
        for column in columns.by_ref().take(op.num() as usize) {
            if reads {
                shown.push(sequence.shown(column, position));
                position += 1;
            } else {
                shown.push(column);
            }
        }
    }
    shown.extend(columns);
    Some(shown)
}

/// The letter that the reports show for the searched letter `letter` (`ncbistdaa_to_ascii`) where
/// the input sequence holds `input` (`ShownLetters::get`): a pyrrolysine (O) that the search read
/// as X is shown as O (`restore_displayed_o`, the query-anchored sections' letters).
#[inline]
fn displayed_letter(letter: char, input: Option<u8>) -> char {
    match input {
        Some(residue) if residue.eq_ignore_ascii_case(&b'O') && letter == 'X' => 'O',
        _ => letter,
    }
}

// No NCBI counterpart: reads the LOSAT_X_BPRECOUNTO switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPRECOUNTO): strict set (on unless `0`). `recount_shown_identities` tests
/// each query's letters for an O, and each subject's for an O or a gap character, once per call, not
/// once per HSP: a sequence's test reads all its letters, and a query with many HSPs, or a subject
/// in many of them, was read again for each (line BPF s12: about 0.15 s of the human x yeast
/// `-task blastp-fast` search at 16 threads, on the search's thread between the blocks, while the
/// pool waited).
fn x_bprecounto() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPRECOUNTO"))
}

// No NCBI counterpart: remembers a test's answer (LOSAT_X_BPRECOUNTO); it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPRECOUNTO): `test()` for the sequence `index`, run the first time that index
/// is asked and remembered in `memo` (0 = not tested, 1 = false, 2 = true). The tests read
/// sequences that do not change during the call, so each answer is the one a new test would give.
fn x_tested_once(memo: &mut Vec<u8>, index: usize, test: impl FnOnce() -> bool) -> bool {
    if memo.len() <= index {
        memo.resize(index + 1, 0);
    }
    if memo[index] == 0 {
        memo[index] = 1 + u8::from(test());
    }
    memo[index] == 2
}

/// Recounts the identities, mismatches and positives of the HSPs whose query or subject has
/// pyrrolysine (O) or a gap character (`-`, stored in a protein database or a PSSM query's gap
/// residue), or whose query shows its Bioseq's letters (`ShownQueries::bioseq`, the letters of a
/// PSSM query's Bioseq where they are not its record's: `PsiSearch::bioseq_letters`), and of every
/// HSP of a search that read other letters than the query
/// record's (`searched_other_letters`: psiblast's later rounds, whose engine counts on the letters
/// that the rounds before masked, `psi_round::RoundQuery`); `ShownLetters::recounted`.
/// The search reads O as X (`encode_protein_sequence`) and a lower-case letter of a PSSM query as
/// its capital, but NCBI's reports compare the letters of the alignment as read, where O is
/// identical only to O (and O and X are not identical), and a lower-case letter to no subject
/// letter (`report::row_letters::QueryCase::Bioseq`).
/// A gap character shows as the letter of a gap, so a column of a gap in one sequence and that
/// residue in the other is an identity (and a gap), and the mismatches (the columns less the
/// identities and the gaps) can be fewer than none.
///
/// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1006-1020
///
/// The Seq-aligns (outfmt 8, 12) show the engine's counts, NCBI's `num_ident` and
/// `num_positives` of the search's letters (blast_seqalign.cpp:1199-1222), and the outfmt 0
/// report reads that `num_ident` back for `-sorthits` and `-sorthsps`
/// (`PairwiseHit::seqalign_num_ident`): the counts of each HSP that the recount changes are
/// kept in `engine_counts` first.
fn recount_shown_identities(
    hit_lists: &mut [Option<BlastpHitList>],
    queries: ShownQueries<'_>,
    subjects: ShownSubjects<'_>,
    matrix: ScoringMatrix,
    searched_other_letters: bool,
    engine_counts: &mut super::seqalign::EngineCounts,
) {
    // EXPERIMENT (LOSAT_X_BPRECOUNTO): the answers of the tests by query and by subject index
    // (`ShownLetters::recounted`, `x_tested_once`).
    let x_once = x_bprecounto();
    let mut x_query_o = Vec::new();
    let mut x_subject_o = Vec::new();
    for (q_idx, hit_list) in hit_lists.iter_mut().enumerate() {
        let Some(hit_list) = hit_list else {
            continue;
        };
        let count = hit_list.hsplist_count.min(hit_list.hsplist_array.len());
        for (list_index, hsp_index, hsp) in hit_list.hsplist_array[..count]
            .iter_mut()
            .enumerate()
            .flat_map(|(list_index, hsp_list)| {
                hsp_list
                    .hsps
                    .iter_mut()
                    .enumerate()
                    .map(move |(hsp_index, hsp)| (list_index, hsp_index, hsp))
            })
        {
            let (Some(query), Some(subject)) = (
                queries.letters(hsp.q_idx as usize),
                subjects.letters(hsp.s_idx as usize),
            ) else {
                continue;
            };
            // LOSAT_X_BPRECOUNTO: each query's and each subject's letters read once per call.
            let recounted = searched_other_letters
                || if x_once {
                    x_tested_once(&mut x_query_o, hsp.q_idx as usize, || query.recounted())
                        || x_tested_once(&mut x_subject_o, hsp.s_idx as usize, || {
                            subject.recounted()
                        })
                } else {
                    query.recounted() || subject.recounted()
                };
            if !recounted {
                continue;
            }
            let implicit = [GapEditOp::Sub(hsp.length as u32)];
            let ops = hsp.gap_info.as_deref().unwrap_or(&implicit);
            let mut query_row = String::new();
            let mut subject_row = String::new();
            let mut q_index = hsp.q_start.saturating_sub(1);
            let mut s_index = hsp.s_start.saturating_sub(1);
            let mut complete = true;
            let compared = query.case().compared_letters();
            for op in ops {
                // A deletion is a gap in the query, an insertion a gap in the subject.
                let (in_query, in_subject) = match *op {
                    GapEditOp::Sub(_) => (true, true),
                    GapEditOp::Del(_) => (false, true),
                    GapEditOp::Ins(_) => (true, false),
                };
                for _ in 0..op.num() as usize {
                    let q = if in_query {
                        query.get(q_index)
                    } else {
                        Some(b'-')
                    };
                    let s = if in_subject {
                        subject.get(s_index)
                    } else {
                        Some(b'-')
                    };
                    let (Some(q), Some(s)) = (q, s) else {
                        complete = false;
                        break;
                    };
                    query_row.push(compared[usize::from(q)] as char);
                    subject_row.push(s.to_ascii_uppercase() as char);
                    q_index += usize::from(in_query);
                    s_index += usize::from(in_subject);
                }
            }
            if !complete || query_row.len() != hsp.length {
                continue;
            }
            let num_ident = query_row
                .bytes()
                .zip(subject_row.bytes())
                .filter(|&(q, s)| q == s)
                .count();
            // The gaps of the Dense-seg (`GetAlignLengths`), not the `-` letters of the rows.
            let gap_letters: usize = ops
                .iter()
                .filter(|op| !op.is_sub())
                .map(|op| op.num() as usize)
                .sum();
            let num_positives = crate::report::row_letters::display_positives(
                query_row.as_bytes(),
                subject_row.as_bytes(),
                matrix,
                query.case(),
            );
            if (num_ident, num_positives) != (hsp.num_ident, hsp.num_positives) {
                engine_counts.save(
                    (q_idx as u32, list_index as u32, hsp_index as u32),
                    (hsp.num_ident, hsp.num_positives),
                );
            }
            hsp.num_ident = num_ident;
            hsp.num_positives = num_positives;
            hsp.mismatch = hsp.length as i64 - (num_ident + gap_letters) as i64;
            hsp.identity = if hsp.length > 0 {
                (num_ident as f64 * 100.0) / hsp.length as f64
            } else {
                0.0
            };
        }
    }
}

/// The columns (0-based) of the rows of `hit` whose letter on `side` is a gap character, a
/// residue `gap` of `sequence` (`row_letters::stored_gap_columns`): the subject's stored residues
/// and `STORED_GAP`, or the query record's letters and `-` (a PSSM query's gap residue).
fn stored_gap_columns(sequence: &[u8], gap: u8, hit: &Hit, side: RowSide) -> Vec<u32> {
    let implicit = [GapEditOp::Sub(hit.length as u32)];
    let ops = hit.gap_info.as_deref().unwrap_or(&implicit);
    crate::report::row_letters::stored_gap_columns(sequence, gap, row_start(hit, side), ops, side)
}

/// The rows of the final HSPs `hits`, which show the letters of the queries as `queries` reads them
/// (a PSSM query's Bioseq letters where they are not its record's: `ShownQueries::letters`).
#[allow(clippy::too_many_arguments)]
fn build_pairwise_hits(
    hits: Vec<(Hit, u8, usize, i32)>,
    query_contexts: &[QueryContext],
    subjects: &[EncodedProtein],
    matrix: ScoringMatrix,
    render_alignment: bool,
    queries: ShownQueries<'_>,
    shown_subjects: ShownSubjects<'_>,
    x_pool: Option<&crate::utils::threading::SearchPool<'_>>,
) -> Result<Vec<PairwiseHit>> {
    let row = |(mut hit, comp_adjustment_method, seqalign_num_ident, num): (
        Hit,
        u8,
        usize,
        i32,
    )|
     -> Result<PairwiseHit> {
        let q_idx = hit.q_idx as usize;
        let s_idx = hit.s_idx as usize;
        let query = query_nomask_sequence(&query_contexts[q_idx]);
        let subject = &subjects[s_idx].aa_seq[1..subjects[s_idx].aa_seq.len() - 1];

        let (query_seq, subject_seq) = if render_alignment {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
            let rendered = build_alignment_view_with_matrix(query, subject, &hit, matrix);
            if let Some(view) = rendered {
                (Some(view.query_aln), Some(view.subject_aln))
            } else {
                let q_start = hit.q_start.saturating_sub(1);
                let q_end = hit.q_end.min(query.len());
                let s_start = hit.s_start.saturating_sub(1);
                let s_end = hit.s_end.min(subject.len());
                let ungapped_span = q_end
                    .saturating_sub(q_start)
                    .min(s_end.saturating_sub(s_start));
                if ungapped_span == hit.length {
                    (
                        decode_alignment(query, q_start, q_start + ungapped_span),
                        decode_alignment(subject, s_start, s_start + ungapped_span),
                    )
                } else {
                    (None, None)
                }
            }
        } else {
            (None, None)
        };
        let query_record = queries.record(q_idx);
        let query_letters = ShownLetters::query(query_record, queries.bioseq);
        let query_seq = restore_displayed_o(query_seq, query_letters, &hit, RowSide::Query);
        let subject_seq = restore_displayed_o(
            subject_seq,
            shown_subjects
                .letters(s_idx)
                .unwrap_or(ShownLetters::Ascii(&[])),
            &hit,
            RowSide::Subject,
        );
        // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:996-1028
        let positives = if render_alignment {
            let (Some(query), Some(subject)) = (&query_seq, &subject_seq) else {
                anyhow::bail!("unimplemented BLASTP report alignment for final HSP");
            };
            crate::report::row_letters::display_positives(
                query.as_bytes(),
                subject.as_bytes(),
                matrix,
                query_letters.case(),
            )
        } else {
            hit.num_positives
        };
        hit.num_positives = positives;
        let gaps = hit.gap_letters();
        let [query_stored_gaps, subject_stored_gaps] = if subject_seq.is_some() {
            [
                stored_gap_columns(query_record, b'-', &hit, RowSide::Query),
                stored_gap_columns(subject, STORED_GAP, &hit, RowSide::Subject),
            ]
        } else {
            [Vec::new(), Vec::new()]
        };

        Ok(PairwiseHit {
            hit,
            query_seq,
            subject_seq,
            // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1049-1093
            query_frame: Some(1),
            subject_frame: Some(1),
            positives: Some(positives),
            gaps: Some(gaps),
            subject_length: Some(subjects[s_idx].aa_len),
            subject_title: shown_subjects
                .title(s_idx)
                .and_then(crate::blastinput::fasta_reader::hit_subject_title),
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:332-342: the HSP's composition adjustment.
            comp_adjust_method: Some(comp_adjustment_method),
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1180-1182
            // The linked set's size with psiblast's `-sum_stats true` (the number in parentheses
            // after the alignment's `Expect`, showalign.cpp:3595-3598); an HSP outside a set has
            // no "sum_n" score.
            sum_n: (num > 1).then_some(num),
            subject_stored_gaps,
            query_stored_gaps,
            query_case: query_letters.case(),
            // The Seq-align's `num_ident` (`collect_hits_and_methods_from_hit_lists`): the
            // engine's, which `hit.num_ident` holds unless the display recount replaced it
            // (`recount_shown_identities`).
            seqalign_num_ident,
        })
    };
    // EXPERIMENT (LOSAT_X_PAIRWISEPAR): the rows of the final list as an indexed map on the
    // search's pool (`report::x_pairwise_par`), the first error in the list order; each row is
    // `row` of its hit, which reads nothing another hit writes.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if let Some(pool) = x_pool {
        use rayon::prelude::*;
        let rows: Vec<Result<PairwiseHit>> =
            pool.install(|| hits.into_par_iter().map(&row).collect());
        let mut kept = Vec::with_capacity(rows.len());
        for row in rows {
            kept.push(row?);
        }
        return Ok(kept);
    }
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    let _ = x_pool;
    // The rows in a vector of the hits' count (a `Result` collect grows its vector from empty,
    // copying the rows of every window at each growth), the first error in the list order.
    let mut rows = Vec::with_capacity(hits.len());
    for hit in hits {
        rows.push(row(hit)?);
    }
    Ok(rows)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:52-61
fn compute_diag_table_shape(query_length: i32, window: i32) -> (i32, i32) {
    let mut diag_array_size = 1i32;
    while diag_array_size < query_length + window {
        diag_array_size <<= 1;
    }

    (diag_array_size, diag_array_size - 1)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3746-3786
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3752-3786
#[cfg(test)]
fn build_restricted_align_array(
    init_hsps: &[InitHSP],
    contexts: &[QueryContext],
    hit_cutoff_scores: &[i32],
    num_queries: usize,
) -> Vec<bool> {
    let mut found = Vec::new();
    let mut restricted = Vec::new();
    build_restricted_align_array_into(
        init_hsps,
        contexts,
        blastp_context_lookup_bounds(contexts),
        hit_cutoff_scores,
        num_queries,
        &mut found,
        &mut restricted,
    );
    restricted
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3748-3786
fn build_restricted_align_array_into(
    init_hsps: &[InitHSP],
    contexts: &[QueryContext],
    context_lookup_bounds: BlastpContextLookupBounds,
    hit_cutoff_scores: &[i32],
    num_queries: usize,
    found: &mut Vec<bool>,
    restricted: &mut Vec<bool>,
) {
    found.clear();
    found.resize(num_queries, false);
    restricted.clear();
    restricted.resize(num_queries, false);
    for init_hsp in init_hsps {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3753-3765
        let context_idx =
            blastp_get_ungapped_hsp_context_with_bounds(contexts, context_lookup_bounds, init_hsp);
        let query_index = usize::try_from(contexts[context_idx].q_idx)
            .expect("NCBI BLAST requires blastp context query indices to fit in usize");
        if query_index >= num_queries || found[query_index] {
            continue;
        }
        found[query_index] = true;
        restricted[query_index] = init_hsp.ungapped_data.score
            < (BLASTP_RESTRICTED_MULT * hit_cutoff_scores[context_idx] as f64) as i32;
    }
}

#[derive(Clone, Copy, Debug, Default)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign_priv.h:153-156
struct BlastpInitHspNode {
    init_index: usize,
    best_score: i32,
    next_index: Option<usize>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3468-3499
fn compare_init_hsps_by_query_offset_score(a: &InitHSP, b: &InitHSP) -> std::cmp::Ordering {
    a.ungapped_data
        .q_start
        .cmp(&b.ungapped_data.q_start)
        .then_with(|| b.ungapped_data.score.cmp(&a.ungapped_data.score))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:274-300
fn compare_init_hsps_by_score_ncbi(a: &InitHSP, b: &InitHSP) -> std::cmp::Ordering {
    b.ungapped_data
        .score
        .cmp(&a.ungapped_data.score)
        .then_with(|| a.ungapped_data.s_start.cmp(&b.ungapped_data.s_start))
        .then_with(|| b.ungapped_data.length.cmp(&a.ungapped_data.length))
        .then_with(|| a.ungapped_data.q_start.cmp(&b.ungapped_data.q_start))
        .then_with(|| b.ungapped_data.length.cmp(&a.ungapped_data.length))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3532-3533
fn ncbi_qsort_init_hsps_by_query_offset_score(init_hsps: &mut [InitHSP]) {
    if init_hsps.len() <= 1 {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3584-3613
    // Stable, as glibc's qsort under the pinned NCBI BLAST+: ungapped alignments
    // with the same query start and score tie, and their order decides which one
    // the chaining links to the other.
    init_hsps.sort_by(compare_init_hsps_by_query_offset_score);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:306-309
fn ncbi_qsort_init_hsps_by_score(init_hsps: &mut [InitHSP]) {
    if init_hsps.len() <= 1 {
        return;
    }
    // Stable, as glibc's qsort under the pinned NCBI BLAST+: score_compare_match
    // ties ungapped alignments of the same segment from other seeds.
    init_hsps.sort_by(compare_init_hsps_by_score_ncbi);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3639-3657
fn compact_chained_init_hsps_ncbi(init_hsps: &mut Vec<InitHSP>, keep: &mut Vec<bool>) {
    while init_hsps
        .last()
        .is_some_and(|_| !keep[keep.len().saturating_sub(1)])
    {
        init_hsps.pop();
        keep.pop();
    }

    if init_hsps.len() <= 1 {
        return;
    }

    let mut index = 0usize;
    while index + 1 < init_hsps.len() {
        if keep[index] {
            index += 1;
            continue;
        }

        let end = init_hsps.len() - 1;
        debug_assert!(end > index);
        init_hsps[index] = init_hsps[end];
        keep[index] = keep[end];
        init_hsps.pop();
        keep.pop();

        while init_hsps.len() > index + 1 && init_hsps.last().is_some_and(|_| !keep[keep.len() - 1])
        {
            init_hsps.pop();
            keep.pop();
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3502-3658
fn chain_blastp_init_hsps(
    mut init_hsps: Vec<InitHSP>,
    contexts: &[QueryContext],
    word_cutoff_scores: &[i32],
    hit_cutoff_scores: &[i32],
    gap_open: i32,
    gap_extend: i32,
    keep: &mut Vec<bool>,
    nodes: &mut Vec<BlastpInitHspNode>,
) -> Vec<InitHSP> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3707-3708
    // A list of one ungapped alignment is chained too: its node keeps its own score
    // and the drop test below can remove it.
    if init_hsps.is_empty() {
        return init_hsps;
    }

    let gap_score = gap_open + gap_extend;
    ncbi_qsort_init_hsps_by_query_offset_score(&mut init_hsps);

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3524-3527
    keep.clear();
    keep.resize(init_hsps.len(), true);
    let mut i = 0usize;
    let mut context_idx = 0usize;
    while i < init_hsps.len() {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3536-3558
        while context_idx < contexts.len()
            && init_hsps[i].q_seed_absolute
                > contexts[context_idx].frame_base + contexts[context_idx].aa_len as i32
        {
            context_idx += 1;
        }
        assert!(
            context_idx < contexts.len(),
            "NCBI BLAST chaining requires InitHSP q_off to resolve to a query context"
        );
        let context_end = contexts[context_idx].frame_base + contexts[context_idx].aa_len as i32;
        let start = i;
        while i < init_hsps.len() && init_hsps[i].q_seed_absolute < context_end {
            i += 1;
        }
        let end = i;

        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3551-3578
        nodes.clear();
        nodes.extend((start..end).map(|init_index| BlastpInitHspNode {
            init_index,
            best_score: init_hsps[init_index].ungapped_data.score,
            next_index: None,
        }));

        for k in (0..nodes.len()).rev() {
            let self_score = nodes[k].best_score;
            for j in (k + 1)..nodes.len() {
                let left = init_hsps[nodes[k].init_index];
                let right = init_hsps[nodes[j].init_index];
                let left_len = left.ungapped_data.length;
                let q_diff = right.ungapped_data.q_start - left.ungapped_data.q_start + left_len;
                let s_diff = right.ungapped_data.s_start - left.ungapped_data.s_start + left_len;
                if s_diff < 0 {
                    continue;
                }

                let mut new_score = self_score + nodes[j].best_score;
                new_score += (q_diff.min(s_diff) * 3).min(word_cutoff_scores[context_idx]);
                new_score -= (q_diff - s_diff).abs().max(1);
                new_score -= gap_open;

                if new_score > nodes[k].best_score {
                    nodes[k].best_score = new_score;
                    nodes[k].next_index = Some(j);
                }
            }
        }

        for node in nodes.iter() {
            if node.best_score - gap_score + word_cutoff_scores[context_idx] - 1
                < hit_cutoff_scores[context_idx]
            {
                keep[node.init_index] = false;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3658
    compact_chained_init_hsps_ncbi(&mut init_hsps, keep);
    ncbi_qsort_init_hsps_by_score(&mut init_hsps);
    init_hsps
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3949-4010
fn should_retry_exact_preliminary_gapped_alignment(
    restricted_alignment: bool,
    raw_score: i32,
    cutoff: i32,
    restricted_cutoff: i32,
) -> bool {
    restricted_alignment && raw_score < cutoff && raw_score >= restricted_cutoff
}

fn preliminary_tree_hsp(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query_context_offset: i32,
) -> TreeHsp {
    TreeHsp {
        query_offset: preliminary_hsp.query_start,
        query_end: preliminary_hsp.query_end,
        subject_offset: preliminary_hsp.subject_start,
        subject_end: preliminary_hsp.subject_end,
        score: preliminary_hsp.raw_score,
        query_frame: 1,
        query_length: preliminary_hsp.query_length as i32,
        query_context_offset,
        subject_frame_sign: preliminary_hsp.subject_frame.signum(),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3919
#[inline]
fn blastp_trace_init_hsp_overlaps_target(target: BlastpTraceHspTarget, init_hsp: &InitHSP) -> bool {
    let q_start = match usize::try_from(init_hsp.ungapped_data.q_start + 1) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let q_end = match usize::try_from(init_hsp.query_end_absolute()) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let s_start = match usize::try_from(init_hsp.ungapped_data.s_start + 1) {
        Ok(value) => value,
        Err(_) => return false,
    };
    let s_end = match usize::try_from(init_hsp.subject_end()) {
        Ok(value) => value,
        Err(_) => return false,
    };

    q_start >= target.q_start
        && q_end <= target.q_end
        && s_start >= target.s_start
        && s_end <= target.s_end
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1330-1353
fn preliminary_score_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    b.raw_score
        .cmp(&a.raw_score)
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| b.query_end.cmp(&a.query_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2267-2315
fn preliminary_query_offset_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    a.query_context
        .cmp(&b.query_context)
        .then_with(|| a.query_start.cmp(&b.query_start))
        .then_with(|| a.subject_start.cmp(&b.subject_start))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.query_end.cmp(&a.query_end))
        .then_with(|| b.subject_end.cmp(&a.subject_end))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2332-2379
fn preliminary_query_end_compare(
    a: &BlastpPreliminaryHsp,
    b: &BlastpPreliminaryHsp,
) -> std::cmp::Ordering {
    a.query_context
        .cmp(&b.query_context)
        .then_with(|| a.query_end.cmp(&b.query_end))
        .then_with(|| a.subject_end.cmp(&b.subject_end))
        .then_with(|| b.raw_score.cmp(&a.raw_score))
        .then_with(|| b.query_start.cmp(&a.query_start))
        .then_with(|| b.subject_start.cmp(&a.subject_start))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2505
fn ncbi_qsort_preliminary_hits(
    hits: &mut [BlastpPreliminaryHsp],
    compar: fn(&BlastpPreliminaryHsp, &BlastpPreliminaryHsp) -> std::cmp::Ordering,
) {
    if hits.len() <= 1 {
        return;
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2495,2504-2520
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2320,2385
    hits.sort_by(compar);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:1355-1369
fn preliminary_hits_is_sorted_by_score_ncbi(hits: &[BlastpPreliminaryHsp]) -> bool {
    if hits.len() <= 1 {
        return true;
    }

    for index in 0..hits.len() - 1 {
        if preliminary_score_compare(&hits[index], &hits[index + 1]) == std::cmp::Ordering::Greater
        {
            return false;
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2478-2502
#[inline]
fn preliminary_hits_share_query_offset_endpoint(
    first: &BlastpPreliminaryHsp,
    candidate: &BlastpPreliminaryHsp,
) -> bool {
    first.query_context == candidate.query_context
        && first.query_start == candidate.query_start
        && first.subject_start == candidate.subject_start
        && first.subject_frame == candidate.subject_frame
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2504-2527
#[inline]
fn preliminary_hits_share_query_end_endpoint(
    first: &BlastpPreliminaryHsp,
    candidate: &BlastpPreliminaryHsp,
) -> bool {
    first.query_context == candidate.query_context
        && first.query_end == candidate.query_end
        && first.subject_end == candidate.subject_end
        && first.subject_frame == candidate.subject_frame
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2493-2500,2519-2525
fn purge_sorted_preliminary_endpoint_duplicates(
    hits: &mut Vec<BlastpPreliminaryHsp>,
    same_endpoint: fn(&BlastpPreliminaryHsp, &BlastpPreliminaryHsp) -> bool,
) -> usize {
    if hits.len() <= 1 {
        return 0;
    }

    let mut write_index = 1usize;
    for read_index in 1..hits.len() {
        if same_endpoint(&hits[write_index - 1], &hits[read_index]) {
            continue;
        }
        if write_index != read_index {
            hits[write_index] = hits[read_index];
        }
        write_index += 1;
    }
    let removed = hits.len().saturating_sub(write_index);
    hits.truncate(write_index);
    removed
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3870-3878
fn blastp_query_index_from_preliminary_context(
    contexts: &[QueryContext],
    query_context: i32,
) -> usize {
    let context_idx = usize::try_from(query_context)
        .expect("NCBI BLAST requires preliminary blastp query context to be non-negative");
    usize::try_from(contexts[context_idx].q_idx)
        .expect("NCBI BLAST requires blastp query index to fit in usize")
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2455-2538
fn purge_preliminary_hits_with_common_endpoints(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> (Vec<BlastpPreliminaryHsp>, usize) {
    if hits.len() <= 1 {
        return (hits, 0);
    }

    ncbi_qsort_preliminary_hits(&mut hits, preliminary_query_offset_compare);
    let mut purged = purge_sorted_preliminary_endpoint_duplicates(
        &mut hits,
        preliminary_hits_share_query_offset_endpoint,
    );

    ncbi_qsort_preliminary_hits(&mut hits, preliminary_query_end_compare);
    purged += purge_sorted_preliminary_endpoint_duplicates(
        &mut hits,
        preliminary_hits_share_query_end_endpoint,
    );

    (hits, purged)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:540-555
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3580-3581
#[cfg(test)]
fn prepare_preliminary_hits_for_kappa(
    mut hits: Vec<BlastpPreliminaryHsp>,
) -> Vec<BlastpPreliminaryHsp> {
    prepare_preliminary_hits_for_kappa_in_place(&mut hits);
    hits
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:540-555
fn prepare_preliminary_hits_for_kappa_in_place(hits: &mut Vec<BlastpPreliminaryHsp>) -> usize {
    let (purged, purged_count) = purge_preliminary_hits_with_common_endpoints(std::mem::take(hits));
    *hits = purged;
    if !preliminary_hits_is_sorted_by_score_ncbi(hits) {
        ncbi_qsort_preliminary_hits(hits, preliminary_score_compare);
    }
    purged_count
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:117-149
#[cfg(test)]
fn split_preliminary_hits_by_query(
    preliminary_hits: Vec<BlastpPreliminaryHsp>,
    contexts: &[QueryContext],
    num_queries: usize,
) -> Vec<Vec<BlastpPreliminaryHsp>> {
    let mut hits_by_query = vec![Vec::new(); num_queries];
    let mut preliminary_hits = preliminary_hits;
    split_preliminary_hits_by_query_into(
        &mut preliminary_hits,
        contexts,
        num_queries,
        &mut hits_by_query,
    );
    hits_by_query
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/hspfilter_collector.c:117-149
fn split_preliminary_hits_by_query_into(
    preliminary_hits: &mut Vec<BlastpPreliminaryHsp>,
    contexts: &[QueryContext],
    num_queries: usize,
    hits_by_query: &mut Vec<Vec<BlastpPreliminaryHsp>>,
) {
    if hits_by_query.len() < num_queries {
        hits_by_query.resize_with(num_queries, Vec::new);
    }
    for query_hits in hits_by_query.iter_mut() {
        query_hits.clear();
    }
    for preliminary_hsp in preliminary_hits.drain(..) {
        let query_index =
            blastp_query_index_from_preliminary_context(contexts, preliminary_hsp.query_context);
        if query_index < num_queries {
            hits_by_query[query_index].push(preliminary_hsp);
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-3997
fn remove_preliminary_hits_for_query(
    preliminary_hits: &mut Vec<BlastpPreliminaryHsp>,
    query_index: usize,
    contexts: &[QueryContext],
) -> usize {
    let before = preliminary_hits.len();
    preliminary_hits.retain(|hsp| {
        blastp_query_index_from_preliminary_context(contexts, hsp.query_context) != query_index
    });
    before.saturating_sub(preliminary_hits.len())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3577-3584
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:151-178
pub(super) fn preliminary_hit_from_local_match_hsp(local_hsp: &BlastpHsp) -> BlastpPreliminaryHsp {
    BlastpPreliminaryHsp {
        query_start: i32::try_from(
            local_hsp
                .q_start
                .checked_sub(1)
                .expect("NCBI BLAST localMatch query start must be 1-based"),
        )
        .expect("NCBI BLAST localMatch query start must fit in Int4"),
        query_end: i32::try_from(local_hsp.q_end)
            .expect("NCBI BLAST localMatch query end must fit in Int4"),
        subject_start: i32::try_from(
            local_hsp
                .s_start
                .checked_sub(1)
                .expect("NCBI BLAST localMatch subject start must be 1-based"),
        )
        .expect("NCBI BLAST localMatch subject start must fit in Int4"),
        subject_end: i32::try_from(local_hsp.s_end)
            .expect("NCBI BLAST localMatch subject end must fit in Int4"),
        gapped_query_start: local_hsp.gapped_q_start,
        gapped_subject_start: local_hsp.gapped_s_start,
        raw_score: local_hsp.raw_score,
        query_context: local_hsp.query_context,
        query_frame: local_hsp.query_frame,
        subject_frame: local_hsp.subject_frame,
        query_length: local_hsp.query_length,
        q_idx: local_hsp.q_idx,
        s_idx: local_hsp.s_idx,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3975-4000
fn rebuild_preliminary_interval_tree_from_list(
    interval_tree: &mut BlastIntervalTree,
    preliminary_hits: &[BlastpPreliminaryHsp],
    contexts: &[QueryContext],
) -> usize {
    interval_tree.reset();
    let mut inserted = 0usize;
    for preliminary_hsp in preliminary_hits {
        let context_idx = usize::try_from(preliminary_hsp.query_context)
            .expect("preliminary blastp query context must be non-negative");
        interval_tree.add_hsp(
            preliminary_tree_hsp(preliminary_hsp, contexts[context_idx].frame_base),
            contexts[context_idx].frame_base,
            IndexMethod::QueryAndSubject,
        );
        inserted += 1;
    }
    inserted
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3777-3786
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_traceback.c:635-688
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3949
// NCBI reference (598d8ae6): c++/src/algo/blast/core/aa_ungapped.c:1054,1083
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3394,3401,3419-3420,3435-3437
// The one-hit extension's best segment can end before the last reset of its running
// sum, so the Int4 length can be negative (a low -threshold with -window_size 0). NCBI
// passes it as Uint4: above HSP_MAX_WINDOW, and `q_start` plus `min(q_length, s_length)` wraps
// back to q_start + length, so only the first window of HSP_MAX_WINDOW letters is scored.
// A window that runs past the context or the subject holds its NULLB, whose
// BLAST_SCORE_MIN makes the sum negative whatever letters follow, so it gives q_start.
fn blastp_get_start_for_gapped_alignment_int4_length(
    query: &[u8],
    subject: &[u8],
    q_start: usize,
    length: i32,
    s_start: usize,
    matrix: ScoringMatrix,
) -> usize {
    if let Ok(length) = usize::try_from(length) {
        return blastp_get_start_for_gapped_alignment(
            query, subject, q_start, length, s_start, length, matrix,
        );
    }
    const HSP_MAX_WINDOW: usize = 11;
    let letter = |sequence: &[u8], index: usize| sequence.get(index).copied().unwrap_or(0);
    let score: i32 = (0..HSP_MAX_WINDOW)
        .map(|offset| {
            standard_score(
                matrix,
                letter(query, q_start + offset),
                letter(subject, s_start + offset),
            )
        })
        .sum();
    if score > 0 {
        q_start + HSP_MAX_WINDOW - 1 - HSP_MAX_WINDOW / 2
    } else {
        q_start
    }
}

fn extend_preliminary_blastp_hit(
    init_hsp: &mut InitHSP,
    s_idx: u32,
    context_idx: usize,
    ctx: &QueryContext,
    query_sequence: &[u8],
    subject_sequence: &[u8],
    matrix: ScoringMatrix,
    gap_open: i32,
    gap_extend: i32,
    x_dropoff: i32,
    mode: BlastpGappedAlignmentMode,
    scratch: &mut GapAlignScratch,
) -> BlastpPreliminaryHsp {
    let query_raw = &query_sequence[..ctx.aa_len];
    let subject_raw = &subject_sequence[..subject_sequence.len() - 1];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3905
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3900-3945
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3936-3945
    let gapped_q_start = blastp_get_start_for_gapped_alignment_int4_length(
        query_sequence,
        subject_sequence,
        usize::try_from(init_hsp.ungapped_data.q_start)
            .expect("NCBI BLAST preliminary q_start must be non-negative"),
        init_hsp.ungapped_data.length,
        usize::try_from(init_hsp.ungapped_data.s_start)
            .expect("NCBI BLAST preliminary s_start must be non-negative"),
        matrix,
    );
    let gapped_q_start = i32::try_from(gapped_q_start)
        .expect("NCBI BLAST preliminary gapped query start must fit in Int4");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3943-3945
    init_hsp.s_seed += gapped_q_start - init_hsp.q_seed_absolute;
    init_hsp.q_seed_absolute = gapped_q_start;
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3943-3949
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:3949-3978
    let aligned = blastp_score_only_gapped_alignment_with_scratch(
        query_raw,
        subject_raw,
        usize::try_from(gapped_q_start)
            .expect("NCBI BLAST preliminary gapped query start must be non-negative"),
        usize::try_from(init_hsp.s_seed)
            .expect("NCBI BLAST preliminary gapped subject start must be non-negative"),
        matrix,
        gap_open,
        gap_extend,
        x_dropoff,
        mode,
        scratch,
    );

    BlastpPreliminaryHsp {
        query_start: aligned.query_start,
        query_end: aligned.query_stop,
        subject_start: aligned.subject_start,
        subject_end: aligned.subject_stop,
        gapped_query_start: gapped_q_start,
        gapped_subject_start: init_hsp.s_seed,
        raw_score: aligned.score,
        query_context: i32::try_from(context_idx)
            .expect("NCBI BLAST preliminary blastp context index must fit in Int4"),
        query_frame: 0,
        subject_frame: 0,
        query_length: ctx.aa_len,
        q_idx: ctx.q_idx,
        s_idx,
    }
}

/// `settings` are the NCBI application settings that LOSAT reproduces
/// (`ncbi_environment::ApplicationSettings`): the readers of the query and the
/// subjects use their data loaders (`ReaderConfig`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_scope_src.cpp:67-72
pub fn run(
    mut args: BlastpArgs,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    args.id_scope = args.parse_deflines.then(Default::default);
    // The run's `SERIAL_VERIFY_DATA_WRITE`, which `-in_pssm`'s reader reads
    // (`BlastpArgs::write_verification`).
    args.write_verification = settings.write_verification.clone();
    use crate::blastinput::app;
    use crate::blastinput::input_files;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:128-133
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3624-3627
    // `ArchiveFormatRequested` parses `-outfmt` (blast_args.cpp:2745-2748) before the
    // option handlers run; a format that NCBI runs and LOSAT does not write is rejected
    // there too. `GetExportSearchStrategyStream` comes first: it creates the search strategy
    // file (`strategy::open_strategy_file`), which an error after it leaves empty; with it a
    // value that does not parse is left to the formatting handler
    // (`app::early_format_choice`).
    let strategy_file = crate::report::archive::strategy::open_strategy_file(
        args.export_search_strategy.as_deref(),
        args.out.as_deref(),
        args.app.name(),
    )?;
    let choice = app::early_format_choice(&args.outfmt, strategy_file.is_some())?;
    let format = choice.report_format(app::Program::Blastp, args.out.as_deref())?;
    // The BLAST archive (outfmt 11): its request, before anything is read or written.
    // The error archives follow the `-outfmt` text (`writes_error_archive`).
    // The run's options handle, which its archive and its search strategy share.
    let options = crate::report::archive::OptionsHandle::default();
    args.archive = super::archive::handle_for(
        &args,
        format == Some(crate::blastinput::app::ReportFormat::ArchiveFormat),
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
    // LOSAT's thread capability is checked before any input or output, as for BLASTN.
    crate::utils::threading::validate_threads(args.num_threads)?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2472-2476
    // With `-db` the handler of the database arguments keeps the name and reads the list: the
    // database is opened after every option handler has run (`search_cli`). With `-subject` it
    // reads the subjects here, before the query and the output are opened.
    let archive = args.archive.clone();
    let mut stderr = crate::cli::CliStderr::process();
    let mut archive_tee;
    let diagnostics: &mut (dyn Write + Send) = match &archive {
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
    let list = match &args.db {
        Some(_) => match args
            .db_filter
            .resolve(&settings, args.app.ncbi_name(), &mut *diagnostics)
        {
            Ok(list) => list,
            Err(error) => {
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        },
        None => None,
    };
    let subjects = match args.db {
        Some(_) => SubjectSet {
            records: std::borrow::Cow::Owned(Vec::new()),
            placements: seq_range::Placements::default(),
        },
        None => match read_cli_subjects(
            &args,
            &settings,
            &mut *diagnostics,
            archive.is_some() || args.strategy.is_some(),
        ) {
            Ok(subjects) => subjects,
            Err(error) => {
                // NCBI reads the `-subject` file inside `Run`, after the diagnostics handler
                // that saves the messages is set: an error archive follows the error.
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        },
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3456-3481
    // The query is opened and the output file created before the options are checked;
    // `-` is standard input and standard output.
    //
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:717-721
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:448-453
    // psiblast `-in_pssm` (`pssm_query`): the query is the PSSM, which the options' handler
    // reads (`check_options`). `-in_pssm` excludes `-query` (blast_args.cpp:1242-1244, the
    // parser's exclusion), so the query is the default `-` (`kDfltArgQuery`,
    // blast_args.cpp:3425-3427), which NCBI takes as `&cin` (no platform refuses it) and never
    // reads. LOSAT opens no query input then, as tblastn's `-in_pssm` (`tblastn_app`).
    let query_file = if pssm_query(&args) {
        None
    } else {
        input_files::check_utf8_file_name(&args.query, "query", args.app.name())?;
        // NCBI opens the query file inside `Run` too: its error is saved and an error archive
        // follows.
        match input_files::open_input(&args.query, "query", args.app.name()) {
            Ok(file) => Some(file),
            Err(error) => {
                write_run_error_archive(&args, archive.as_ref(), &error);
                return Err(error);
            }
        }
    };
    let out_file = open_out_file(&args)?;
    // Outfmt 0 to 4 (`blastp_pairwise_view_of`; psiblast refuses 1 to 4 before its search): the
    // stream of NCBI's `CBlastFormat`.
    let pairwise =
        format.is_some_and(|format| blastp_pairwise_view_of(args.app, format as i32).is_some());
    // The organism report's `-out` file is written as NCBI's buffer writes it.
    let tax_report = format == Some(crate::blastinput::app::ReportFormat::TaxFormat);
    // LOSAT's reports read the format number and the custom fields that NCBI keeps; a custom
    // tabular report of a delimiter, which `normalized` drops, reads the value as given, as the
    // library's: whether the app's formatter gets the delimiter (blastp's does, psiblast's and
    // rpsblast's do not, so their `6 delim=@` is the default outfmt 6) is the spec's, the one
    // that the reports read (`tabular::custom_tabular_spec`, `TabularSpec::of_app`).
    // A value not read yet (`EarlyFormat::Deferred`) fails at the options' handlers below.
    let outfmt = match choice.parsed() {
        Some(parsed)
            if parsed.delimiter.is_empty()
                || super::tabular::custom_tabular_spec(&args.outfmt, args.app).is_none() =>
        {
            parsed.normalized()
        }
        _ => args.outfmt.clone(),
    };
    let mut stream = crate::cli::ReportStream::of_cli(out_file, tax_report, args.strategy.as_ref());
    let strategy = args.strategy.clone();
    let result = {
        let mut outputs =
            ReportOutputs::single(&outfmt, OutputSink::Writer(&mut stream), &mut *diagnostics);
        search_cli(
            args,
            list,
            query_file,
            &subjects,
            &mut outputs,
            choice,
            settings,
        )
    };
    // A post of the search strategy left unwritten is LOSAT's internal error (`end_run`).
    let result = crate::report::archive::strategy::end_run(strategy.as_ref(), result);
    // NCBI writes an error archive after the archives of the batches before when the
    // application stops with one of its own errors, and when the run ends with messages that
    // no archive carried (`PrintErrorArchive` after its `CATCH_ALL`, whatever the status).
    if let Some(handle) = &archive {
        match &result {
            Err(error) => {
                // psiblast's `-in_pssm` file, which NCBI opens inside `Run` too.
                let in_run = crate::algorithm::pssm_search::in_pssm::file_error_in_run(error);
                if let Some(lines) = super::archive::error_lines(error, in_run) {
                    let _ = crate::report::archive::write_error_archive_after_post(
                        handle,
                        error,
                        &lines,
                        &mut stream,
                    );
                } else if let Some(lines) = crate::report::archive::silent_end_lines(handle, error)
                {
                    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:573-576
                    // A run that `CATCH_ALL` ended without a post (`cli::NativeError::silent_end`:
                    // psiblast's `-in_pssm` reader, NCBI's object manager): the archive of the
                    // messages saved before it, with no text before it, as tblastn's end writes.
                    let _ =
                        crate::report::archive::write_error_archive(handle, &lines, &mut stream);
                }
            }
            Ok(()) => {
                if handle.has_messages() {
                    let _ = crate::report::archive::write_error_archive(handle, &[], &mut stream);
                }
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:252-255
    // The stream of the outfmt 0 to 4 formatter throws when a write fails, as the organism
    // report's where NCBI's buffer writes it (`cli::ReportStream::finish`).
    stream.finish(result, pairwise || tax_report)
}

/// The error archive of a run that stopped in `Run` with `error` before its search wrote
/// anything, when NCBI writes one (`PrintErrorArchive`): into `-out`, or nothing when that
/// cannot be opened (NCBI's call swallows that).
fn write_run_error_archive(
    args: &BlastpArgs,
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
    let Some(lines) = super::archive::error_lines(error, true) else {
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
fn open_out_file(args: &BlastpArgs) -> Result<Option<std::io::BufWriter<std::fs::File>>> {
    args.out
        .as_deref()
        .filter(|path| path.as_os_str() != "-")
        .map(|path| {
            crate::cli::create_out_file(path, Some(args.app.name())).map(std::io::BufWriter::new)
        })
        .transpose()
}

/// The subjects of the command line's `-subject`, read where NCBI's handler of the
/// database arguments reads them (`run`). `diagnostics` takes the reader's messages;
/// `keep_deflines` is whether the BLAST archive or the search strategy is written (they need the
/// deflines).
fn read_cli_subjects(
    args: &BlastpArgs,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    diagnostics: &mut (dyn Write + Send),
    keep_deflines: bool,
) -> Result<SubjectSet<'static>> {
    use crate::blastinput::app;
    use crate::blastinput::fasta_reader::{read_subjects, FastaInputSource, ReaderConfig};
    use crate::blastinput::input_files;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2553-2562
    // The handler of the database arguments reads the subjects (an empty subject set
    // fails there), before the query and the output are opened.
    let Some(subject_path) = args.subject.as_deref() else {
        return Err(app::missing_subject_error());
    };
    input_files::check_utf8_file_name(subject_path, "subject", args.app.name())?;
    let subject_file = input_files::open_input(subject_path, "subject", args.app.name())?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2537-2545
    // The subject range is read after the subject file is opened, before it is read.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
    )?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_input_aux.cpp:242-246
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_fasta_input.cpp:328-330
    // The subjects of blastp are read as proteins (`IsProtein()`). NCBI reads the records
    // one at a time, writes each one's messages as it reads it, and checks each one's range
    // after reading it: a range that starts past the end of a record stops the reading
    // there (`read_subjects`). Records without residues (and intervals without letters) are
    // read without a message and reported when the formatter sets up the subjects
    // (`search`).
    let mut subject_source = FastaInputSource::from_file(
        subject_file,
        ReaderConfig::subject(args.app.name(), true, settings.data_loaders)
            .with_defline(keep_deflines)
            .with_parse_deflines(args.parse_deflines),
    );
    subject_source.set_id_scope(args.id_scope.clone());
    // Diagnostic (LOSAT_TIMING): the subject FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let mut read_subject_records = read_subjects(
        &mut subject_source,
        subject_range.as_ref(),
        &mut |message: &[u8]| std::io::Write::write_all(&mut *diagnostics, message),
    )?;
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(subject_source);
    if args.lcase_masking && subject_range.is_some() {
        x_soft_subject::x_keep_whole_letters(&mut read_subject_records);
    }
    let subjects = match seq_range::cut_subjects(&read_subject_records, subject_range.as_ref())? {
        Some((cut, placements)) => SubjectSet {
            records: std::borrow::Cow::Owned(cut),
            placements,
        },
        None => SubjectSet {
            records: std::borrow::Cow::Owned(read_subject_records),
            placements: seq_range::Placements::default(),
        },
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/objmgr_query_data.cpp:375-380
    // `Reset` of `m_Subjects` with a new `blast::CObjMgr_QueryFactory` of `subjects` (blast_args.cpp:2557)
    // fails for a subject input without records.
    if subjects.records.is_empty() {
        return Err(app::empty_subjects_error());
    }
    Ok(subjects)
}

/// The checks of the options alone, without the inputs (the `validate` of web ABI v2):
/// NCBI's processing of the options (`BlastpArgs::check_options`) and LOSAT's limits.
pub fn check_options(args: &BlastpArgs) -> Result<()> {
    // `-db` is checked as the command line checks it: the database is opened by the run
    // (`run_local_database`, which LOSAT Web's host gives the run's files; line WQA D12a).
    let resolved = args.resolve()?;
    validate_requested_blastp_support(&resolved)?;
    // Without the inputs: IDENTITY with the compressed lookup table, which `search` refuses where it
    // searches its first query batch.
    check_identity_compressed_lookup_of(&resolved)?;
    // The host validates the argv it runs, with its `-num_threads` (docs/web/abi_v2.md).
    crate::utils::threading::validate_threads(args.num_threads)
}

/// What `search_cli` and `run_local_database` do after the options are processed and before the
/// query: the opening of the database of `-db` (`list`: its taxid or seqid list, `None` for a
/// search without `-db`), and what the requested reports ask of it.
fn open_search_database(
    args: &BlastpArgs,
    resolved: &mut ResolvedBlastpArgs,
    list: Option<&crate::seqdb::DbList>,
    settings: &crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<Option<database::BlastpDatabase>> {
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:194-199
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:213-219
    // The database of `-db` is opened here: after the option handlers (the query file is
    // open and the output file exists) and before the test for an empty query.
    let database = match resolved.db.as_deref() {
        Some(name) => {
            let database = database::BlastpDatabase::open_restricted(
                name,
                settings,
                args.app.name(),
                list,
                resolved.db_mask.as_ref(),
                outputs.diagnostics,
            )?
            .showing_gis(args.show_gis);
            let rps = resolved.rps_search.is_some();
            let database = if rps {
                rps_hooks::with_units(resolved, database)?
            } else {
                database
            };
            // The requests of the archive and of the search strategy, which a `-in_pssm` run that
            // the object manager ends as its query joins the scope never writes
            // (`PsiSearch::joins_scope`; NCBI saves the strategy in `DoIterations`, after it).
            // rpsblast: each batch's `EffectiveSearchSpace` option is its last search unit's own
            // (`rps_hooks::request_space`, `run_call`'s `set_call_searched`), and its units read
            // `-dbsize` themselves (`rps_search::volumes::RpsUnits`): no program of RPS-BLAST takes
            // a list.
            if resolved
                .psi_search
                .as_deref()
                .is_none_or(crate::algorithm::pssm_search::PsiSearch::joins_scope)
            {
                crate::algorithm::common::database::note_archive_database_filter(
                    args.archive.as_ref(),
                    args.strategy.as_ref(),
                    &args.db_filter,
                    database.list_taxids(),
                    &database.statistics(),
                    false,
                    rps || args.searchsp.is_some_and(|space| space > 0),
                    !rps && args.dbsize.is_some_and(|size| size > 0),
                    args.app.name(),
                )?;
            }
            // rpsblast: what the application and each batch's RPS set-up set in the options
            // handle, which the requests of the archive and of the search strategy hold.
            if let (Some(rps), Some(handle)) = (
                resolved.rps_search.as_deref(),
                args.archive
                    .as_ref()
                    .map(crate::report::archive::ArchiveHandle::options)
                    .or_else(|| args.strategy.as_ref().map(|strategy| strategy.options())),
            ) {
                rps_hooks::note_handle_sets(handle, rps, resolved, &database)?;
            }
            crate::algorithm::common::database::check_formats_report_sums(
                || database.pairwise_database(),
                outputs
                    .formats
                    .iter_mut()
                    .map(|format| (format.outfmt, &mut format.refusal)),
                args.app.name(),
            )?;
            Some(database)
        }
        None => None,
    };
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:126-133
    // A custom tabular field list that looks every taxid up (the plural names of the taxonomy
    // database, the organism report) asks whether `taxdb.bti` is on the search path: the plural
    // names' formatters warn without it. A `-db` search asks its database (`BlastpDatabase::
    // taxdb_found`, at the search); a `-subject` search the search path of NCBI's databases, as
    // BLASTN's (`subject_taxonomy`).
    if database.is_none() {
        subject_taxonomy(
            args,
            resolved,
            Some(&crate::algorithm::common::database::TaxdbSetting::of(
                settings,
            )),
            outputs,
        )?;
    }
    // LOSAT's refusal of a database name for the Seq-annot formats, before the search
    // (`report::serial::seqalign::check_database`).
    if let Some(database) = database.as_ref() {
        for format in &outputs.formats {
            let choice = crate::blastinput::app::parse_formatting_string(format.outfmt)?;
            if let Some(serial) = crate::blastinput::app::ReportFormat::from_number(choice.number)
                .and_then(crate::report::serial::SerialFormat::of)
            {
                crate::report::serial::seqalign::check_database(
                    serial,
                    database.name(),
                    crate::blastinput::app::Program::Blastp,
                )?;
            }
        }
    }
    Ok(database)
}

/// What the reports of `outputs` of a `-subject` search ask of the taxonomy database, through
/// `setting` (`common::database::subject_taxonomy`): asked by the command line
/// (`open_search_database`) and by a search of records (`run_local_formats`) at the same point,
/// after the options are processed.
fn subject_taxonomy(
    args: &BlastpArgs,
    resolved: &mut ResolvedBlastpArgs,
    setting: Option<&crate::algorithm::common::database::TaxdbSetting>,
    outputs: &ReportOutputs<'_>,
) -> Result<()> {
    let taxonomy = crate::algorithm::common::database::subject_taxonomy(
        setting,
        args.app.name(),
        outputs
            .formats
            .iter()
            .filter_map(|format| super::tabular::custom_tabular_spec(format.outfmt, args.app)),
    )?;
    resolved.subject_taxdb = taxonomy.taxdb;
    resolved.subject_taxdb_found = taxonomy.found.unwrap_or(false);
    Ok(())
}

/// The part of `run` after the files are opened: NCBI's processing of the options,
/// `Query is Empty!`, the reading of the queries and the search.
fn search_cli(
    mut args: BlastpArgs,
    list: Option<crate::seqdb::DbList>,
    query_file: Option<crate::blastinput::input_files::InputFile>,
    subjects: &SubjectSet<'_>,
    outputs: &mut ReportOutputs<'_>,
    choice: crate::blastinput::app::EarlyFormat,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
) -> Result<()> {
    use crate::blastinput::fasta_reader::{
        read_queries, FastaInputSource, QueryEnd, QueryRecords, ReaderConfig,
    };
    // An ungapped search's `-matrix` that NCBI's tables do not have is looked for as a file with
    // the settings' search when the options are resolved (`BlastpArgs::matrix_search`); the first
    // query batch is refused then, and a run without a batch reads no file.
    args.matrix_search = Some(settings.matrix_search.clone());
    let mut resolved = args.check_options(outputs.diagnostics)?;
    // rpsblast: NCBI's check of the database's frequency ratios comes between the options and
    // the database (`rps_hooks::check_freq_ratio_file`).
    if resolved.rps_search.is_some() {
        rps_hooks::check_freq_ratio_file(&mut resolved, &settings, outputs.diagnostics)?;
    }
    let database = open_search_database(&args, &mut resolved, list.as_ref(), &settings, outputs)?;
    // rpsblast `-num_threads 0`: one thread for each volume of the database.
    if let (Some(rps), Some(database)) = (resolved.rps_search.as_deref(), database.as_ref()) {
        if rps.auto_threads {
            resolved.num_threads = rps_hooks::volume_count(database);
        }
    }
    // psiblast sets its subjects up here, before the query (`psiblast_check_subject_letters`).
    if resolved.app == BlastpApp::Psiblast && database.is_none() {
        psiblast_check_subject_letters(&subjects.records)?;
    }
    // psiblast `-in_pssm`: the query is the PSSM's (read with the options, `check_options`);
    // NCBI reads no query input (c++/src/app/blast/psiblast_app.cpp:452-471). LOSAT's limits of the PSSM search
    // come when the options are processed and the subjects set up (`search_pssm_query`).
    if let Some(psi) = resolved.psi_search.clone() {
        let choice_to_standard_output = resolved
            .out
            .as_deref()
            .is_none_or(|path| path.as_os_str() == "-");
        return search_pssm_query(
            resolved,
            &psi,
            args.id_scope.as_deref(),
            subjects,
            database.as_ref(),
            "",
            outputs,
            choice.parsed().filter(|_| choice_to_standard_output),
        );
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:202-216
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:856-860
    // The query of blastp is a protein (`QueryIsProtein()`, `fAssumeProt`) and its data
    // loaders are those of the subjects (`SDataLoaderConfig`). The position is taken on the
    // opened file, before it is read (`stream_is_empty`). When the subjects were read from
    // standard input too, `cin` has reached its end (a failed stream), so NCBI gets no
    // position for the query either, and its reader is at its end: no query batch is read.
    // The query input is opened for every query but the PSSM of `-in_pssm` (`pssm_query`), which
    // the branch above searched (`psi_search` is that PSSM, read with the options).
    let Some(query_file) = query_file else {
        return Err(anyhow::anyhow!(
            "LOSAT internal error: a FASTA query search without its query input"
        ));
    };
    let shared_standard_input = args.query.as_os_str() == "-"
        && args
            .subject
            .as_deref()
            .is_some_and(|path| path.as_os_str() == "-");
    let mut query_source = FastaInputSource::from_file(
        query_file,
        ReaderConfig::query(args.app.name(), true, settings.data_loaders)
            .with_defline(args.archive.is_some() || args.strategy.is_some())
            .with_parse_deflines(args.parse_deflines),
    );
    query_source.set_id_scope(args.id_scope.clone());
    if !shared_standard_input && query_source.stream_is_empty() {
        // psiblast_app.cpp:461-464 is the same test, with psiblast's prefix.
        writeln!(
            outputs.diagnostics,
            "Warning: [{}] Query is Empty!",
            args.app.ncbi_name()
        )?;
        // NCBI returns from `Run` here, before its `PrintErrorArchive`: the message is not an
        // archive's.
        if let Some(handle) = &args.archive {
            handle.take_messages();
        }
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:256-259
    // The records are read here, and their messages are written batch by batch
    // (`read_query_batches`). What ends the reading (the end of the input, an `eEOF`, a
    // reader error or LOSAT's rejection of a Seq-id line) comes in the batch that NCBI reads
    // it in, after the reports of the batches before.
    // Diagnostic (LOSAT_TIMING): the query FASTA for the common line (utils/xtiming.rs). No NCBI
    // counterpart.
    let t_read_fasta = crate::utils::xtiming::clock();
    let QueryRecords { records, end } = read_queries(&mut query_source);
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::ReadFasta, t_read_fasta);
    drop(query_source);
    // The XInclude formats fail only on standard output; to a file they take its name as the
    // base of their side files.
    let choice_to_standard_output = resolved
        .out
        .as_deref()
        .is_none_or(|path| path.as_os_str() == "-");
    let xinclude = choice.parsed().filter(|_| choice_to_standard_output);
    // The refusal of a PSSM file that LOSAT cannot write as NCBI does (one of the report
    // streams, a named pipe) comes before the first byte (`check_psiblast_pssm_files`). (A query
    // title with a byte that NCBI's stream replaces is written as NCBI writes it:
    // `pssm::file_posts`.)
    check_psiblast_pssm_files(&resolved, &records, &pssm_check_scope(&args), xinclude)?;
    let reading = match end {
        QueryEnd::Input => QueryReading::End,
        QueryEnd::EofError => QueryReading::BlankLines,
        QueryEnd::Error { error, warnings } => QueryReading::Error { error, warnings },
    };
    search(
        resolved,
        &records,
        reading,
        subjects,
        database.as_ref(),
        "",
        outputs,
        xinclude,
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:195-295
// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2891-2978
/// Runs one BLASTP search over records read as NCBI's reader reads them
/// (`blastinput::fasta_reader`; their messages are written where NCBI reads them) and
/// writes every requested output format from the same result (the shared entry of the
/// CLI, web ABI v1 and v2). Queries without a record are NCBI's `Query is Empty!`.
///
/// Each requested `-outfmt` is parsed as on the command line, before the search starts.
/// The search options are resolved once per format and must not differ between formats:
/// NCBI derives the hitlist size from the output format only through `-num_descriptions`
/// and `-num_alignments`, which BLASTP does not accept, so a difference is an error here.
///
/// The BLAST archive (outfmt 11) among the formats writes each record's DefLine as the reader
/// kept it (`ReaderConfig::with_defline(true)`). Records without a kept defline get `>` and
/// their title, which differs from NCBI for a defline with extra spaces or control bytes.
pub fn run_local(
    mut args: BlastpArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[FastaRecord],
    query_label: &str,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTP_TASK.check(&mut args.task)?;
    // `run_local` searches subject records (the library, `api::run`; LOSAT Web's adapter; ABI v1);
    // a `-db` search runs through `run_local_database` (or the command line's `run`).
    if args.db.is_some() {
        bail!(crate::api::local_blast::RUN_LOCAL_DB_REFUSAL);
    }
    report_formats(&args, outputs)?;
    // The BLAST archive (outfmt 11) among the formats: its request and a diagnostics stream that
    // saves the messages.
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| blastp_is_archive(format.outfmt));
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    if archive_requested {
        args.archive = super::archive::handle_for(
            &args,
            true,
            &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
        )
        .map(|handle| outputs.with_archive_outfmt(handle));
    }
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input() && !pssm_query(&args);
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        run_local_formats(
            args,
            query_records,
            query_end,
            subject_records,
            query_label,
            subject_label,
            outputs,
        )
    })
}

// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:534-543
/// psiblast `-in_pssm`: the search of the PSSM's query (the file read with the options,
/// `check_options`), from where the entries would read their FASTA queries: the command line's
/// (`search_cli`) and those of the records (`run_local_formats`, `search_database`, the library
/// API's `Inputs::pssm`). A request whose PSSM file has members that the reader skipped, and the
/// input PSSM's files that NCBI would write back in a way LOSAT does not (`-out_pssm`,
/// `-out_ascii_pssm`, written in round 1, which a run stopped before its search never reaches),
/// are refused here, once the database is opened (NCBI's error of a database it cannot open comes
/// first) and before any output; so are LOSAT's limits of the PSSM search for each requested
/// format (`check_in_pssm_options`) and the PSSM files of the rounds that LOSAT cannot write as NCBI does
/// (`check_pssm_files`, checked as a FASTA query's).
///
/// The query is the PSSM's Bioseq as the reports show it (`PsiSearch::query_record`); with
/// `-parse_deflines` its ids join `scope`, the scope of the `-subject` records (`IdScope::
/// add_record`: NCBI shows a subject with one of them as the query, which LOSAT refuses as for a
/// FASTA query). Its messages are those of NCBI's object manager about the query's annotations
/// (`PsiSearch::annot`), written after the prolog as a FASTA query's batch writes its reader's.
/// The library's observer receives it (`FormatObserver::pssm_query`). One search, no query
/// batches.
#[allow(clippy::too_many_arguments)]
fn search_pssm_query(
    resolved: ResolvedBlastpArgs,
    psi: &std::sync::Arc<crate::algorithm::pssm_search::PsiSearch>,
    scope: Option<&crate::seqid::id_scope::IdScope>,
    subjects: &SubjectSet<'_>,
    database: Option<&database::BlastpDatabase>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    xinclude: Option<&crate::blastinput::app::FormatChoice>,
) -> Result<()> {
    let program = resolved.app.name();
    // The PSSM files of round 1 (written back from the input PSSM, or the engine's), which a run
    // that never searches does not write (NCBI BLAST+ 2.17.0 leaves `-out_pssm` and
    // `-out_ascii_pssm` alone when the object manager ends the run: line RG D3-1 fix 3 S1, fix 4).
    // The query set-up's stop comes after round 1's files (`PssmQueryStop::after_round_files`),
    // so their refusals apply to it too.
    let searches = psi.searches(resolved.strategy.is_some(), database.is_some());
    let writes_round_files = psi.joins_scope()
        && psi
            .stop_before_search(resolved.strategy.is_some(), database.is_some())
            .is_none_or(|stop| stop.after_round_files());
    if writes_round_files {
        if let Some(input) = resolved
            .psi_rounds
            .as_deref()
            .and_then(|plan| plan.input.as_ref())
        {
            input.refuse_written_back()?;
        }
    }
    let outfmts = outputs
        .formats
        .iter()
        .map(|format| {
            crate::blastinput::app::parse_formatting_string(format.outfmt)
                .map(|choice| choice.number)
        })
        .collect::<Result<Vec<i32>>>()?;
    // An XInclude format without a file name ends the run with NCBI's error where the formatter
    // is made (below, psiblast_app.cpp:503), before the prolog and the search: LOSAT's refusals of
    // the search's reports wait for it (`xinclude_error_comes`; NCBI 2.17.0, line RG F/d36 probe
    // ps.pssm.o14: `Please provide a file name for outfmt 14.`).
    if !xinclude_error_comes(xinclude) {
        crate::algorithm::psiblast::in_pssm::check_in_pssm_options(&resolved, psi, &outfmts)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:482-484 (the formatter's
    // `believe_query` is `-parse_deflines`, with `-in_pssm` too)
    // With `-parse_deflines` the reports believe the PSSM query's local Seq-id as they believe
    // a FASTA query's (NCBI BLAST+ 2.17.0: `Query= Query_1 <title>`, `qseqid` `Query_1`); its
    // other ids they show without it (`PsiSearch::query_record`). The ids the reports show
    // join the search's scope as a FASTA query's do: NCBI shows a subject or a database
    // sequence with one of them as the PSSM's query (its title and letters), refused as for a
    // FASTA query (`IdScope::add_record` of the `-subject` records that `-parse_deflines`
    // reads; the database's, `query_batch::QueryScopes`).
    let mut record = psi.query_record_in(resolved.parse_deflines, scope, program, searches)?;
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:506, 537 (cited at
    // `PsiSearch::write_annot_posts`)
    // The messages of NCBI's object manager about the query's annotations come after the
    // prolog, where a FASTA query's batch has its reader's messages: they are the PSSM query's
    // (`psiblast_stop_at_pssm_query` marks the stream they fail).
    record.warnings = psi.annot.posts.text();
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux.cpp:1155-1172
    // With the archive as the first format (no prolog before them) each of these posts is saved
    // as one message of its severity, as tblastn's (`PsiSearch::write_annot_posts`), where the
    // diagnostics stream would save each of its lines as a warning: they are written here.
    if let (Some(11), Some(handle)) = (outfmts.first(), resolved.archive.as_ref()) {
        psi.annot
            .posts
            .write(&mut *outputs.diagnostics, Some(handle))?;
        record.warnings.clear();
    }
    let records = vec![record];
    // The PSSM files of `-in_pssm`'s rounds (item PS5b), checked as a FASTA query's
    // (`check_psiblast_pssm_files`: not where NCBI's error of the format comes first).
    if let Some(request) = resolved
        .psi_rounds
        .as_deref()
        .and_then(crate::algorithm::psiblast::rounds::PsiRoundPlan::file_request)
        .filter(|_| searches && !xinclude_error_comes(xinclude))
    {
        // The strategy's request before the PSSM file (`check_psiblast_pssm_files`).
        crate::report::archive::strategy::check_report_stream(resolved.strategy.as_ref())?;
        crate::algorithm::psiblast::pssm::check_pssm_files(
            request,
            records.len(),
            resolved
                .out
                .as_deref()
                .filter(|path| path.as_os_str() != "-"),
        )?;
    }
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.pssm_query(&records[0]);
    }
    search(
        resolved,
        &records,
        QueryReading::End,
        subjects,
        database,
        subject_label,
        outputs,
        xinclude,
    )
}

/// Whether the query of `args` is a PSSM (psiblast's `-in_pssm` or `-in_msa`): a run of the records' entries
/// with no query record then has a query (`ReportOutputs::finish_archive`).
fn pssm_query(args: &BlastpArgs) -> bool {
    args.psiblast_pssm
        .as_ref()
        .is_some_and(|options| options.in_pssm.is_some() || options.in_msa.is_some())
}

/// Before the search of the FASTA queries `records`: the refusal of a PSSM file that LOSAT cannot
/// write as NCBI does (one of the report streams, a named pipe) comes before the first byte, for the
/// PSSM after the last round (`psi_pssm`) or those of the rounds (`PsiRoundPlan::file_request`). (A
/// query title with a byte that NCBI's stream replaces is written as NCBI writes it,
/// `pssm::file_posts`.) The command line's entry (`search_cli`) makes it where it has read the queries, the
/// records' entries (`run_local_formats`, `search_database`) where it would, for records that the
/// command line reads from a file that is not empty (no records: the file of `Query is Empty!`,
/// before the queries are read). `scope` says which queries the check covers: the call's
/// `records`, or the run's in the first of the library's query batches (none in the later ones).
/// Where the run ends with NCBI's error of an XInclude format without a file name (`xinclude`,
/// `xinclude_error_comes`) NCBI writes no PSSM file (its formatter is made before `DoIterations`,
/// psiblast_app.cpp:503 before :517), so nothing is refused.
fn check_psiblast_pssm_files(
    resolved: &ResolvedBlastpArgs,
    records: &[FastaRecord],
    scope: &crate::algorithm::psiblast::pssm::PssmCheckScope,
    xinclude: Option<&crate::blastinput::app::FormatChoice>,
) -> Result<()> {
    use crate::algorithm::psiblast::pssm::PssmCheckScope;
    if xinclude_error_comes(xinclude) {
        return Ok(());
    }
    let queries = match scope {
        PssmCheckScope::Call => records.len(),
        PssmCheckScope::Run { queries } => *queries,
        PssmCheckScope::Checked => return Ok(()),
    };
    if let Some(request) = resolved.psi_pssm.as_ref().or_else(|| {
        resolved
            .psi_rounds
            .as_deref()
            .and_then(crate::algorithm::psiblast::rounds::PsiRoundPlan::file_request)
    }) {
        // NCBI writes the strategy's request (`SaveSearchStrategy`, psiblast_app.cpp:340) before
        // the PSSM file (`SavePssmToFile`, :385): LOSAT's refusal of the strategy file comes first.
        crate::report::archive::strategy::check_report_stream(resolved.strategy.as_ref())?;
        crate::algorithm::psiblast::pssm::check_pssm_files(
            request,
            queries,
            resolved
                .out
                .as_deref()
                .filter(|path| path.as_os_str() != "-"),
        )?;
    }
    Ok(())
}

/// Whether the run ends with NCBI's error of the XInclude format `xinclude` (13 or 14 without a file
/// name, which `search` gives where NCBI makes its formatter, `xinclude_check`): LOSAT's refusals of
/// what NCBI writes only after that point (psiblast's PSSM files, `-in_pssm`'s options) wait for
/// it. The command line hands over its choice whenever it writes to standard output, the records'
/// entries their first format without a file (`xinclude_format`).
fn xinclude_error_comes(xinclude: Option<&crate::blastinput::app::FormatChoice>) -> bool {
    xinclude.is_some_and(|choice| crate::blastinput::app::xinclude_check(choice).is_err())
}

/// The queries that the call's checks of its PSSM files cover (`PsiblastArgs::pssm_check`).
fn pssm_check_scope(args: &BlastpArgs) -> crate::algorithm::psiblast::pssm::PssmCheckScope {
    args.psiblast_pssm
        .as_ref()
        .map(|options| options.check_scope.clone())
        .unwrap_or_default()
}

fn run_local_formats(
    args: BlastpArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    subject_records: &[FastaRecord],
    _query_label: &str,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The thread count is checked before the inputs are read, as on the command line (`run`).
    crate::utils::threading::validate_threads(args.num_threads)?;
    // The subject range and its record checks come where NCBI reads the subjects (`run`):
    // the messages of each record read, up to a record whose range starts past its end.
    let subject_range = seq_range::parse_optional_range(
        args.subject_loc.as_deref(),
        seq_range::RangeRole::Subject,
    )?;
    for record in
        &subject_records[..seq_range::subjects_read(subject_records, subject_range.as_ref())]
    {
        outputs.diagnostics.write_all(&record.warnings)?;
    }
    let x_whole = (args.lcase_masking && subject_range.is_some()).then(|| {
        let mut records = subject_records.to_vec();
        x_soft_subject::x_keep_whole_letters(&mut records);
        records
    });
    let subjects = match seq_range::cut_subjects(
        x_whole.as_deref().unwrap_or(subject_records),
        subject_range.as_ref(),
    )? {
        Some((cut, placements)) => SubjectSet {
            records: std::borrow::Cow::Owned(cut),
            placements,
        },
        None => SubjectSet {
            records: std::borrow::Cow::Borrowed(subject_records),
            placements: seq_range::Placements::default(),
        },
    };
    if subject_records.is_empty() {
        return Err(crate::blastinput::app::empty_subjects_error());
    }
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    let mut resolved = resolve_for_formats(&args, outputs)?;
    // The taxonomy database of the reports, where the command line asks (`open_search_database`).
    subject_taxonomy(&args, &mut resolved, args.taxdb_setting.as_ref(), outputs)?;
    let xinclude = xinclude_format(&args, outputs)?;
    let scope = pssm_check_scope(&args);
    let args = resolved;
    if args.app == BlastpApp::Psiblast {
        psiblast_check_subject_letters(&subjects.records)?;
    }
    // psiblast `-in_pssm` (the library API's `Inputs::pssm`): the query is the PSSM's, whose ids
    // join the scope of the subject records read (`PsiSearch::scope_of_subjects`).
    if let Some(psi) = args.psi_search.clone() {
        let id_scope = crate::algorithm::pssm_search::scope_of_subjects(
            args.parse_deflines,
            subject_records,
            args.app.name(),
        )?;
        return search_pssm_query(
            args,
            &psi,
            id_scope.as_ref(),
            &subjects,
            None,
            subject_label,
            outputs,
            xinclude.as_ref(),
        );
    }
    if !query_records.is_empty() {
        check_psiblast_pssm_files(&args, query_records, &scope, xinclude.as_ref())?;
    }
    search(
        args,
        query_records,
        query_end.reading(),
        &subjects,
        None,
        subject_label,
        outputs,
        xinclude.as_ref(),
    )
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:194-199
/// Runs one BLASTP search of the records `query_records` against the database of `args` and
/// writes every requested output format from the same result (the library API's entry; the
/// command line's is `run`). `settings` are the NCBI application settings, which find the
/// database. The steps are `run`'s and `search_cli`'s, in their order; what differs is what they
/// read: the records are read already, and the outputs are the caller's. An archive (outfmt 11)
/// ends as the command line's does (`ReportOutputs::finish_archive`).
pub fn run_local_database(
    mut args: BlastpArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // NCBI's options factory checks the task first (`TaskParser::check`); the command line's
    // parser has checked it, a caller of `run_local` has not.
    super::args::BLASTP_TASK.check(&mut args.task)?;
    report_formats(&args, outputs)?;
    args.write_verification = settings.write_verification.clone();
    let archive_requested = outputs
        .formats
        .iter()
        .any(|format| blastp_is_archive(format.outfmt));
    // The archive shares the options handle of the run's search strategy (LOSAT Web's
    // `-export_search_strategy`, `web_adapter::give_strategy`), as on the command line (`run`).
    args.archive = super::archive::handle_for(
        &args,
        archive_requested,
        &crate::report::archive::strategy::options_of(args.strategy.as_ref()),
    )
    .map(|handle| outputs.with_archive_outfmt(handle));
    crate::utils::threading::validate_threads(args.num_threads)?;
    let archive = args.archive.clone();
    let strategy = args.strategy.clone();
    let no_query = query_records.is_empty() && query_end.is_input() && !pssm_query(&args);
    outputs.run_records(strategy.as_ref(), archive.as_ref(), no_query, |outputs| {
        search_database(args, query_records, query_end, settings, outputs)
    })
}

/// `run_local_database` from the taxid or seqid list on.
fn search_database(
    mut args: BlastpArgs,
    query_records: &[FastaRecord],
    query_end: &QueryEnd,
    settings: crate::blastinput::ncbi_environment::ApplicationSettings,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    // The handler of the database arguments reads the list before the other handlers (`run`).
    let list =
        args.db_filter
            .resolve(&settings, args.app.ncbi_name(), &mut *outputs.diagnostics)?;
    // `-out` (`run_to_file`) is made where the command line opens it (`run`).
    outputs.open_out()?;
    args.matrix_search = Some(settings.matrix_search.clone());
    let mut resolved = resolve_for_formats(&args, outputs)?;
    let database = open_search_database(&args, &mut resolved, list.as_ref(), &settings, outputs)?
        .ok_or_else(crate::blastinput::app::missing_subject_error)?;
    let xinclude = xinclude_format(&args, outputs)?;
    let no_subjects = SubjectSet {
        records: std::borrow::Cow::Borrowed(&[]),
        placements: seq_range::Placements::default(),
    };
    // psiblast `-in_pssm` (the library API's `Inputs::pssm_database`): the query is the PSSM's. No
    // subject record holds an id (the database's sequences are checked against the query's,
    // `query_batch::QueryScopes`).
    if let Some(psi) = resolved.psi_search.clone() {
        return search_pssm_query(
            resolved,
            &psi,
            None,
            &no_subjects,
            Some(&database),
            "",
            outputs,
            xinclude.as_ref(),
        );
    }
    if !query_records.is_empty() {
        check_psiblast_pssm_files(
            &resolved,
            query_records,
            &pssm_check_scope(&args),
            xinclude.as_ref(),
        )?;
    }
    search(
        resolved,
        query_records,
        query_end.reading(),
        &no_subjects,
        Some(&database),
        "",
        outputs,
        xinclude.as_ref(),
    )
}

/// The first requested format that is an XInclude format (13 or 14) without a file
/// (`ReportOutputs::xinclude_choice`), which `search` refuses where the command line refuses its
/// `-outfmt` (after `Query is Empty!` and the subjects without letters).
fn xinclude_format(
    args: &BlastpArgs,
    outputs: &ReportOutputs<'_>,
) -> Result<Option<crate::blastinput::app::FormatChoice>> {
    outputs.xinclude_choice(args.out.as_deref())
}

/// LOSAT's rejection of a protein subject title of punctuation that NCBI's
/// `x_CleanAndCompress` reads past the end of in the outfmt 0 titles of the subjects that
/// the report shows (NCBI BLAST+ 2.17.0 crashes there). The approved exception 2 of
/// PD-LOSAT-NCBI-DEFECTS covers the nucleotide subjects of BLASTN, TBLASTX and TBLASTN
/// only (DW-17), so LOSAT's BLASTP keeps the rejection.
///
/// NCBI reference (598d8ae6): c++/src/objmgr/util/create_defline.cpp:4092-4095
fn check_shown_subject_titles(
    shown: impl IntoIterator<Item = usize>,
    subject_records: &[FastaRecord],
) -> Result<()> {
    for index in shown {
        let title = &subject_records[index].title;
        let reads_past_end = [false, true].into_iter().any(|leave_prefix| {
            crate::report::defline::generate_defline(title, true, leave_prefix)
                .is_ok_and(|defline| defline.reads_past_end)
        });
        if reads_past_end {
            bail!(
                "subject record {} has a title of punctuation that NCBI BLAST+'s x_CleanAndCompress reads past the end of (NCBI BLAST+ 2.17.0 crashes in outfmt 0); this is not supported by LOSAT's BLASTP",
                index + 1
            );
        }
    }
    Ok(())
}

/// The titles that a pairwise report shows, checked where the subjects are a database's
/// (`ListedTitles::check_shown_titles`, which takes the subjects whose load throws from `loads`)
/// or a FASTA set's (`check_shown_subject_titles`).
fn check_pairwise_shown_titles(
    listed_titles: Option<&database::ListedTitles<'_>>,
    database: Option<&database::BlastpDatabase>,
    shown: &std::collections::BTreeSet<usize>,
    subject_records: &[FastaRecord],
    loads: &crate::report::db_subject::HitLoads,
) -> Result<()> {
    match listed_titles {
        Some(titles) => titles.check_shown_titles(
            database.expect("listed titles are a database's"),
            shown,
            loads,
        ),
        None => check_shown_subject_titles(shown.iter().copied(), subject_records),
    }
}

/// The subjects of a search as it searches them (each cut to its `-subject_loc` interval)
/// and where each lies in its record.
struct SubjectSet<'a> {
    records: std::borrow::Cow<'a, [FastaRecord]>,
    placements: seq_range::Placements,
}

/// What a search reads its subjects from: the records of a FASTA subject set (`-subject`;
/// every record is encoded and kept for the whole search), or a database (`-db`; read in
/// blocks, and only the sequences of the hit lists are kept, `database`).
#[derive(Clone, Copy)]
enum BlastpSubjectSource<'a> {
    Fasta(&'a [FastaRecord]),
    Database(&'a database::BlastpDatabase),
}

impl<'a> BlastpSubjectSource<'a> {
    /// The records of a FASTA subject set; none for a database.
    fn fasta_records(self) -> &'a [FastaRecord] {
        match self {
            Self::Fasta(records) => records,
            Self::Database(_) => &[],
        }
    }

    fn database(self) -> Option<&'a database::BlastpDatabase> {
        match self {
            Self::Fasta(_) => None,
            Self::Database(database) => Some(database),
        }
    }

    /// The number of sequences and of letters that the `Database:` lines of the reports
    /// show: the real totals, also with `-dbsize`, as NCBI's `Int8` (`PairwiseDatabase`).
    fn shown_totals(self) -> (i64, i64) {
        match self {
            Self::Fasta(records) => (
                records.len() as i64,
                records.iter().map(|record| record.seq().len() as i64).sum(),
            ),
            Self::Database(database) => {
                let shown = database.pairwise_database();
                (shown.num_sequences, shown.total_letters)
            }
        }
    }
}

/// One subject as the search kernel reads it: the residues in NCBIstdaa with the
/// sentinel byte before and after them (`EncodedProtein::aa_seq`, or the bytes of a
/// database block), and their number.
#[derive(Clone, Copy)]
struct BlastpSubjectRef<'a> {
    aa_seq: &'a [u8],
    aa_len: usize,
}

impl<'a> BlastpSubjectRef<'a> {
    fn of(subject: &'a EncodedProtein) -> Self {
        Self {
            aa_seq: &subject.aa_seq,
            aa_len: subject.aa_len,
        }
    }
}

/// The ranges of a search as the search and the report use them: the query input (every
/// record read, and where each searched query lies in its record) and where each
/// subject's searched letters lie in its record.
struct BlastpRanges<'a> {
    input_records: &'a [FastaRecord],
    query_input: &'a seq_range::QueryInput,
    subject_placements: &'a seq_range::Placements,
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the input position of the first query of a block (0 for a
    // search of every query), added to `query_input`'s positions where a query is named; and
    // psiblast's rounds of one query (`psi_round::run_queries`): its input position.
    x_ordinal_base: usize,
    /// The index among the run's searched queries of the call's first query: where the search
    /// strategy's requests of the call's batches start (`CallRequests`). It is the input position
    /// of `x_ordinal_base` unless `-query_loc` skips records before the call's queries (a range
    /// that starts past the end of its record: psiblast's rounds of a query after such a record).
    searched_base: usize,
    /// A query chunk's search (`search_split_query_batches`): the options' search space of each
    /// context of the chunk's whole batch, which NCBI sets as the chunk's options
    /// (`SplitQuery_SetEffectiveSearchSpace`: 0 for a context that is not valid, unless
    /// `-searchsp` gives it) and the chunk reads by its own context index
    /// (`split_option_space`); `None` for every other search.
    split_search_spaces: Option<&'a [i64]>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:676-697
/// The option space that a query chunk's context `context_index` reads of its batch's
/// (`BlastpRanges::split_search_spaces`, one per batch query): the batch's space of the same
/// index, the one space of a batch of one query, or the user's `-searchsp` (`options_searchsp`,
/// which `SplitQuery_SetEffectiveSearchSpace` leaves in the options).
fn split_option_space(batch: &[i64], context_index: usize, options_searchsp: i64) -> i64 {
    if options_searchsp != 0 {
        return options_searchsp;
    }
    match batch {
        [one] => *one,
        _ => batch[context_index],
    }
}

/// Each requested format as NCBI parses `-outfmt` before its option handlers run, and LOSAT's
/// rejection of a format that its BLASTP does not write (`ReportOutputs::early_choices`): the
/// command line's checks before it opens `-out` (`run`). The formatting handler's come after
/// (`resolve_for_formats`), and so does the error of a value left to it with a search strategy
/// file.
fn report_formats(args: &BlastpArgs, outputs: &ReportOutputs<'_>) -> Result<()> {
    outputs
        .early_choices(
            crate::blastinput::app::Program::Blastp,
            args.out.as_deref(),
            args.export_search_strategy.is_some(),
        )
        .map(drop)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:2894-2978
fn resolve_for_formats(
    args: &BlastpArgs,
    outputs: &mut ReportOutputs<'_>,
) -> Result<ResolvedBlastpArgs> {
    // Each format is checked as the CLI checks its single `-outfmt`; NCBI's warnings of the
    // option handlers are written once, for the first format.
    let mut resolved_formats = Vec::with_capacity(outputs.formats.len());
    for (index, format) in outputs.formats.iter().enumerate() {
        let mut requested = args.clone();
        requested.outfmt = format.outfmt.to_string();
        let choice = crate::blastinput::app::early_format_choice(
            &requested.outfmt,
            args.export_search_strategy.is_some(),
        )?;
        let program = crate::blastinput::app::Program::Blastp;
        let out = format.side_out(args.out.as_deref());
        // A format without a report: NCBI's formatting handler refuses SAM, AIRR and FASTA here;
        // the XInclude formats without a file come after `Query is Empty!`, in `search`
        // (`xinclude_format`). A value left to the formatting handler
        // (`app::EarlyFormat::Deferred`) fails in `check_options`, after the query options.
        if let Some(choice) = choice.parsed() {
            if crate::blastinput::app::report_format(choice, program, out)?.is_none() {
                crate::blastinput::app::formatting_handler_check(choice, program)?;
            }
        }
        // The options handlers of every format (`check_options`, also through `resolve`) make
        // the formatting handler's checks, the custom delimiter's among them.
        resolved_formats.push(if index == 0 {
            requested.check_options(&mut *outputs.diagnostics)?
        } else {
            requested.resolve()?
        });
    }
    // The numbers of descriptions and alignments belong to the format (only the pairwise
    // format reads them); the hit list size is the search's and must be one.
    let counts: Vec<crate::blastinput::report_counts::ReportCounts> = resolved_formats
        .iter()
        .map(|resolved| crate::blastinput::report_counts::ReportCounts {
            hitlist_size: resolved.max_target_seqs,
            num_descriptions: resolved.num_descriptions,
            num_alignments: resolved.num_alignments,
        })
        .collect();
    let numbers: Vec<i32> = outputs
        .formats
        .iter()
        .map(|format| {
            crate::blastinput::app::parse_formatting_string(format.outfmt).map(|c| c.number)
        })
        .collect::<Result<_>>()?;
    let shared = crate::blastinput::report_counts::common_report_counts(&numbers, &counts)?;
    // The hit lists are cut once for every report of the search (`report_prep::report_values`),
    // so its formats' cuts must be one (`BlastpApp::list_cut`): the library searches formats of
    // other cuts apart (`api::run::search_groups`).
    if numbers
        .iter()
        .any(|&number| args.app.list_cut(number) != args.app.list_cut(numbers[0]))
    {
        bail!(
            "output formats that cut the hit lists differently in one search (outfmt 11 keeps the search's lists, psiblast's outfmt 5 cuts them before each round's PSSM) are not supported by LOSAT's {}; request them in separate runs",
            args.app.name()
        );
    }
    let search_options = |resolved: &ResolvedBlastpArgs| {
        let mut options = resolved.clone();
        options.outfmt.clear();
        options.num_descriptions = 0;
        options.num_alignments = 0;
        // The `-html` links belong to the format as well: only the reports of outfmt 0 to 4 read
        // them (`blastinput::app::html_links`), so a run of outfmt 0 and 6 has them for 0 alone.
        options.html = None;
        format!("{options:?}")
    };
    let mut resolved_formats = resolved_formats.into_iter();
    let mut first = resolved_formats
        .next()
        .ok_or_else(|| anyhow::anyhow!("no output format was requested"))?;
    let first_options = search_options(&first);
    for resolved in resolved_formats {
        if search_options(&resolved) != first_options {
            bail!("the requested output formats resolve to different search options");
        }
        // The run's `-html` links are those of the format that reads them (outfmt 0 to 4, in
        // any place among the formats).
        first.html = first.html.or(resolved.html);
    }
    first.num_descriptions = shared.num_descriptions;
    first.num_alignments = shared.num_alignments;
    Ok(first)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:173-183
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:310-314
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:277-289
// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:131-137
/// psiblast's batches up to the first query whose Karlin-Altschul parameters cannot be computed:
/// blastp's `CLocalBlast` does not search such a batch (`CheckInternalData`, local_blast.cpp:
/// 177-180), but psiblast runs its preliminary search, which fails with `BLAST engine error:
/// Could not calculate ungapped Karlin-Altschul parameters ...` (exit 3) after the reports of the
/// queries before and the reader's messages of the batch. A query long enough to be split into
/// chunks (`split_protein_batch`) takes the split search instead, which ignores the chunks' error,
/// and NCBI's psiblast crashes in the traceback (exit 139; 2.17.0, 20,000 X letters at the default
/// CHUNK_SIZE): LOSAT rejects it. NCBI's psiblast also crashes, after the prolog, on such a PSSM
/// query (`-in_pssm`, which is never split): LOSAT rejects that too. Each batch of `batches` holds
/// one query
/// (`read_psiblast_query_batches`); its validity is that of the search's contexts (the query with
/// its SEG masks, or its letters with `-soft_masking true`: `blastp_query_frames`), and for a PSSM
/// query also that of `kbp_psi`, which a negative ungapped value of its file makes invalid (NCBI
/// 2.17.0: the same error, exit 3; line PS DECISIONS 30).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_parameters.c:183-184
/// Before that, the setup of the first batch's search reads the length of the database's first
/// sequence when the database has no letters: a database without sequences has none, and
/// psiblast stops with `BLAST Database error: OID not found` (exit 2), after the prolog, whatever
/// its query (NCBI 2.17.0 with `tk.faa` and an all-X query; probes `n.ps.empty.*` of the PS5c review;
/// blastp's `CLocalBlast` fails earlier, `run_call`; line PS DECISIONS 30).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:932-936
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbimpl.cpp:425
/// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:199-201
/// A database whose sequence file cannot be opened fails the first read of the preliminary
/// search of a valid query in the same way: psiblast's `CBlastPrelimSearch` runs without
/// `CLocalBlast`, whose catch makes it `Database memory map file error` for blastp
/// (local_blast.cpp:229-252, cited at `seqdb::SeqDbError`), so the reader's own text reaches
/// the application (`BLAST Database error: Error: File (<path>.psq) not found.`, exit 2, after
/// the prolog; NCBI 2.17.0, probe `y.g.ps.tk.nopsq.0` of the PS5c review; line PS DECISIONS 31).
///
/// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbfile.cpp:123-129
/// LOSAT's refusal of a PSSM search whose query letters give no Karlin-Altschul parameters
/// (`kbp_std`): NCBI's psiblast crashes where `PsiBlastSetupScoreBlock` reads the `kbp_psi` block
/// that `Blast_ScoreBlkKbpUngappedCalc` did not make (2.17.0, exit 139): a PSSM query
/// (`psiblast_stop_at_invalid_query`), and a later round whose masks leave no letter
/// (`psi_round::run_queries`).
const INVALID_PSSM_QUERY_REFUSAL: &str = "a PSSM query whose letters give no Karlin-Altschul parameters (with which NCBI BLAST+'s psiblast crashes) is not supported by LOSAT's PSIBLAST";

fn psiblast_stop_at_invalid_query(
    args: &ResolvedBlastpArgs,
    searched_queries: &[FastaRecord],
    database: Option<&database::BlastpDatabase>,
    mut batches: crate::blastinput::query_batch::QueryBatches,
) -> crate::blastinput::query_batch::QueryBatches {
    // The PSSM search's own stops (`PsiSearch::stop_before_search`: NCBI's check of the PSSM, its
    // scores, the query's set-up, blast_aux_priv.cpp:149-151) come before the score block whose
    // Karlin-Altschul parameters and the database read that this function checks
    // (`psiblast_stop_at_pssm_query` makes them; a query of invalid residues only has no valid
    // context here, M4a-1c E probe `old_ps_lower_ncbieaa.ps.strat`).
    if args.psi_search.as_deref().is_some_and(|psi| {
        psi.stop_before_search(args.strategy.is_some(), database.is_some())
            .is_some()
    }) {
        return batches;
    }
    let frames: Vec<Vec<QueryFrame>> = batches
        .starts
        .iter()
        .flat_map(|&(_, query)| {
            blastp_query_frames(args, std::slice::from_ref(&searched_queries[query]))
        })
        .collect();
    let contexts = pssm_hooks::psi_contexts(args, &frames);
    // `kbp_psi` of a PSSM query, valid as `s_BlastKarlinBlkIsValid` holds it
    // (`PsiSearch::ungapped_karlin_valid`).
    let kbp_psi_valid = |context: &QueryContext| {
        args.psi_search
            .as_deref()
            .is_none_or(|psi| psi.ungapped_karlin_valid(&context.karlin_params))
    };
    // The first read of the database (`PsiSearchReadFailure`, PSI-TBLASTN's too): no OID, or a
    // sequence file that cannot be opened.
    use crate::algorithm::common::database::PsiSearchReadFailure;
    let read_failure = database.and_then(|database| PsiSearchReadFailure::of(database.seqdb()));
    // NCBI's order (`PsiSearchReadFailure::precedes_setup_checks`): the empty database before the
    // PSSM's Karlin-Altschul check, for a query with a valid context (the PSSM query whose letters
    // give no Karlin-Altschul parameters crashes NCBI before, below); a sequence file that cannot
    // be opened stops the first batch when its query and `kbp_psi` are valid (the first read of the
    // preliminary search; volumes in the order NCBI reads them).
    let read_stop = match (&read_failure, contexts.first()) {
        (Some(failure), Some(context)) if failure.precedes_setup_checks() => {
            args.psi_search.is_none() || context.is_valid
        }
        (Some(_), Some(context)) => context.is_valid && kbp_psi_valid(context),
        _ => false,
    };
    let first_stop = if read_stop {
        Some(0)
    } else {
        contexts
            .iter()
            .position(|context| !context.is_valid || !kbp_psi_valid(context))
    };
    let Some(batch) = first_stop else {
        return batches;
    };
    let (input_start, searched_start) = batches.starts[batch];
    let span = batches.searched_span(searched_start);
    let query_length = searched_queries[searched_start].seq().len();
    let error = match protein_query_split_sizes(10000) {
        _ if read_stop => read_failure
            .expect("the first read's failure")
            .error(args.app.name()),
        // A PSSM query (`-in_pssm`, never split: `search_split_query_batches`) whose letters give no
        // Karlin-Altschul parameters: NCBI's psiblast crashes (2.17.0, exit 139) where
        // `PsiBlastSetupScoreBlock` reads the `kbp_psi` block that `Blast_ScoreBlkKbpUngappedCalc`
        // did not make (psiblast_aux_priv.cpp:89-107, blast_stat.c:2783-2799).
        Ok(_) if args.psi_search.is_some() && !contexts[batch].is_valid => {
            anyhow::anyhow!(INVALID_PSSM_QUERY_REFUSAL)
        }
        // Its letters give them, but the file's ungapped values make `kbp_psi` invalid.
        Ok(_) if args.psi_search.is_some() => crate::blastinput::app::engine_error(
            crate::report::query_warnings::INVALID_QUERY_MESSAGE,
        ),
        Ok(sizes)
            if crate::algorithm::blastp::protein_query_split::split_protein_batch(
                &[query_length],
                sizes,
            )
            .is_none() =>
        {
            crate::blastinput::app::engine_error(
                crate::report::query_warnings::INVALID_QUERY_MESSAGE,
            )
        }
        Ok(_) => anyhow::anyhow!(
            "a query whose Karlin-Altschul parameters cannot be computed and that is split into chunks (with which NCBI BLAST+'s psiblast crashes) is not supported by LOSAT's PSIBLAST"
        ),
        Err(error) => error,
    };
    let warnings = std::mem::take(&mut batches.warnings[searched_start]);
    batches.starts.truncate(batch);
    batches.input_done = input_start;
    batches.searched_done = searched_start;
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:337-340
    // NCBI made the batch's query factory and wrote its search strategy request
    // (`SaveSearchStrategy` at the start of `DoIterations`) before the search that fails
    // (`QueryBatches::stop_batch`; NCBI 2.17.0 `p9.ps.err.x_pol`, `p9.ps.err.nopsq`); a refusal of
    // LOSAT's writes none.
    batches.stop_batch = error.is::<crate::cli::NativeError>().then_some(span);
    batches.stop = Some((warnings, error));
    // The stop of a batch before LOSAT's refusal of a later batch's two queries with one id
    // replaces it (`QueryBatches::scope_refusal`): NCBI 2.17.0 stops at the query without
    // Karlin-Altschul parameters (line WQA D10 round 4, 2026-11-20, `r3-psiblast-probe`).
    batches.scope_refusal = false;
    batches
}

/// The stops of the query batches of LOSAT's refusals of two queries of one scope with one id
/// (`QueryBatches::resolve_refusals`), in the one form of every program
/// (`query_batch::refused_batch_error`; line WQA D10 round 5, 2026-11-20): with
/// `-export_search_strategy` the refusal (no request of the batch is written: `stop_batch` stays
/// `None`), else NCBI's error of the set-up of the batch's records that are first in the scope
/// (`refused_batch_setup_error`), else, where that set-up passes and the run's XML failure for its
/// matrix name takes the refusal (`ScopeRefusal::xml_takes`), the batch's search, else the refusal.
/// NCBI's errors of the batches before come first: their searches fail before the stop is acted
/// on.
fn refused_batch_stop(
    args: &ResolvedBlastpArgs,
    database: Option<&database::BlastpDatabase>,
    subject_records: &[FastaRecord],
    searched: &[FastaRecord],
    input: &seq_range::QueryInput,
    mut batches: crate::blastinput::query_batch::QueryBatches,
) -> crate::blastinput::query_batch::QueryBatches {
    use crate::blastinput::query_batch::BatchSetUp;
    batches.resolve_refusals(
        args.strategy.is_some(),
        |(input_start, searched_start), first| {
            // psiblast's batch holds the repeated query alone; rpsblast's set-up reads its
            // database's files first (`rps_hooks::RpsEngineError`): both refuse.
            if args.app != BlastpApp::Blastp {
                return BatchSetUp::Unknown;
            }
            let batch: Vec<usize> = (searched_start..searched.len())
                .take_while(|&query| input.ordinal(query) < input_start + first.len())
                .collect();
            let lengths: Vec<usize> = batch
                .iter()
                .map(|&query| searched[query].seq().len())
                .collect();
            let queries: Vec<FastaRecord> = batch
                .iter()
                .filter(|&&query| first[input.ordinal(query) - input_start])
                .map(|&query| searched[query].clone())
                .collect();
            refused_batch_setup_error(args, database, subject_records, &lengths, &queries)
        },
    );
    batches
}

/// BLASTP's set-up of the query records `queries` of a batch whose searched queries have the
/// lengths `batch_lengths` (`refused_batch_stop`): NCBI's error where a context of theirs is valid
/// and the set-up fails for the subjects, as `run_call` fails the first searched batch (a database
/// without sequences, `Database memory map file error`; the average subject length against
/// `-searchsp` or `OLD_FSC`), or where the search's first read of the database then fails (a
/// sequence file that cannot be opened, `BlastpDatabase::sequence_file_read_error`: `Database
/// memory map file error`; line WQA D10 round 6, review 4 M2); `Searched` where it passes with a
/// valid context; `Unknown` where no context is valid, and for a batch that NCBI splits into
/// chunks (each chunk is set up apart, `search_split_query_batches`), which is refused.
fn refused_batch_setup_error(
    args: &ResolvedBlastpArgs,
    database: Option<&database::BlastpDatabase>,
    subject_records: &[FastaRecord],
    batch_lengths: &[usize],
    queries: &[FastaRecord],
) -> crate::blastinput::query_batch::BatchSetUp {
    use crate::blastinput::query_batch::BatchSetUp;
    let split = protein_query_split_sizes(10000).map_or(true, |sizes| {
        crate::algorithm::blastp::protein_query_split::split_protein_batch(batch_lengths, sizes)
            .is_some()
    });
    if split {
        return BatchSetUp::Unknown;
    }
    // `call_setup::call_lookup`'s validity: a context without letters is invalid.
    let contexts = pssm_hooks::psi_contexts(args, &blastp_query_frames(args, queries));
    if !contexts
        .iter()
        .any(|context| context.is_valid && context.aa_len > 0)
    {
        return BatchSetUp::Unknown;
    }
    if database.is_some_and(database::BlastpDatabase::first_oid_read_fails) {
        return BatchSetUp::Fails(crate::blastinput::app::engine_error(
            "Database memory map file error",
        ));
    }
    let encoded: Vec<EncodedProtein> = match database {
        Some(_) => Vec::new(),
        None => subject_records
            .iter()
            .map(|record| encode_protein_sequence(record.seq()))
            .collect(),
    };
    let Ok(extent) = call_setup::db_extent(args, database, &encoded) else {
        return BatchSetUp::Unknown;
    };
    if crate::stats::setup::setup_subject_length(
        extent.total_db_len,
        extent.num_subjects,
        (!args.old_fsc).then_some((extent.cutoff_subject_length, false)),
    )
    .is_none()
    {
        return BatchSetUp::Fails(crate::blastinput::app::engine_error(
            "The average subject length is too short",
        ));
    }
    let Some(read_error) =
        database.and_then(|database| database.sequence_file_read_error("BLASTP"))
    else {
        return BatchSetUp::Searched;
    };
    BatchSetUp::Fails(
        match read_error.downcast::<database::PreliminaryReadError>() {
            Ok(read_error) => crate::blastinput::app::engine_error(read_error.0),
            Err(error) => error,
        },
    )
}

/// psiblast: LOSAT's refusal of a subject set (`-subject`) none of whose sequences has letters
/// (records without residues, `-subject_loc` intervals that start just past their records' ends).
/// NCBI's psiblast checks its subjects as one sequence block when it sets them up, after the
/// options are processed (the PSSM read too) and before the query is read or anything is written
/// (`InitializeSubject`, psiblast_app.cpp:444-445). The block of such a set is not made, and the
/// check reads through a null pointer: NCBI 2.17.0 exits 139 with empty standard output and error
/// (with `-query` and with `-in_pssm`; one empty record, two, or a cut to no letters), while a set
/// with one sequence that has letters runs with the warning of each empty one, as blastp runs
/// every set (`CBl2Seq`, without this check).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_db_adapter.cpp:73-76
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_aux_priv.cpp:651-663
#[cold]
#[inline(never)]
fn psiblast_check_subject_letters(subjects: &[FastaRecord]) -> Result<()> {
    if !subjects.is_empty() && subjects.iter().all(|record| record.seq().is_empty()) {
        bail!(
            "a subject set whose sequences all have no letters (with which NCBI BLAST+'s psiblast crashes) is not supported by LOSAT's PSIBLAST"
        );
    }
    Ok(())
}

/// psiblast `-in_pssm`: what stops the run at its one query batch before the PSSM's search, after
/// the outfmt 0 prolog, with NCBI's error (`PsiSearch::stop_before_search`, the rule tblastn's
/// `-in_pssm` shares): with `-export_search_strategy` and `-subject` (no database, `database`) the
/// search strategy's "Remote PSI-BL2SEQ is not supported" (the options error, exit status 1; without
/// the strategy NCBI runs the search), else a PSSM that NCBI's `CPsiBlastValidate::Pssm` refuses
/// (the engine error; the PSSM is checked before the query's Karlin-Altschul parameters,
/// `psiblast_stop_at_invalid_query`), else the PSSM engine's error for the scores of a PSSM of
/// frequency ratios only. The first two write no request (`QueryBatches::stop_batch` cleared); the
/// last stops after the request (`stop_batch`: the PSSM query's batch), as NCBI's
/// `SaveSearchStrategy` comes before `CPsiBlast` (psiblast_app.cpp:340; the local rounds make
/// `CPsiBlast` later, psiblast_app.cpp:239).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/psiblast_impl.cpp:58-67 (cited at
/// `PsiSearch::stop_before_search`)
#[cold]
#[inline(never)]
fn psiblast_stop_at_pssm_query(
    args: &ResolvedBlastpArgs,
    database: Option<&database::BlastpDatabase>,
    mut batches: crate::blastinput::query_batch::QueryBatches,
) -> crate::blastinput::query_batch::QueryBatches {
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:537-540 (the object manager's
    // index of the query's annotations, `scope->AddTopLevelSeqEntry`, before `DoIterations` and its
    // `SaveSearchStrategy` and `CPsiBlast`)
    // The PSSM query's batch carries the object manager's messages (`PsiSearch::annot`); a
    // diagnostic stream that they fail silences the run's later messages from here on, and an
    // exception of the index stops the run at its query, after the prolog and those messages,
    // with no search strategy's request.
    if let Some(psi) = args.psi_search.as_deref() {
        psi.fail_diag_stream();
        if let Some(abort) = psi.annot_abort() {
            // Before the search's own stops (`psiblast_stop_at_invalid_query`, which took the
            // batch's messages already).
            let warnings = match batches.stop.take() {
                Some((warnings, _)) => warnings,
                None => batches
                    .warnings
                    .first_mut()
                    .map(std::mem::take)
                    .unwrap_or_default(),
            };
            batches.starts.clear();
            batches.input_done = 0;
            batches.searched_done = 0;
            batches.stop = Some((warnings, abort));
            batches.stop_batch = None;
            // The query's ids never join the scope (`PsiSearch::joins_scope`): no refusal of
            // another record with one of them (NCBI BLAST+ 2.17.0, line RG D3-1 fix 4 d01-d06).
            batches.scope_refusal = false;
            return batches;
        }
    }
    let Some(stop) = args
        .psi_search
        .as_deref()
        .and_then(|psi| psi.stop_before_search(args.strategy.is_some(), database.is_some()))
    else {
        return batches;
    };
    let request = stop
        .after_request()
        .then(|| batches.first_request_end().map(|end| 0..end))
        .flatten();
    let warnings = batches
        .warnings
        .first_mut()
        .map(std::mem::take)
        .unwrap_or_default();
    batches.starts.clear();
    batches.input_done = 0;
    batches.searched_done = 0;
    batches.stop = Some((warnings, stop.error()));
    batches.stop_batch = request;
    batches.scope_refusal = false;
    batches
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:121, :196 (cited at
// `blast_xml2::reports_held`)
/// Where blastp's (and psiblast's) run of the formats of `outputs` stops for a matrix name that
/// NCBI's BLAST XML report data does not know (`blast_xml::matrix::xml_formatting_failure`), with
/// the facts of the search: no database (`-subject`, `m_IsBl2Seq`) and psiblast (`m_IsIterative`).
/// The one decision behind the search strategy's last request (`search`) and the formats' error
/// (`run_call`, at the run's query of `blastp_xml_failing_query`).
fn blastp_xml_matrix_failure(
    args: &ResolvedBlastpArgs,
    outputs: &ReportOutputs<'_>,
    database: bool,
) -> Option<crate::report::blast_xml::matrix::XmlMatrixFailure> {
    crate::report::blast_xml::matrix::xml_formatting_failure(
        outputs.formats.iter().map(|format| format.outfmt),
        &args.matrix_name,
        !database,
        args.app == BlastpApp::Psiblast,
    )
}

/// The run's query at which blastp's formats of `outputs` fail for the matrix name
/// (`blastp_xml_matrix_failure`), of the run's scopes of the batches added (`scopes`, under
/// `-parse_deflines`): `QueryScopes::xml_failing_query`.
fn blastp_xml_failing_query(
    args: &ResolvedBlastpArgs,
    outputs: &ReportOutputs<'_>,
    database: bool,
    scopes: Option<&crate::blastinput::query_batch::QueryScopes>,
) -> Option<crate::blastinput::query_batch::XmlFailingQuery> {
    blastp_xml_matrix_failure(args, outputs, database).map(|failure| {
        crate::blastinput::query_batch::QueryScopes::xml_failing_query(scopes, failure)
    })
}

/// The search of `run_local` and of the CLI, after NCBI's processing of the options and
/// the reading of the subjects: `Query is Empty!` for records without a query, the query
/// batch size, the subjects without letters, the XInclude formats (`xinclude`, the CLI's
/// format), LOSAT's limits and environment, then the query batches
/// (`read_query_batches`), the search of the batches that NCBI searches and their reports.
#[allow(clippy::too_many_arguments)]
fn search(
    mut args: ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    reading: QueryReading,
    subjects: &SubjectSet<'_>,
    database: Option<&database::BlastpDatabase>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    xinclude: Option<&crate::blastinput::app::FormatChoice>,
) -> Result<()> {
    // What the search reads its subjects from: the database of `-db` (then `subjects` has
    // no record), or the subject records.
    let subject_source = match database {
        Some(database) => BlastpSubjectSource::Database(database),
        None => {
            // LOSAT_X_DBIDXREAD acts on a database's index files only; the new side is marked
            // where a volume is opened (`seqdb::Volume::open`).
            if crate::seqdb::x_dbidxread() != crate::seqdb::XDbIdxRead::Off {
                crate::utils::xreach::x_reach!("LOSAT_X_DBIDXREAD", false);
            }
            // LOSAT_X_BPDBKEEP acts on the preliminary stage of `-db` (`database/x_db_keep.rs`).
            if crate::seqdb::x_bp_db_keep() {
                crate::utils::xreach::x_reach!("LOSAT_X_BPDBKEEP", false);
            }
            BlastpSubjectSource::Fasta(&subjects.records)
        }
    };
    // NCBI reads the subjects before the queries (`run`); the records are those of NCBI's
    // reader (`fasta_reader`; `run_local`'s callers read them so). The CLI has already
    // reported an empty query file and read its queries; the records of the other callers
    // get NCBI's warning for an input without records.
    if query_records.is_empty() && matches!(reading, QueryReading::Records) {
        writeln!(
            outputs.diagnostics,
            "Warning: [{}] Query is Empty!",
            args.app.ncbi_name()
        )?;
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:215-216
    // The batch size is read before the formatter is made and before the queries, which
    // NCBI reads one batch at a time after the outfmt 0 prolog (`read_query_batches`).
    let batch_size = args.query_batch_size()?;
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:140-153 (cited at
    // `blastinput::db_mask::DbMaskOption::resolve`)
    // The formatter resolves the database masking (`-db_soft_mask`, `-db_hard_mask`; series NP,
    // M4b) before its other work (`shared::resolve_database_masking`).
    if let (Some(database), Some(option)) = (database, args.db_mask.as_ref()) {
        let diagnostics = &mut *outputs.diagnostics;
        let archive = args.archive.as_ref();
        args.db_mask_resolved = crate::algorithm::common::database::resolve_database_masking(
            database.seqdb(),
            database.name(),
            option,
            args.app.name(),
            archive,
            args.strategy.as_ref(),
            &mut |text| match archive {
                Some(handle) => handle.write_warning(diagnostics, text),
                None => diagnostics.write_all(text),
            },
        )?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:225-227
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_setup_cxx.cpp:773-789
    // The formatter, made after the query batch size and before the outfmt 0 prolog, sets
    // up the subjects: a subject without letters (a record without residues, or an
    // interval that starts just past its record's end) gets its warning, stays in the
    // database statistics with no letters, and is never searched.
    crate::blastinput::input_files::write_empty_subject_warnings(
        &subjects.records,
        args.app.ncbi_name(),
        outputs.diagnostics,
    )?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:250-253
    // An XInclude format written to standard output fails after the formatter is made.
    // (`run_local` checks its formats with the options, `resolve_for_formats`.)
    if let Some(choice) = xinclude {
        crate::blastinput::app::xinclude_check(choice)?;
    }
    // `PrintProlog`: outfmt 20's formatter warns here when its field list asks for a name of the
    // taxonomy database that it cannot find (`tabular::write_tabular_prolog_warning`); its line
    // of the field names comes with the first report (or `write_pairwise_prologs`), after LOSAT's
    // refusals.
    for format in outputs.formats.iter() {
        if let Some(spec) = super::tabular::custom_tabular_spec(format.outfmt, args.app)
            .filter(|spec| spec.checks_taxdb())
        {
            let taxdb_found = match database {
                Some(database) => database.taxdb_found(),
                None => args.subject_taxdb_found,
            };
            super::tabular::write_tabular_prolog_warning(
                spec,
                taxdb_found,
                args.app.ncbi_name(),
                &mut *outputs.diagnostics,
            )?;
        }
    }
    // The query range applies to every query record as NCBI's batch reader reads it
    // (`seq_range::cut_queries`): a record is searched cut, and one whose interval starts
    // just past its end is a query without data.
    let (searched_queries, query_input) = match args.query_range.as_ref() {
        Some(range) => {
            let ranged = seq_range::cut_queries(query_records, range);
            (std::borrow::Cow::Owned(ranged.records), ranged.input)
        }
        None => (
            std::borrow::Cow::Borrowed(query_records),
            seq_range::QueryInput::whole(query_records),
        ),
    };
    // LOSAT's refusals of the masked search (`ResolvedMask::check_losat_limits`) come where NCBI
    // starts the search, after its checks.
    if let (Some(database), Some(mask)) = (database, args.db_mask_resolved.as_ref()) {
        mask.check_losat_limits(
            database.seqdb(),
            args.app.name(),
            crate::blastinput::db_mask::HardScan::Protein,
        )?;
    }
    check_blastp_environment(args.app, database.is_none())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259
    // What the run writes before its query batches ends here.
    if let Some(observer) = outputs.observer.as_deref_mut() {
        observer.queries_begin();
    }
    // Under `-parse_deflines`, the scopes of the batches' queries (`QueryScopes`); with
    // `-in_pssm`, of the PSSM's query when the reports show its ids (`PsiSearch::query_record`).
    let believed_queries = crate::blastinput::query_batch::believed_queries(
        args.parse_deflines,
        args.psi_search.is_some(),
        query_records,
    );
    let query_scopes = crate::blastinput::query_batch::QueryScopes::of_search(
        believed_queries,
        database.is_none(),
        believed_queries
            && crate::blastinput::query_batch::outfmts_have_xml(
                outputs.formats.iter().map(|format| format.outfmt),
            ),
        args.app.name(),
    )
    .map(|scopes| {
        scopes
            .with_refusal(args.scope_refusal.clone())
            .with_xml_matrix_failure(blastp_xml_matrix_failure(
                &args,
                outputs,
                database.is_some(),
            ))
    });
    let batches = match args.app {
        // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:158-159
        // NCBI never splits an RPS-BLAST query, so it does not read OVERLAP_CHUNK_SIZE for it
        // (`SplitQuery_ShouldSplit`, split_query_aux_priv.cpp:73-97).
        BlastpApp::Blastp | BlastpApp::Rpsblast => {
            crate::blastinput::query_batch::read_query_batches(
                query_records,
                &searched_queries,
                &query_input,
                reading,
                batch_size,
                args.app != BlastpApp::Rpsblast,
                query_scopes,
            )
        }
        BlastpApp::Psiblast => psiblast_stop_at_invalid_query(
            &args,
            &searched_queries,
            database,
            crate::blastinput::query_batch::read_psiblast_query_batches(
                query_records,
                &searched_queries,
                &query_input,
                reading,
                query_scopes,
            ),
        ),
    };
    // psiblast `-in_pssm`: the object manager's stop as the query joins the scope, NCBI's check of
    // the PSSM, or the search strategy's request without a database, stops the run at its query
    // (`psiblast_stop_at_pssm_query`); a query that never joins NCBI's scope meets no other record
    // there.
    // LOSAT's refusal of two queries of one scope with one id is the stop of their batch
    // (`QueryBatches::scope_refusal`): the run writes the reports of the batches before, unless
    // NCBI stops at one of them first (line WQA D10 round 4, 2026-11-20).
    let batches = psiblast_stop_at_pssm_query(&args, database, batches);
    // The batch of LOSAT's refusal stops with NCBI's error of its set-up where NCBI's set-up of the
    // batch fails for the subjects (`refused_batch_stop`; line WQA D10 round 5, 2026-11-20).
    let batches = refused_batch_stop(
        &args,
        database,
        &subjects.records,
        &searched_queries,
        &query_input,
        batches,
    );
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-262
    // LOSAT's refusal of a strategy file on a report stream comes where NCBI writes its first
    // request: once its first batch is read, before it is searched (`first_request_end`; this
    // search writes the requests where its call ends, after the searches, which would hold the
    // reports). A batch that stops with LOSAT's refusal of an id is read first too.
    if batches.first_request_end().is_some() {
        crate::report::archive::strategy::check_report_stream(args.strategy.as_ref())?;
    }
    // LOSAT's limits come where NCBI starts the search, after its checks: a run that stops at its
    // first batch (`batches.starts` empty: a reader's error, NCBI's check of a PSSM, the object
    // manager's stop as a PSSM's query joins the scope) never searches, and gives NCBI's error
    // (NCBI BLAST+ 2.17.0 with `-evalue 1e999`, `-use_sw_tback`, `-gapopen 9 -gapextend 2` and
    // `-matrix BLOSUM45`: line RG D3-1 fix 4 `r4/probe`). NCBI crashes with IDENTITY's compressed
    // lookup table once it searches a query batch, so that refusal follows the same rule.
    if !batches.starts.is_empty() {
        validate_requested_blastp_support(&args)?;
        check_identity_compressed_lookup_of(&args)?;
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513-523
    // psiblast's rounds of FASTA queries (`-num_iterations` other than 1): the round driver,
    // one search per round per query (`psi_round`).
    if let Some(plan) = args.psi_rounds.clone() {
        return psi_round::search_rounds(
            &args,
            &plan,
            query_records,
            &searched_queries,
            &query_input,
            batches,
            subjects,
            subject_source,
            subject_label,
            outputs,
        );
    }
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:294-295
    // A run that stops at a query batch writes the reports of the batches before it
    // without the epilog.
    let epilog = batches.stop.is_none();
    // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:274
    // IDENTITY's reset of the composition mode (`ResolvedBlastpArgs::identity_cbs_reset`) writes
    // through the options memento into the options that every query batch reuses: the first
    // batch resets the mode and warns its queries, and the later batches start in mode 0 without
    // the warning. These are the searched queries of the first batch.
    let identity_warned = if args.identity_cbs_reset {
        batches
            .starts
            .get(1)
            .map_or(batches.searched_done, |&(_, searched)| searched)
    } else {
        0
    };
    // The searches and reports of the batches up to the input record `input_done` and
    // the searched query `searched_done`; `stopped`: the PSSM files of a PSSM's stop at the last
    // of them (`BlastpReports::stopped_pssms`).
    let search_batches_once = |args: ResolvedBlastpArgs,
                               outputs: &mut ReportOutputs<'_>,
                               pool: &crate::utils::threading::SearchPool<'_>,
                               input_done: usize,
                               searched_done: usize,
                               epilog: bool,
                               stopped: Option<crate::algorithm::psiblast::pssm::PssmFiles>,
                               next: Option<NextRequest<'_>>|
     -> Result<()> {
        let done_input = seq_range::QueryInput {
            placements: seq_range::Placements(
                query_input
                    .placements
                    .0
                    .iter()
                    .take(searched_done)
                    .copied()
                    .collect(),
            ),
            ordinals: query_input.ordinals[..searched_done].to_vec(),
            input_lengths: query_input.input_lengths[..input_done].to_vec(),
            skipped: query_input.skipped[..input_done].to_vec(),
            range_length: query_input.range_length,
        };
        let ranges = BlastpRanges {
            input_records: &query_records[..input_done],
            query_input: &done_input,
            subject_placements: &subjects.placements,
            x_ordinal_base: 0,
            searched_base: 0,
            split_search_spaces: None,
        };
        let stopped = std::sync::Mutex::new(stopped);
        let reports = BlastpReports {
            reader_warnings: &batches.warnings[..searched_done],
            identity_warned,
            epilog,
            x_prolog: true,
            x_queries_before: 0,
            x_xml2_carry: None,
            last_call: true,
            psi_round: None,
            stopped_pssms: Some(&stopped),
            query_scopes: batches
                .query_scopes
                .as_ref()
                .map(|scopes| scopes.numbered_from(0)),
            prefetch_scope: None,
            strategy_next: next,
            call_requests: None,
            x_plan_blocks: true,
        };
        run_resolved_in_pool(
            args,
            &searched_queries[..searched_done],
            subject_source,
            "",
            subject_label,
            outputs,
            pool,
            None,
            &ranges,
            &reports,
            x_lookup_ahead::XAheadBlock::none(),
        )
    };
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:513-523
    // psiblast makes each query's PSSM after its report: a PSSM that fails
    // (`psiblast::pssm::PssmStop`) stops the run there, after the reports of the queries up to
    // it (without the epilog) and the PSSMs of those before it. LOSAT searches all queries
    // together, so those queries are searched again and reported alone.
    let search_batches = |args: ResolvedBlastpArgs,
                          outputs: &mut ReportOutputs<'_>,
                          pool: &crate::utils::threading::SearchPool<'_>,
                          input_done: usize,
                          searched_done: usize,
                          epilog: bool,
                          next: Option<NextRequest<'_>>|
     -> Result<()> {
        if args.psi_pssm.is_none() {
            return search_batches_once(
                args,
                outputs,
                pool,
                input_done,
                searched_done,
                epilog,
                None,
                next,
            );
        }
        // The search strategy's requests (`CallRequests`): the first search releases those of
        // every batch and of the batch after the call (`next`) when no PSSM stops it; a PSSM's stop
        // comes before the call marks its end, so the first search releases none, and the search
        // again of the queries up to the failing one releases theirs (NCBI wrote each query's
        // before its search, and none after the failing query's).
        let again_args = args.clone();
        match search_batches_once(
            args,
            outputs,
            pool,
            input_done,
            searched_done,
            epilog,
            None,
            next,
        ) {
            Err(error) if error.is::<crate::algorithm::psiblast::pssm::PssmStop>() => {
                let crate::algorithm::psiblast::pssm::PssmStop {
                    query,
                    error,
                    files,
                } = error
                    .downcast::<crate::algorithm::psiblast::pssm::PssmStop>()
                    .expect("the error was checked to be a PssmStop");
                let after = batches
                    .starts
                    .partition_point(|&(_, searched)| searched <= query);
                let input_end = batches
                    .starts
                    .get(after)
                    .map_or(input_done, |&(input, _)| input.min(input_done));
                // The queries up to the failing one, searched and reported again with the PSSMs
                // made before it (`BlastpReports::stopped_pssms`): the reads of each PSSM (up to
                // the failing query's) look their taxids up after the query's report, and their
                // messages follow it, before the next query's warnings (the failing query's after
                // its report), as where no PSSM stops (review 1 of line DBID pssmread, L1; review
                // 3, H1). The files are written after the reports, then the run stops with the
                // PSSM's error.
                search_batches_once(
                    again_args,
                    outputs,
                    pool,
                    input_end,
                    query + 1,
                    false,
                    Some(files),
                    None,
                )?;
                Err(error)
            }
            other => other,
        }
    };
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the queries in blocks of whole batches, each block a search
    // of its own whose report is written before the next block is searched
    // (`x_query_blocks.rs`, exactness argument in its header). These reports take it: tabular
    // format 6, outfmt 0 or 7 (the outfmt 0 prolog with the first block, the epilog and outfmt 7's
    // count of the queries with the last), outfmt 8, 9, 11 or 12, outfmt 5 (the header with the
    // first block, the end with the last), outfmt 13 to 16 (as outfmt 5, the framing state handed
    // from block to block, `Xml2Carry`), also several of them, the hit records and an observer
    // (`run_local`; not outfmt 11 or 13 to 16 in such a call, see `x_query_block_plan`). One block is the search of all queries below, which needs none of
    // the checks: a refusal in it has no block's bytes before it. Where a check fails, the search
    // of all queries below refuses as NCBI does. With `-db` the blocks are searched only when no
    // definition line of the database can fail a batch where its titles are read
    // (`database::x_titles_allow_blocks`, or `x_titles_read_allow_blocks` for outfmt 8, 9, 12 and
    // 5, `x_titles_taxids_allow_blocks` for 15 and 16 (the titles and the taxids in one pass), for
    // outfmt 0, 5, 6 and 7 with `taxdb.bti` on the search path also the taxids that each batch's
    // prefetch reads and looks up, read in the block search's pool before the first block), and
    // with outfmt 0 and `-subject` only when no subject's title is refused where the report shows
    // it (`x_block_shown_titles_allowed`), with outfmt 5, 15 and 16 only when no title of the
    // queries and the subjects is refused (`names::Refusal`; `x_query_blocks.rs`
    // `x_block_xml_*_allowed`).
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/data4xmlformat.cpp:213-218 (the BLAST
    // XML report data of the first result NCBI formats: the matrix of its formatting table,
    // found by the name as typed)
    // NCBI fails when it formats a batch's results in BLAST XML with a matrix name that its XML
    // data does not know ("unsupported score matrix", `run_resolved_in_pool`'s format loop): the
    // batch of the run's first searched query, or where the report is held
    // (`blastp_xml_matrix_failure`; `XmlMatrixFailure::SecondQueryOrEnd`) of NCBI's next query of
    // another Seq-id (`QueryScopes::xml_failing_query` of the batches' scopes); the
    // search strategy's file then holds the requests up to that
    // batch's (probe 7 of RG-F1, probe K4f of line RG), whichever call writes them. A held report
    // that fails at the run's end, or a run that stops before, writes every request it makes.
    if let Some(strategy) = args.strategy.as_ref() {
        if let Some(end) = blastp_xml_matrix_failure(&args, outputs, database.is_some())
            .and_then(|failure| batches.xml_matrix_request_end(failure))
        {
            strategy.stop_after(end);
        }
    }
    let x_plan = x_query_block_plan(
        &args,
        &searched_queries,
        &query_input,
        &batches,
        subject_source,
        outputs,
        batch_size,
    );
    if x_plan.is_none() {
        super::x_query_blocks::x_reach_block_search(false);
    }
    if let Some((x_kinds, blocks)) = x_plan {
        use super::x_query_blocks::{XBlockReport, XBlockSink, XOffsetObserver};
        let reader_warnings = &batches.warnings;
        let ReportOutputs {
            other_outfmts,
            formats,
            diagnostics,
            later_warnings,
            hits,
            observer,
            html_posts,
        } = &mut *outputs;
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the outputs' `-outfmt` values, in request order.
        let outfmts: Vec<&str> = formats.iter().map(|format| format.outfmt).collect();
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the reports the blocks write, whose titles are checked
        // below.
        let x_pairwise = x_kinds.contains(&XBlockReport::Pairwise);
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the hit records of all blocks, each with the index of its
        // query among the call's queries (the block's index plus the block's first query), handed to
        // the caller's sink once, after the last block: one call with the final list of all
        // queries, as the search of all queries makes it.
        let x_want_records = hits.is_some();
        let mut x_records: Option<Vec<PairwiseHit>> = x_want_records.then(Vec::new);
        // EXPERIMENT (LOSAT_X_BPQBLOCK): the framing state that the blocks of an outfmt 13 to 16
        // run hand on (such a run has one output, `x_query_block_plan`).
        let x_xml2_carry = std::sync::Mutex::new(crate::report::blast_xml2::Xml2Carry::default());
        let searched_in_blocks =
            crate::utils::threading::with_search_pool(args.num_threads, "blastp", |pool| {
                // The last checks of the block search, on the pool; when one fails, nothing has
                // been written and the search of all queries below runs: the lookup tables' words
                // (`x_word_blocks`, which also closes the blocks of a compressed table on their
                // words), then with `-db` the definition lines
                // (`database::x_titles_allow_blocks`; for outfmt 8, 9, 12 and 5
                // `database::x_titles_read_allow_blocks`; for outfmt 0, 5, 6 and 7 with `taxdb.bti`
                // found also the taxids and lookups of the batches' prefetch), with outfmt 0 and
                // `-subject` the
                // subjects' titles (`x_block_shown_titles_allowed`), and with outfmt 5 the
                // subjects' names (`x_block_xml_subjects_allowed`, `-subject`) and the queries'
                // titles (`x_block_xml_run_allowed`); with several outputs, the
                // checks of every kind of report among them.
                // A plan of one block of letters whose standard table's words pass the most
                // (`x_aa_one_block_splits`) is split here; one block stays the search of all queries.
                let Some((blocks, x_block_words)) = x_word_blocks(
                    &args,
                    &searched_queries,
                    &query_input,
                    batch_size,
                    &blocks,
                    pool,
                )
                .filter(|(blocks, _)| blocks.len() > 1) else {
                    super::x_query_blocks::x_reach_block_search(false);
                    return Ok(false);
                };
                // Outfmt 0, 5, 6 and 7 (psiblast's too) prefetch the subjects of each
                // batch, which with `taxdb.bti` on the search path reads their taxids and looks
                // them up in the taxonomy database (`ListedTitles::prefetch_messages`); a refusal
                // there would come in a later block, so their checks then also read every
                // sequence's taxid and ask whether any lookup can be refused
                // (`shared::taxdb_may_refuse`, without the scientific names, which the prefetch
                // does not ask for). The refusal of NCBI's taxonomy cache
                // (`shared::TaxdbRun::cache_collision`) can come at any batch that prefetches or
                // reads (`shared::TaxdbRun::cache_may_refuse`), in a later block: such a run takes
                // the search of all queries, which refuses before any row. So does a run whose first
                // format is a Seq-annot that the batches prefetch for
                // (`PRE_FETCH_SEQS_LIMIT=2147483647`, `shared::seq_annots_prefetch`), whose
                // prefetch can refuse a subject with parsed Seq-ids too
                // (`ListedTitles::check_prefetched`). A table of fields that read nothing of the
                // subject has the same refusals of the prefetch, only on a database with parsed
                // Seq-ids, whose definition lines the checks below never take
                // (`database::x_titles_allow_blocks`).
                if database.is_some()
                    && outfmts.first().is_some_and(|outfmt| {
                        crate::algorithm::common::database::seq_annots_prefetch(
                            crate::blastinput::app::parse_formatting_string(outfmt)
                                .ok()
                                .map(|choice| choice.number),
                        )
                    })
                {
                    super::x_query_blocks::x_reach_block_search(false);
                    return Ok(false);
                }
                let x_prefetch_taxids = |x_kind: XBlockReport| {
                    matches!(
                        x_kind,
                        XBlockReport::Tabular | XBlockReport::Pairwise | XBlockReport::Xml
                    ) && database.is_some_and(|database| database.taxdb_found())
                };
                // Outfmt 13 to 16 ask for the same refusal of NCBI's taxonomy cache at the first
                // batch that reports a subject with a positive taxid (`read_reported_taxonomy`),
                // which can be in a later block too: the search of all queries, which refuses
                // before any byte.
                let titles_allowed = |x_kind: XBlockReport| match (database, x_kind) {
                    (Some(database), _)
                        if (x_prefetch_taxids(x_kind)
                            || matches!(x_kind, XBlockReport::Xml2) && database.taxdb_found())
                            && database.cache_may_refuse() =>
                    {
                        false
                    }
                    // Outfmt 5 whose batches read the taxids: the check of outfmt 15 and 16.
                    // The prefetch asks for no scientific name.
                    (Some(database), XBlockReport::Xml) if x_prefetch_taxids(x_kind) => {
                        database::x_titles_taxids_allow_blocks(database, pool)
                    }
                    // Outfmt 8, 9, 12 and 5 write no tabular id: only the definition lines' reading
                    // and bytes can fail (`x_titles_read_allow_blocks`).
                    (Some(database), XBlockReport::Seqalign | XBlockReport::Xml) => {
                        database::x_titles_read_allow_blocks(database, pool)
                    }
                    // Outfmt 15 and 16 also read each subject's taxid and the taxonomy
                    // database's name of it, in the same pass over the headers
                    // (`x_titles_taxids_allow_blocks`).
                    (Some(database), XBlockReport::Xml2) => {
                        database::x_titles_taxids_allow_blocks(database, pool)
                    }
                    (Some(database), _) => database::x_titles_allow_blocks(
                        database,
                        x_pairwise,
                        x_prefetch_taxids(x_kind),
                        pool,
                    ),
                    (None, XBlockReport::Pairwise) => {
                        x_block_shown_titles_allowed(&subjects.records, pool)
                    }
                    (None, XBlockReport::Xml | XBlockReport::Xml2) => {
                        super::x_query_blocks::x_block_xml_subjects_allowed(&subjects.records, pool)
                    }
                    (None, _) => true,
                };
                if !x_kinds
                    .iter()
                    .enumerate()
                    .filter(|&(index, x_kind)| !x_kinds[..index].contains(x_kind))
                    .all(|(_, &x_kind)| titles_allowed(x_kind))
                {
                    super::x_query_blocks::x_reach_block_search(false);
                    return Ok(false);
                }
                // EXPERIMENT (LOSAT_X_BPQBLOCK): the XML report's refusals of the queries'
                // titles (`names::Refusal`), which it makes before its first byte.
                if x_kinds.iter().any(|x_kind| x_kind.is_xml())
                    && !super::x_query_blocks::x_block_xml_run_allowed(query_records, pool)
                {
                    super::x_query_blocks::x_reach_block_search(false);
                    return Ok(false);
                }
                super::x_query_blocks::x_reach_block_search(true);
                // EXPERIMENT (LOSAT_X_BPQBLOCK): the subjects encoded once for all blocks
                // (`x_query_blocks::XBlockSubjects`), on the search's pool.
                // Diagnostic (LOSAT_TIMING): the blocks' shared subject encoding for the common
                // line (utils/xtiming.rs). No NCBI counterpart.
                let t_encode = crate::utils::xtiming::clock();
                let mut x_subjects =
                    x_bp_encode_block_subjects(&subjects.records, args.num_threads, pool)?;
                crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, t_encode);
                // Each block's words, for its table's build (`x_query_blocks::x_block_appends`).
                x_subjects.x_set_block_words(
                    blocks
                        .iter()
                        .zip(&x_block_words)
                        .map(|(block, &words)| (&searched_queries[block.clone()], words)),
                );
                let pool = &pool.x_with_ext(&x_subjects);
                // The caller's receivers of the side files of outfmt 13 and 14, which every
                // block's output keeps.
                let side_sinks: Vec<_> = formats
                    .iter()
                    .map(|format| format.sink.side_files().cloned())
                    .collect();
                let mut sinks: Vec<XBlockSink<'_>> = formats
                    .iter_mut()
                    .map(|format| XBlockSink::new(&mut format.sink))
                    .collect();
                // The HSPs of the blocks before: the first index of the block's HSPs in the list
                // of all queries.
                let mut x_hsps_before: HspIndex = 0;
                // EXPERIMENT (LOSAT_X_BPLKPAHEAD): each block's search builds the next block's lookup
                // table beside it, on the pool (`blast_engine/x_lookup_ahead.rs`).
                let mut x_ahead = x_lookup_ahead::XLookupAhead::new(&args, pool, batch_size);
                for (index, block) in blocks.iter().enumerate() {
                    let block_queries = &searched_queries[block.clone()];
                    let x_next_queries = blocks
                        .get(index + 1)
                        .map(|next| &searched_queries[next.clone()]);
                    let block_input = seq_range::QueryInput::whole(block_queries);
                    let records_before = x_records.as_ref().map_or(0, Vec::len);
                    let x_first_query =
                        u32::try_from(block.start).expect("a query index of the call fits q_idx");
                    // The block's rows, shared with its reports until the block has ended.
                    let mut x_found = Vec::new();
                    let mut x_collect =
                        |found: std::sync::Arc<Vec<PairwiseHit>>| x_found.push(found);
                    // EXPERIMENT (LOSAT_X_BPQBLOCK): the caller's observer with the block's HSP
                    // and query indices moved past the blocks before.
                    let mut x_offset = observer.as_deref_mut().map(|inner| XOffsetObserver {
                        inner,
                        base: x_hsps_before,
                        first_query: block.start,
                    });
                    let mut block_outputs = ReportOutputs {
                        other_outfmts,
                        formats: sinks
                            .iter_mut()
                            .zip(&outfmts)
                            .zip(&side_sinks)
                            .map(|((output, &outfmt), side)| {
                                crate::api::local_blast::FormatOutput {
                                    outfmt,
                                    // The block search runs only where no format's check can
                                    // refuse (`x_query_block_plan`): a refusal fails the call.
                                    refusal: crate::api::local_blast::FormatRefusal::FailsCall,
                                    sink: match side {
                                        Some(side) => OutputSink::Files(output, side.clone()),
                                        None => OutputSink::Writer(output),
                                    },
                                }
                            })
                            .collect(),
                        diagnostics: &mut **diagnostics,
                        later_warnings: later_warnings.as_deref_mut(),
                        html_posts: html_posts
                            .as_deref_mut()
                            .map(|posts| posts as &mut (dyn std::io::Write + Send)),
                        hits: x_want_records.then_some(
                            &mut x_collect as &mut crate::api::local_blast::HitsSink<'_>,
                        ),
                        observer: x_offset
                            .as_mut()
                            .map(|offset| offset as &mut (dyn FormatObserver + Send)),
                    };
                    let x_later_blocks = blocks.len() - index - 1;
                    x_ahead.run_block(x_next_queries, x_later_blocks, |x_block| {
                        run_resolved_in_pool(
                            args.clone(),
                            block_queries,
                            subject_source,
                            "",
                            subject_label,
                            &mut block_outputs,
                            pool,
                            None,
                            &BlastpRanges {
                                input_records: &query_records[block.clone()],
                                query_input: &block_input,
                                subject_placements: &subjects.placements,
                                x_ordinal_base: block.start,
                                searched_base: block.start,
                                split_search_spaces: None,
                            },
                            &BlastpReports {
                                reader_warnings: &reader_warnings[block.clone()],
                                identity_warned,
                                epilog: epilog && index + 1 == blocks.len(),
                                x_prolog: index == 0,
                                x_queries_before: block.start,
                                x_xml2_carry: x_kinds
                                    .contains(&XBlockReport::Xml2)
                                    .then_some(&x_xml2_carry),
                                last_call: index + 1 == blocks.len(),
                                psi_round: None,
                                stopped_pssms: None,
                                query_scopes: batches
                                    .query_scopes
                                    .as_ref()
                                    .map(|scopes| scopes.numbered_from(block.start)),
                                prefetch_scope: None,
                                strategy_next: None,
                                call_requests: None,
                                x_plan_blocks: false,
                            },
                            x_block,
                        )
                    })?;
                    drop(block_outputs);
                    if let Some(records) = x_records.as_mut() {
                        for found in x_found {
                            // The block has ended, its reports hold no row: the rows move here.
                            let found = std::sync::Arc::try_unwrap(found)
                                .unwrap_or_else(|shared| shared.as_ref().clone());
                            records.extend(found.into_iter().map(|mut record| {
                                record.hit.q_idx += x_first_query;
                                record
                            }));
                        }
                    }
                    x_hsps_before += x_records.as_ref().map_or(0, Vec::len) - records_before;
                }
                for sink in &mut sinks {
                    sink.flush()?;
                }
                Ok(true)
            })?;
        if searched_in_blocks {
            // EXPERIMENT (LOSAT_X_BPQBLOCK): the hit records of all queries, once.
            if let (Some(sink), Some(records)) = (hits.as_mut(), x_records.take()) {
                sink(std::sync::Arc::new(records));
            }
            return Ok(());
        }
    }
    // The search strategy's request of a batch outside every call (`-export_search_strategy`): the
    // batch that stops or fails the run when no batch before it was searched, written before any
    // output of the run (`CallRequests` releases the others); NCBI's posts of the bytes replaced
    // in it follow the batch's reader messages, before its error (`write_request_posts`).
    let write_alone = |span: std::ops::Range<usize>| -> Result<()> {
        match args.strategy.as_ref() {
            Some(strategy) => super::archive::requests(
                subject_source
                    .database()
                    .is_none()
                    .then(|| subject_source.fasta_records()),
                args.lcase_masking,
                args.psi_search.as_ref().and_then(|psi| psi.request_pssm()),
            )
            .write(
                strategy,
                span.start,
                &searched_queries[span.clone()],
                query_input.placements.get(span.start),
            ),
            None => Ok(()),
        }
    };
    // A batch whose request follows the reports of the batches searched before it (`NextRequest`):
    // the batch that stops the run (every query without letters), or the batch that fails with a
    // header file NCBI cannot read; NCBI writes its request once those are formatted.
    let next_request = |span: std::ops::Range<usize>| NextRequest {
        first_query: span.start,
        records: &searched_queries[span.clone()],
        placement: query_input.placements.get(span.start),
    };
    if batches.searched_done > 0 {
        crate::utils::threading::with_search_pool(args.num_threads, "blastp", |pool| {
            let result = search_batches(
                args.clone(),
                outputs,
                pool,
                batches.input_done,
                batches.searched_done,
                epilog,
                batches.stop_batch.clone().map(next_request),
            );
            // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
            // A header file that NCBI cannot read (or rpsblast's profile sequences, read by a
            // batch's traceback, or blastp's `-use_sw_tback` redo whose batch's last redone match
            // fails, `Traceback failed`, `redo_run::sw_batch_failed`) fails the batch whose
            // hit lists reach it first (`database::HeaderReadError`), after the outfmt 0 prolog,
            // the reports of the batches before (without the epilog) and that batch's reader
            // messages; so do LOSAT's refusals of a step where NCBI's search of the batch stops or
            // never ends or crashes (`redo_run::sw_endless`, `redo_run::sw_batch_failed` where the
            // batch's first query holds a redone match, a query chunk split again or without a
            // query in `search_split_query_batches`), as TBLASTX's and TBLASTN's
            // (`blastinput::app::ErrorAfterBatchReports`). LOSAT searches all batches together: the
            // batches before are searched again alone. A failure that the first search meets
            // before its redo (unified P-values refused for a batch,
            // `redo_run::evalue_adjust_contexts`) may come after a batch whose failure only the
            // redo finds (`-use_sw_tback`): the search again then fails at that earlier batch,
            // which NCBI reached first, and the batches before it are searched again in turn
            // (review R3 of line RG item X5, M1).
            let as_failure = |result: Result<()>| match result {
                Err(error) if error.is::<database::HeaderReadError>() => Ok(error
                    .downcast::<database::HeaderReadError>()
                    .expect("the error was checked to be a HeaderReadError")),
                other => Err(other),
            };
            let mut failure = match as_failure(result) {
                Ok(failure) => failure,
                Err(other) => return other,
            };
            let (searched_start, failed_end) = loop {
                let (input_start, searched_start) = batches.start_of(failure.query);
                // The failed batch's search strategy request, which NCBI writes before its
                // search: after the requests of the batches before, which their search again
                // releases with it (`NextRequest`), or alone, before any output.
                let span = batches.searched_span(failure.query);
                let failed_end = span.end;
                if searched_start == 0 {
                    write_alone(span)?;
                    write_pairwise_prologs(&args, subject_source, subject_label, outputs)?;
                    break (searched_start, failed_end);
                }
                // The first search looked the taxonomy up for the batches before (their
                // prefetch and reports): the search again writes those texts again, from the
                // state of a run that has looked nothing up (review 1 of line DBID mixdb2, H1).
                if let Some(database) = database {
                    database.restart_taxonomy();
                }
                match as_failure(search_batches(
                    args.clone(),
                    outputs,
                    pool,
                    input_start,
                    searched_start,
                    false,
                    Some(next_request(span)),
                )) {
                    Ok(earlier) => failure = earlier,
                    Err(Ok(())) => break (searched_start, failed_end),
                    Err(other) => return other,
                }
            };
            outputs
                .diagnostics
                .write_all(&batches.warnings[searched_start])?;
            write_request_posts(&args, &mut *outputs.diagnostics)?;
            // NCBI's warning of its first read of a UTF8String header, where the failed batch
            // made it (before its failure; a batch before posted its own in its reports).
            if let Some(posts) = database.and_then(|database: &database::BlastpDatabase| {
                database::Utf8Read::post_in(
                    failure.utf8.as_ref(),
                    0..failed_end,
                    database.seqdb(),
                    args.app.ncbi_name(),
                )
            }) {
                posts.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
            }
            Err(failure.error)
        })?;
    } else if batches.stop.is_some() {
        // NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253
        // The first batch fails after the outfmt 0 prolog. A batch whose set-up fails (every query
        // without letters) has its search strategy request: NCBI made its query factory before
        // the search (`QueryBatches::stop_batch`); it is written before the prolog.
        if let Some(span) = batches.stop_batch.clone() {
            // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:506, 340
            // NCBI's stream that throws while it writes the request (an empty descriptor of a
            // PSSM query, `encoder::NullElement`) ends the run after the prolog (NCBI BLAST+
            // 2.17.0, probe D3c-3b `probe5` `t_dn1.strat0`).
            if let Err(error) = write_alone(span) {
                if error.is::<crate::cli::NativeError>() {
                    write_pairwise_prologs(&args, subject_source, subject_label, outputs)?;
                }
                return Err(error);
            }
        }
        write_pairwise_prologs(&args, subject_source, subject_label, outputs)?;
    } else {
        // A reader at the end of its input reads no batch (its `End()` is true): the reports of
        // no query (an empty pipe, or `-query -` after `-subject -`).
        write_reports_without_queries(&args, subject_source, subject_label, outputs)?;
    }
    // The stopping batch's reader messages, NCBI's posts of the bytes replaced in its request
    // (released after the reports of the batches before, or written alone), its error.
    if let Some((warnings, error)) = batches.stop {
        outputs.diagnostics.write_all(&warnings)?;
        write_request_posts(&args, &mut *outputs.diagnostics)?;
        return Err(error);
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1426
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the encoding of `run_call`
/// (`prepare_blastp_subjects_preserving_order`), made once for all blocks.
fn x_bp_encode_block_subjects(
    subject_records: &[FastaRecord],
    requested_threads: usize,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Result<super::x_query_blocks::XBlockSubjects> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let prepared = prepare_blastp_subjects_preserving_order(
        subject_records,
        requested_threads,
        pool.threads(),
        pool,
    )?;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let prepared = {
        let _ = pool;
        prepare_blastp_subjects_preserving_order(subject_records, requested_threads)?
    };
    Ok(super::x_query_blocks::XBlockSubjects::new(
        subject_records,
        prepared
            .into_iter()
            .map(|prepared| prepared.encoded)
            .collect(),
    ))
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the blocks of the block search, or `None` for the search of all
// queries. The block search is taken only where the blocks' reports written one after the other
// are the report of the search of all queries: tabular format 6, outfmt 0, 7, 8, 9, 12, 5 or 13 to
// 16 (`x_query_blocks.rs`; the kind of report of each output is returned with the blocks), also
// several outputs, the hit records and an observer (`run_local`, ABI v2 of the Web build) where
// every output is a writer of the caller, no `-db` is searched and an observer comes with the hit
// records (an outfmt 13 to 16 output is not taken in such a call: see the arm below), no
// `-query_loc` (every record is searched whole), every batch read, more than one block (one block
// is the search of all queries), no subject whose tabular ID the reference path rejects before it
// writes any row (`tabular_subject_ids`, for outfmt 6 and 7; with `-db`, the definition lines are
// checked where the blocks are searched, `database::x_titles_allow_blocks` (outfmt 8, 9, 12 and 5:
// `x_titles_read_allow_blocks`); outfmt 0's titles also there, `x_block_shown_titles_allowed`;
// outfmt 5's names too, `x_query_blocks::x_block_xml_subjects_allowed` and
// `x_block_xml_run_allowed`; each check where any output is of its kind),
// no query whose effective search space in NCBI's `Int8` the reference path rejects
// (`blastp_database_search_space`, with `-db` or `-dbsize`), no database whose setup count the
// reference path refuses (`call_setup::setup_count_refused`), no query that is also a subject under
// `-parse_deflines` where outfmt 0 can fail alone (its refusal of such a record's masks is outfmt
// 0's alone, `believed::refuse_shared_masks_in`), and none of the other refusals of the
// search of all queries that come before its rows (`x_block_search_refuses_nothing`): a block
// search would have written the reports of the blocks before. Both lookup tables
// (`AaLookupTable`, and `CompressedAaLookupTable` of blastp-fast) take it.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
#[allow(clippy::too_many_arguments)]
fn x_query_block_plan(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batches: &crate::blastinput::query_batch::QueryBatches,
    subject_source: BlastpSubjectSource<'_>,
    outputs: &ReportOutputs<'_>,
    batch_size: u32,
) -> Option<(
    Vec<super::x_query_blocks::XBlockReport>,
    Vec<std::ops::Range<usize>>,
)> {
    use super::x_query_blocks::XBlockReport;
    super::x_query_blocks::x_bp_qblock()?;
    // psiblast's PSSMs (`psiblast_ascii_pssms`) are written by the search of all queries. So is a
    // `-use_sw_tback` search: a later batch's Smith-Waterman redo can fail NCBI's traceback
    // (`redo_run::sw_batch_failed`; NCBI's crash there is LOSAT's refusal) or meet the redo that
    // NCBI never ends (LOSAT's refusal, `redo_run::sw_endless`), which come after the reports of
    // the batches before in the search of all queries (`search`). So is a search with unified
    // P-values without a search space given (`-searchsp`): the redo's set-up may refuse a query
    // batch (`redo_run::evalue_adjust_contexts`: a query before one without a search space, known
    // only once the queries are set up), after the reports of the batches before
    // (`database::HeaderReadError`). An ungapped search plans no blocks (`args.ungapped`).
    if args.query_range.is_some()
        || args.psi_pssm.is_some()
        || args.psi_search.is_some()
        || args.rps_search.is_some()
        || args.ungapped
        || args.use_sw_tback
        || (args.comp_based_stats.unified_p && args.searchsp.unwrap_or(0) == 0)
        || outputs.formats.is_empty()
        || batches.stop.is_some()
        || batches.searched_done != queries.len()
        || batches.input_done != queries.len()
    {
        return None;
    }
    // Several outputs, the hit records or an observer. Each
    // output of several is a writer of the caller: the search of all queries writes each output's
    // whole report in turn and opens its sink just before (`OutputSink::open`), while the block
    // search opens every sink once and writes them block by block, so two outputs on standard
    // output would get each block's formats in turn, two outputs to one file would write over each
    // other through two handles (the search of all queries leaves the last report), and a file that
    // cannot be created would stop the run after the first block's part of the outputs before it
    // instead of after their whole reports (test
    // `several_file_outputs_are_written_one_after_the_other`); no `-db` (the records' `s_idx` of a
    // database search rank the subjects that the call's lists name); an observer only with the hit
    // records, which count the HSPs of each block for the observer's indices (`XOffsetObserver`).
    let x_web = outputs.formats.len() > 1 || outputs.hits.is_some() || outputs.observer.is_some();
    if x_web
        && (subject_source.database().is_some()
            || (outputs.observer.is_some() && outputs.hits.is_none())
            || (outputs.formats.len() > 1
                && !outputs.formats.iter().all(|format| format.sink.is_writer())))
    {
        return None;
    }
    // A query that is also a subject under `-parse_deflines`, in a call whose outfmt 0 can fail
    // alone (LOSAT Web's outfmt 0 when the argv asks for another format): outfmt 0's refusal of
    // such a record's masks is that format's alone (`believed::refuse_shared_masks_in`), which the
    // blocks' outputs, whose refusals fail the call, cannot make (line WQA's sweep, D7, 2026-11-19).
    if args.parse_deflines
        && outputs.formats.iter().any(|format| {
            matches!(
                format.refusal,
                crate::api::local_blast::FormatRefusal::FailsAlone
            ) && matches!(
                blastp_text_format(format.outfmt),
                Ok(OutputFormat::Pairwise)
            )
        })
        && queries
            .iter()
            .any(|record| record.believed().is_some_and(|ids| ids.shared.is_some()))
    {
        return None;
    }
    // Outfmt 1 to 4, 8, 9, 12, 5, 15 and 16 are the ones `blastp_text_format` does not know.
    // The custom tabular reports of outfmt 6 and 7 (a field list or a delimiter, `tabular`) as
    // the default ones: each query's lines depend only on its final list and fixed inputs (the
    // subjects' ids, lengths and titles, a `-subject` search's taxonomy database warning), outfmt
    // 7's last line counts the queries before (`BlastpReports::x_queries_before`). With `-db`, a
    // field list that reads more of the definition lines than the default report (`sall*`) or
    // the taxonomy (`staxid`, `staxids`, the names) can be refused after the rows of some blocks
    // (`tabular::write_blastp_tabular`, `ListedTitles::tabular_taxonomy`,
    // `BlastpDatabase::tabular_names`): those take the reference path, as outfmt 10 and 20 (whose
    // header line is the run's).
    let custom_kind = |spec: &crate::report::tabular::TabularSpec| {
        let needs = spec.needs();
        let db_reads_more =
            subject_source.database().is_some() && (needs.all_deflines || needs.taxonomy);
        (!db_reads_more).then_some(spec.kind)
    };
    let default_kind = |outfmt: &str| match blastp_text_format(outfmt) {
        Ok(OutputFormat::Tabular | OutputFormat::TabularWithComments) => {
            Some(XBlockReport::Tabular)
        }
        Ok(OutputFormat::Pairwise) => Some(XBlockReport::Pairwise),
        // Outfmt 1 to 4 (BLASTP's, `blastp_pairwise_view_of`): outfmt 0's report with each
        // query's alignments as one section, which reads only that query's hits and the
        // sequences (`report::query_anchored`).
        Err(_)
            if crate::blastinput::app::parse_formatting_string(outfmt)
                .is_ok_and(|choice| blastp_pairwise_view_of(args.app, choice.number).is_some()) =>
        {
            Some(XBlockReport::Pairwise)
        }
        Err(_) if blastp_seqalign_format(outfmt).is_some() => Some(XBlockReport::Seqalign),
        // Outfmt 11 (the BLAST archive: one document per NCBI batch, `x_query_blocks.rs`). Only in a
        // call of one output without hit records and observer, as outfmt 15 and 16: the run-level
        // `ArchiveHandle` and the batches' messages were argued for the CLI's shape (Owner-delegated
        // 2026-10-23, DECISIONS 44).
        Err(_) if blastp_is_archive(outfmt) => (!x_web).then_some(XBlockReport::Seqalign),
        Err(_) if blastp_is_xml(outfmt) => Some(XBlockReport::Xml),
        // Outfmt 13 to 16 (`x_xml2_kind_takes_blocks`; 13 and 14 write side files). Only in a
        // call of one output without hit records and observer (the CLI's shape): the run-level
        // frame and carry of the XML2 document (`Xml2Carry`) were checked for that call only
        // (Owner-delegated 2026-10-22).
        Err(_) if blastp_xml2_kind(outfmt).is_some_and(x_xml2_kind_takes_blocks) => {
            (!x_web).then_some(XBlockReport::Xml2)
        }
        _ => None,
    };
    // A custom tabular report of outfmt 6 or 7 takes the default one's kind.
    let kind_of = |outfmt: &str| match super::tabular::custom_tabular_spec(outfmt, args.app) {
        Some(spec) => match custom_kind(spec) {
            Some(crate::report::tabular::TabularKind::Tabular) => default_kind("6"),
            Some(crate::report::tabular::TabularKind::Commented) => default_kind("7"),
            _ => None,
        },
        None => default_kind(outfmt),
    };
    let kinds = outputs
        .formats
        .iter()
        .map(|format| kind_of(format.outfmt))
        .collect::<Option<Vec<_>>>()?;
    // Outfmt 1 to 4 refuse a section that NCBI BLAST+ aborts on, and a database sequence that NCBI
    // shows as a query, in the call that reports the query, before any byte of it (`run_call`, as
    // outfmt 0 refuses the masks of a shared record, `believed::refuse_shared_masks_in`): in a block
    // search, after the bytes of the blocks before (NCBI's own bytes before its abort depend on
    // its buffers, so no prefix is NCBI's); the blocks are taken as for outfmt 0. (The second
    // has parsed Seq-ids, which the block search never takes, `database::x_titles_allow_blocks`:
    // it refuses in the search of all queries, before any byte.)
    // A subject set whose tabular id lookup fails (`check_failed_lookups`) or, for outfmt 0 to 4,
    // whose title of punctuation `check_shown_subject_titles` refuses takes the search of all
    // queries, which refuses before any byte: there LOSAT Web's format that the argv does not ask
    // for fails alone (`FormatRefusal`), and the block search's outputs fail the call.
    let tabular = kinds.contains(&XBlockReport::Tabular);
    let pairwise = kinds.contains(&XBlockReport::Pairwise);
    match subject_source {
        BlastpSubjectSource::Fasta(subjects) => {
            if tabular
                && subjects
                    .iter()
                    .any(|subject| !matches!(subject.subject_tabular_id(true), Ok(Some(_))))
            {
                return None;
            }
            if pairwise && check_shown_subject_titles(0..subjects.len(), subjects).is_err() {
                return None;
            }
        }
        // A database without sequences fails in the first batch (`run_call`). A setup count that
        // LOSAT refuses (`call_setup::setup_count_refused`, with the setup count of a search
        // without RPS, `db_extent`'s: the statistics' count) is refused by the call that holds a
        // valid query: the call of all queries refuses before any output, a block only after the
        // reports of the blocks before, whose queries were invalid (review D of line RG D1a-4,
        // SF1).
        BlastpSubjectSource::Database(database) => {
            if database.num_sequences() == 0
                || call_setup::setup_count_refused(args, {
                    let statistics = database.statistics();
                    (statistics.setup_total, statistics.setup_count)
                })
            {
                return None;
            }
        }
    }
    let lengths: Vec<usize> = queries.iter().map(|record| record.seq().len()).collect();
    // The effective search space is at most `(db + |n| q) q` for a database of `db` letters and
    // `n` sequences (NCBI's `Int4`) and a query of `q` letters: the length adjustment `a` is
    // between 0 and `q`, so the effective database length `db - n a` is at most `db` for `n >= 0`
    // and `db + |n| q` for a negative `n` (an alias file's `NSEQ`), and the effective query length
    // at most `q`. No query's product can pass `Int8` when that bound of the longest one does not.
    if subject_source.database().is_some() || args.dbsize > 0 {
        let (db_length, num_seqs) = match subject_source.database() {
            Some(database) if args.dbsize <= 0 => (
                i128::from(database.total_letters()),
                database.statistics_num_sequences(),
            ),
            Some(database) => (i128::from(args.dbsize), database.statistics_num_sequences()),
            None => (i128::from(args.dbsize.max(0)), 1),
        };
        let longest = lengths.iter().copied().max().unwrap_or(0) as i128;
        let negative_seqs = i128::from(num_seqs.min(0).unsigned_abs());
        if (db_length + negative_seqs * longest) * longest > i64::MAX as i128 {
            return None;
        }
    }
    let query_batches = query_input.batches(batch_size);
    let blocks = x_batch_blocks(args, queries, &query_batches)?;
    x_block_search_refuses_nothing(args, queries, &query_batches, &lengths)
        .then_some((kinds, blocks))
}

// No NCBI counterpart: chooses LOSAT's blocks; it does not change any value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the blocks of whole batches of `queries` (the searched queries of
/// `query_batches`, which are batches of the whole input) closed on their letters
/// (`x_query_blocks::x_query_blocks` with `x_bp_qblock`), where there is more than one, or one
/// block of letters whose standard table's words pass the most of a block (`x_aa_one_block_splits`;
/// the blocks are then closed on their words on the search's pool, `x_word_blocks`); `None`
/// otherwise (the switch off, a batch the blocks cannot take). The blocks of the block search
/// (`x_query_block_plan`, with the checks of its reports) and of the preliminary stage of the
/// search of all queries (`x_prelim_block_plan`).
fn x_batch_blocks(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_batches: &[seq_range::RangedBatch],
) -> Option<Vec<std::ops::Range<usize>>> {
    let block_letters = super::x_query_blocks::x_bp_qblock()?;
    let lengths: Vec<usize> = queries.iter().map(|record| record.seq().len()).collect();
    let blocks =
        super::x_query_blocks::x_query_blocks(query_batches, &lengths, block_letters, None)?;
    (blocks.len() > 1 || x_aa_one_block_splits(args, queries, query_batches)).then_some(blocks)
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:215-219
// EXPERIMENT (LOSAT_X_BPQBLOCK): the blocks of the preliminary stage of the search of all queries
// (`run_call`'s `x_blocks`): the blocks of whole batches of `x_batch_blocks`, closed on their words
// on the search's pool (`x_word_blocks`, which leaves a compressed table that could fill NCBI's
// overflow banks to the one table of all queries, where the refusal is made), where there is more
// than one. Any report: the blocks write nothing, and the hit lists of all queries are joined
// before the redo and the reports, so the reports, their order and every error after the set-up
// are those of one table. The set-up's errors and refusals come before the first block: the call
// of all queries makes them, but it builds no table, so a refusal of a table's build (a compressed
// table that fills NCBI's overflow banks) cannot come first in a block only because
// `x_word_blocks` plans no blocks wherever a batch's table could fill them. Only blastp: psiblast
// searches one query per batch (its `-in_pssm`, `-sum_stats` and rounds), and rpsblast scans with
// the database's table.
// No NCBI counterpart: chooses LOSAT's blocks; it does not change any value NCBI computes.
fn x_prelim_block_plan(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batch_size: u32,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<Vec<std::ops::Range<usize>>> {
    // An ungapped search plans no blocks (`args.ungapped`: its queries are not split either,
    // `search_split_query_batches`).
    if args.app != BlastpApp::Blastp
        || args.rps_search.is_some()
        || args.psi_search.is_some()
        || args.sum_stats
        || args.ungapped
    {
        return None;
    }
    let query_batches = query_input.batches(batch_size);
    let blocks = x_batch_blocks(args, queries, &query_batches)?;
    let (blocks, _) = x_word_blocks(args, queries, query_input, batch_size, &blocks, pool)?;
    (blocks.len() > 1).then_some(blocks)
}

// No NCBI counterpart: the input records of a block of whole batches; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the input records of `block`, whole batches of `query_batches`
/// (`x_prelim_block_plan`; the batches' searched queries follow each other, none empty): from the
/// input start of the batch whose searched queries start the block to the input end of the batch
/// whose searched queries end it.
fn x_block_input(
    query_batches: &[seq_range::RangedBatch],
    block: &std::ops::Range<usize>,
) -> std::ops::Range<usize> {
    let first =
        &query_batches[query_batches.partition_point(|batch| batch.searched.start < block.start)];
    let last =
        &query_batches[query_batches.partition_point(|batch| batch.searched.end < block.end)];
    assert!(
        first.searched.start == block.start && last.searched.end == block.end,
        "a block is whole batches of the call"
    );
    first.input.start..last.input.end
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): whether the search of all queries passes the checks that refuse
// it before any row is written and that a block search would make only in the block that holds
// the cause, after the rows of the blocks before:
// - the diagonal table of all queries (`check_diag_table_window` in `call_setup::scan_shape`, with
//   the concatenated length of all queries; a block's is shorter);
// - the split of a long query batch (the refusals of `search_split_query_batches`);
// - with the compressed lookup table, a table of one of NCBI's query batches (or of a query
//   chunk of a split batch) that fills NCBI's overflow banks: `x_word_blocks`, which the block
//   search tests on its pool before the first block (and leaves the search to the reference path
//   when it fails).
// A search that passes them is refused by none of them in any block either: a block's batches are
// those of the whole input (`x_query_blocks.rs`).
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
fn x_block_search_refuses_nothing(
    args: &ResolvedBlastpArgs,
    _queries: &[FastaRecord],
    query_batches: &[seq_range::RangedBatch],
    lengths: &[usize],
) -> bool {
    use crate::algorithm::blastp::protein_query_split::{
        protein_chunk_would_be_split, split_protein_batch,
    };
    // `blastp_concatenated_query_length`: each query and the sentinel after it, but the last.
    let query_length = lengths
        .iter()
        .map(|&length| length as u64 + 1)
        .sum::<u64>()
        .saturating_sub(1);
    let Ok(query_length) = i32::try_from(query_length) else {
        return false;
    };
    let window = args.window_size as i32;
    if crate::blastinput::app::check_diag_table_window(query_length, window, "BLASTP").is_err() {
        return false;
    }
    let Ok(sizes) = protein_query_split_sizes(10000) else {
        return false;
    };
    let searched_batches = || {
        query_batches
            .iter()
            .map(|batch| batch.searched.clone())
            .filter(|searched| !searched.is_empty())
    };
    for searched in searched_batches() {
        if let Some(chunks) = split_protein_batch(&lengths[searched], sizes) {
            if sizes.negative_chunk_size()
                || protein_chunk_would_be_split(&chunks, sizes)
                || chunks.iter().any(|chunk| chunk.queries.is_empty())
            {
                return false;
            }
        }
    }
    true
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_aux_priv.cpp:206-215
// NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:208-211
// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_aalookup.h:189,221-231
// EXPERIMENT (LOSAT_X_BPQBLOCK): the last check of `x_block_search_refuses_nothing`, made on the
// search's pool before the first block, and the blocks' ends. With the compressed lookup table:
// - whether no table of one of NCBI's query batches (or of a query chunk of a split batch) can fill
//   NCBI's overflow banks (`try_build_blosum62_compressed_lookup`). A table takes at most a quarter
//   of its words as overflow cells, and a batch's tables have no more words than the batch's
//   offsets add (`CompressedWordCount`): a batch whose words are at most four times the banks'
//   cells cannot fill them; a batch with more leaves the search to the reference path (`None`),
//   which tests the batch exactly.
// - the blocks closed also on their words (`x_query_blocks`). NCBI builds the table of one batch,
//   searches the batch and frees the table before the next batch (above), so its largest table is
//   that of the batch with the most words: a backbone cell has slots for four offsets whatever the
//   threshold (24 bytes with its count), and the words past them go to overflow cells of four (24
//   bytes). A block's table (`x_clut_direct.rs`: an entry of 4 bytes per word and per backbone cell)
//   is kept to the slots of that table: at most four words per backbone cell plus the words of the
//   batch with the most, so a block's finished table takes about the memory of NCBI's largest
//   table, not a multiple of it (the builder's own vectors, and at more than one thread the next
//   block's table built ahead, `x_lookup_ahead.rs`, come on top). At the default thresholds of
//   blastp-fast (word size 5, 20) and of word size 6 (21) a block of `LOSAT_X_BPQBLOCK` letters
//   holds fewer words (the backbone's slots dominate), and the blocks are those of the letters; at
//   those of word sizes 5 (19.3) and 7 (20.25) of `-task blastp` the words of a block of 50,000
//   letters come near the most, and a block can end one batch earlier; with a lowered threshold a
//   block holds one batch or a few.
// The batches are counted one by one, on the pool's threads when it has more than one, with one
// counter whose kept counts serve them all. Each block comes with its words (the sum of its
// batches'), the most appends of its table's build (`x_query_blocks::x_block_appends`).
// The standard table (`AaLookupTable`) has no overflow banks; its blocks are closed on their words
// too (`x_aa_word_blocks`, item SC8b), and come with no words (its builds reserve nothing).
// No NCBI counterpart: chooses LOSAT's path and blocks; it does not change any value NCBI computes.
fn x_word_blocks(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batch_size: u32,
    blocks: &[std::ops::Range<usize>],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<(Vec<std::ops::Range<usize>>, Vec<u64>)> {
    use crate::algorithm::tblastx::lookup::compressed::{
        CompressedWordCount, COMPRESSED_HITS_PER_BACKBONE_CELL,
    };
    if args.lookup_table_type != BlastpLookupTableType::CompressedAaLookupTable {
        return Some((
            x_aa_word_blocks(args, queries, query_input, batch_size, pool)
                .unwrap_or_else(|| blocks.to_vec()),
            Vec::new(),
        ));
    }
    // The threshold as the table is built (`run_call`), with the words of the search's
    // matrix (`blastp_matrix_compressed_lookup` builds another matrix's table).
    let counter = CompressedWordCount::new_of(
        args.scoring.matrix,
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(100.0 * args.threshold),
    )?;
    // The words of all letters of the queries: the lookup query
    // (`encode_protein_query_frame_with_seg`) has the same letters, and its lookup segments
    // (`build_lookup_segments_for_context`) leave out those that SEG or `-lcase_masking` mask, so
    // a block's table holds at most these words (all of them without masks). Counting the masked
    // letters too only closes such a block earlier and reserves room that its table does not fill
    // (`x_block_appends`) and never touches.
    let words = |searched: &std::ops::Range<usize>| {
        let encoded: Vec<EncodedProtein> = queries[searched.clone()]
            .iter()
            .map(|record| encode_protein_sequence(record.seq()))
            .collect();
        let mut letters = encoded
            .iter()
            .map(|query| &query.aa_seq[1..query.aa_seq.len() - 1]);
        counter.batch_words(&mut letters)
    };
    let query_batches = query_input.batches(batch_size);
    let batch_words: Vec<u64> = x_each_batch_on(&query_batches, pool, &words)
        .into_iter()
        .collect::<Option<_>>()?;
    let most = (COMPRESSED_HITS_PER_BACKBONE_CELL as u64)
        .saturating_mul(counter.backbone_cells())
        .saturating_add(batch_words.iter().copied().max().unwrap_or(0));
    let lengths: Vec<usize> = queries.iter().map(|record| record.seq().len()).collect();
    let blocks = super::x_query_blocks::x_query_blocks(
        &query_batches,
        &lengths,
        super::x_query_blocks::x_bp_qblock()?,
        Some((&batch_words, most)),
    )?;
    // A block is whole batches, in order (`x_query_blocks`).
    let mut block_words = vec![0u64; blocks.len()];
    let mut block = 0;
    for (batch, &words) in query_batches.iter().zip(&batch_words) {
        while blocks[block].end < batch.searched.end {
            block += 1;
        }
        block_words[block] = block_words[block].saturating_add(words);
    }
    Some((blocks, block_words))
}

// No NCBI counterpart: where the counts of `x_word_blocks` run; it changes no value NCBI computes.
/// `each` of every batch's searched queries, in batch order: on the pool's threads when it has
/// more than one (with one counter, whose kept counts serve all batches), else on this thread.
fn x_each_batch_on<T: Send>(
    query_batches: &[seq_range::RangedBatch],
    pool: &crate::utils::threading::SearchPool<'_>,
    each: &(dyn Fn(&std::ops::Range<usize>) -> T + Sync),
) -> Vec<T> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        return pool.install(|| {
            query_batches
                .par_iter()
                .map(|batch| each(&batch.searched))
                .collect()
        });
    }
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let _ = pool;
    query_batches
        .iter()
        .map(|batch| each(&batch.searched))
        .collect()
}

// No NCBI counterpart: a floor of the most words of a block of the standard table
// (`x_aa_word_blocks`); it changes no value NCBI computes.
/// The words a block of the standard table may always hold: 2^22 (16 MiB of `Int4` offsets).
const X_AA_BLOCK_WORDS_FLOOR: u64 = 1 << 22;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:283-358 (BlastAaLookupFinalize:
// a cell of at most AA_HITS_PER_CELL offsets keeps them in the backbone, a cell of more puts all of
// its offsets in the overflow array)
// EXPERIMENT (LOSAT_X_BPQBLOCK), item SC8b: the blocks of the standard table (`AaLookupTable`) closed
// also on their words, as the compressed table's (`x_word_blocks`). A query offset adds the cells
// T(s) of its word to the table (`XWideWordCount`: about 35 per letter for BLOSUM62's words of 3 at
// threshold 11, over a thousand for its words of 4), NCBI builds and frees the table of one batch at
// a time, and its largest table takes its backbone (cells of AA_HITS_PER_CELL offsets) and the
// overflow of the batch with the most words. A block's table is kept to that: at most
// AA_HITS_PER_CELL words per backbone cell plus the words of the batch with the most, or
// `X_AA_BLOCK_WORDS_FLOOR` where that is more, so that a block of `LOSAT_X_BPQBLOCK` letters whose
// table takes at most 16 MiB of offsets ends on its letters (BLOSUM62's words of 2 and 3 at their
// default thresholds: below 2 million words in 50,000 letters, 7-9 MiB); a table of 4-letter words
// takes one batch or a few. A plan of one block of letters comes here too where its words pass the
// most (`x_aa_one_block_splits`). `None` where the word size has no counter (none of the standard
// table's).
// No NCBI counterpart: chooses LOSAT's blocks; it does not change any value NCBI computes.
fn x_aa_word_blocks(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_input: &seq_range::QueryInput,
    batch_size: u32,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> Option<Vec<std::ops::Range<usize>>> {
    let counter = x_aa_word_counter(args)?;
    let query_batches = query_input.batches(batch_size);
    let batch_words = x_each_batch_on(&query_batches, pool, &|searched| {
        x_aa_batch_words(&counter, &queries[searched.clone()])
    });
    let lengths: Vec<usize> = queries.iter().map(|record| record.seq().len()).collect();
    super::x_query_blocks::x_query_blocks(
        &query_batches,
        &lengths,
        super::x_query_blocks::x_bp_qblock()?,
        Some((&batch_words, x_aa_most_words(&counter, &batch_words))),
    )
}

// No NCBI counterpart: the counts of `x_aa_word_blocks`; they change no value NCBI computes.
/// The words of `queries`' tables (`XWideWordCount::batch_words`): of all their letters, as the
/// compressed table's count in `x_word_blocks` (SEG and `-lcase_masking` only leave words out).
fn x_aa_batch_words(
    counter: &crate::algorithm::tblastx::lookup::x_pair_lut::x_pair_lut_wide::XWideWordCount,
    queries: &[FastaRecord],
) -> u64 {
    let encoded: Vec<EncodedProtein> = queries
        .iter()
        .map(|record| encode_protein_sequence(record.seq()))
        .collect();
    let mut letters = encoded
        .iter()
        .map(|query| &query.aa_seq[1..query.aa_seq.len() - 1]);
    counter.batch_words(&mut letters)
}

// No NCBI counterpart: the most words of a block (`x_aa_word_blocks`); it changes no value NCBI
// computes.
/// AA_HITS_PER_CELL words per backbone cell plus the words of the batch with the most (NCBI's
/// largest table), or `X_AA_BLOCK_WORDS_FLOOR` where that is more.
fn x_aa_most_words(
    counter: &crate::algorithm::tblastx::lookup::x_pair_lut::x_pair_lut_wide::XWideWordCount,
    batch_words: &[u64],
) -> u64 {
    (crate::algorithm::tblastx::lookup::AA_HITS_PER_CELL as u64)
        .saturating_mul(counter.backbone_cells())
        .saturating_add(batch_words.iter().copied().max().unwrap_or(0))
        .max(X_AA_BLOCK_WORDS_FLOOR)
}

// No NCBI counterpart: the counter of the standard table's blocks, kept for the process; it changes
// no value NCBI computes.
/// The counter of `args`' standard table (its matrix, word size and threshold as
/// `blastp_aa_lookup` builds the table), kept with its counts for the next call of the same
/// search (`x_aa_one_block_splits` counts in the plan, `x_aa_word_blocks` again on the pool).
fn x_aa_word_counter(
    args: &ResolvedBlastpArgs,
) -> Option<Arc<crate::algorithm::tblastx::lookup::x_pair_lut::x_pair_lut_wide::XWideWordCount>> {
    use crate::algorithm::tblastx::lookup::x_pair_lut::x_pair_lut_wide::XWideWordCount;
    type Kept = Option<((ScoringMatrix, usize, i32), Arc<XWideWordCount>)>;
    static KEPT: std::sync::Mutex<Kept> = std::sync::Mutex::new(None);
    let key = (
        args.scoring.matrix,
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
    );
    let mut kept = KEPT
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner);
    if let Some((kept_key, counter)) = kept.as_ref() {
        if *kept_key == key {
            return Some(Arc::clone(counter));
        }
    }
    let counter = Arc::new(XWideWordCount::new(key.0, key.1, key.2)?);
    *kept = Some((key, Arc::clone(&counter)));
    Some(counter)
}

// No NCBI counterpart: whether the block search's plan goes on with one block of letters; it
// changes no value NCBI computes.
/// Whether the queries' one block of letters (`x_query_block_plan`) would be split on its words
/// (`x_aa_word_blocks`): the standard table's words of all the queries pass the most. Not counted
/// where they cannot: a word adds at most every cell of the alphabet's words of its length (the
/// cells T(s) are distinct), so letters whose number times those cells is at most the floor
/// (`X_AA_BLOCK_WORDS_FLOOR`) never pass it.
fn x_aa_one_block_splits(
    args: &ResolvedBlastpArgs,
    queries: &[FastaRecord],
    query_batches: &[seq_range::RangedBatch],
) -> bool {
    if args.lookup_table_type != BlastpLookupTableType::AaLookupTable {
        return false;
    }
    let letters: u64 = queries.iter().map(|record| record.seq().len() as u64).sum();
    let cells = (crate::algorithm::tblastx::lookup::LOOKUP_ALPHABET_SIZE as u64)
        .saturating_pow(args.word_size as u32);
    if letters.saturating_mul(cells) <= X_AA_BLOCK_WORDS_FLOOR {
        return false;
    }
    let Some(counter) = x_aa_word_counter(args) else {
        return false;
    };
    let batch_words: Vec<u64> = query_batches
        .iter()
        .map(|batch| {
            queries
                .get(batch.searched.clone())
                .map_or(0, |queries| x_aa_batch_words(&counter, queries))
        })
        .collect();
    batch_words.iter().sum::<u64>() > x_aa_most_words(&counter, &batch_words)
}

// EXPERIMENT (LOSAT_X_BPQBLOCK): the outfmt 0 report of the search of all queries refuses a
// subject title that LOSAT does not write as NCBI where the report shows the subject
// (`check_shown_subject_titles`), before it writes anything; a block search would have written the reports of the blocks before. So with outfmt 0
// and `-subject` the blocks are searched only when no subject's title is refused (with `-db`:
// `database::x_titles_allow_blocks`). Tested on the pool before the first block, in groups of
// consecutive subjects on its threads when it has more than one; the outcome is whether all pass.
// No NCBI counterpart: chooses LOSAT's path; it does not change any value NCBI computes.
fn x_block_shown_titles_allowed(
    subject_records: &[FastaRecord],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> bool {
    let allowed = |subjects: std::ops::Range<usize>| {
        check_shown_subject_titles(subjects, subject_records).is_ok()
    };
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    if pool.enabled() && pool.threads() > 1 {
        let group = subject_records.len().div_ceil(pool.threads()).max(1);
        return pool.install(|| {
            (0..subject_records.len())
                .into_par_iter()
                .step_by(group)
                .all(|start| allowed(start..(start + group).min(subject_records.len())))
        });
    }
    let _ = pool;
    allowed(0..subject_records.len())
}

/// The `Database:` text of the pairwise reports of a local `-subject` search
/// (`blastp_scan_mode_definition` of `blastp_subject_tag`).
fn blastp_database_name(args: &ResolvedBlastpArgs, subject_label: &str) -> String {
    blastp_scan_mode_definition(&blastp_subject_tag(args, subject_label))
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:242-243
// NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:497-499
/// The subject tag of the formatter (`m_SubjectTag`): `subject_label`, or the `-subject`
/// argument; psiblast makes its formatter without one (the default `kEmptyStr`).
fn blastp_subject_tag(args: &ResolvedBlastpArgs, subject_label: &str) -> String {
    match args.app {
        // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:182-184
        BlastpApp::Psiblast | BlastpApp::Rpsblast => String::new(),
        BlastpApp::Blastp if subject_label.is_empty() => args.subject.display().to_string(),
        BlastpApp::Blastp => subject_label.to_string(),
    }
}

// NCBI reference (598d8ae6): c++/src/objtools/align_format/align_format_util.cpp:418-424
/// The definition of a subject set in the pairwise report (`FillScanModeBlastDbInfo`), without
/// its period, which the writers add (`ensure_trailing_period`).
fn blastp_scan_mode_definition(tag: &str) -> String {
    if tag.is_empty() {
        "User specified sequence set".to_string()
    } else {
        format!("User specified sequence set (Input: {tag})")
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:118-119
/// Writes and flushes the prolog of every outfmt 0 to 4, 15, 16 and 20 format (`PrintProlog`;
/// 1 to 4 where BLASTP writes them, `blastp_pairwise_view_of`), for a run that stops at its
/// first query batch (outfmt 20's line of the field names; its formatter's warning came before
/// the search, `search`).
fn write_pairwise_prologs(
    args: &ResolvedBlastpArgs,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let (num_sequences, total_letters) = subject_source.shown_totals();
    let pairwise_database = subject_source
        .database()
        .map(|database| database.pairwise_database_masked(args.db_mask_resolved.as_ref()));
    for format in outputs.formats.iter_mut() {
        let number = crate::blastinput::app::parse_formatting_string(format.outfmt)?.number;
        if blastp_pairwise_view_of(args.app, number).is_some() {
            let mut writer = format.sink.open()?;
            crate::report::pairwise::write_blastp_pairwise_prolog(
                &mut writer,
                crate::wording::REPORT_VERSION,
                args.app.ncbi_name(),
                blastp_database_name(args, subject_label).as_bytes(),
                num_sequences,
                total_letters,
                pairwise_database.as_ref(),
                args.html.as_deref().map(|links| {
                    let posts = crate::report::query_warnings::html_posts_sink(
                        outputs.html_posts.as_deref_mut(),
                        &mut *outputs.diagnostics,
                    );
                    (links, posts as &mut dyn Write)
                }),
            )?;
            writer.flush()?;
        } else if let Some(kind) = blastp_xml2_kind(format.outfmt) {
            // The XML2 and JSON documents start with their prolog too.
            let mut writer = format.sink.open()?;
            writer.write_all(&crate::report::blast_xml2::framing::prolog_bytes(kind))?;
            writer.flush()?;
        } else if let Some(spec) = super::tabular::custom_tabular_spec(format.outfmt, args.app)
            .filter(|spec| spec.kind == crate::report::tabular::TabularKind::CsvWithHeader)
        {
            let mut writer = format.sink.open()?;
            super::tabular::write_tabular_prolog(spec, &mut writer)?;
            writer.flush()?;
        }
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:253-259,294-295
// NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:1322-1325
/// The reports of a run without a query batch: the outfmt 0 prolog and epilog, the outfmt
/// 7 count of 0 queries, nothing for outfmt 6.
fn write_reports_without_queries(
    args: &ResolvedBlastpArgs,
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
) -> Result<()> {
    let (num_sequences, total_db_len) = subject_source.shown_totals();
    let pairwise_database = subject_source
        .database()
        .map(|database| database.pairwise_database_masked(args.db_mask_resolved.as_ref()));
    for format in outputs.formats.iter_mut() {
        // A Seq-annot format writes nothing without a query. BLAST XML writes the epilog's end
        // of line alone: `m_Outfile` gets `m_BlastXMLIncremental->m_SerialXmlEnd` and `endl`, with the
        // end of the document empty (blast_format.cpp:2242-2243).
        // The archive writes one document per batch: none without a query batch.
        if blastp_seqalign_format(format.outfmt).is_some() || blastp_is_archive(format.outfmt) {
            format.sink.open()?.flush()?;
            continue;
        }
        if blastp_is_xml(format.outfmt) {
            let mut writer = format.sink.open()?;
            writer.write_all(b"\n")?;
            writer.flush()?;
            continue;
        }
        // The custom tabular reports: outfmt 20's line of the field names, outfmt 7's count of 0
        // queries (`tabular::write_tabular_without_queries`).
        if let Some(spec) = super::tabular::custom_tabular_spec(format.outfmt, args.app) {
            let taxdb_found = spec.checks_taxdb()
                && match subject_source.database() {
                    Some(database) => database.taxdb_found(),
                    None => args.subject_taxdb_found,
                };
            let mut writer = format.sink.open()?;
            super::tabular::write_tabular_without_queries(
                spec,
                taxdb_found,
                args.app.ncbi_name(),
                &mut writer,
                &mut *outputs.diagnostics,
            )?;
            writer.flush()?;
            continue;
        }
        // The XML2 and JSON documents write the prolog and the epilog without a query
        // (`PrintProlog`, `PrintEpilog`: an empty pipe), and nothing for a query file that is
        // empty (`Query is Empty!`, before the formatter).
        if let Some(kind) = blastp_xml2_kind(format.outfmt) {
            let sink = format.sink.side_files().cloned();
            let mut writer = format.sink.open()?;
            let mut bytes = Vec::new();
            crate::report::blast_xml2::Xml2Document::for_run(
                kind,
                false,
                args.out.as_deref(),
                sink.as_ref(),
            )
            .epilog(&mut bytes)?;
            writer.write_all(&bytes)?;
            writer.flush()?;
            continue;
        }
        // Outfmt 1 to 4 write what outfmt 0 writes without a query.
        let view = blastp_query_anchored_view(args.app, format.outfmt);
        let outfmt = match view {
            Some(_) => OutputFormat::Pairwise,
            None => blastp_text_format(format.outfmt).map_err(anyhow::Error::msg)?,
        };
        let mut writer = format.sink.open()?;
        match outfmt {
            OutputFormat::Pairwise => {
                let report = BlastpPairwiseReport {
                    version: crate::wording::REPORT_VERSION.to_string(),
                    database_name: blastp_database_name(args, subject_label).into(),
                    database_num_sequences: num_sequences,
                    database_total_letters: total_db_len,
                    // NCBI reference (598d8ae6): c++/src/app/blast/rpsblast_app.cpp:194-229
                    // rpsblast's options name a matrix once a query batch's search has set it up
                    // (`CreateRpsStructures`): without a batch the epilog's matrix name is null
                    // (`write_matrix_line`).
                    matrix_name: args.rps_search.is_none().then(|| args.matrix_name.clone()),
                    gap_open: args.scoring.gap_open,
                    gap_extend: args.scoring.gap_extend,
                    word_threshold: args.threshold,
                    window_size: args.window_size as i32,
                    // None for an ungapped search, whose gap costs NCBI does not read
                    // (`call_setup::gapped_block`).
                    gapped_karlin: call_setup::gapped_block(args).unwrap_or_default(),
                    gumbel: if args.has_gumbel() {
                        Some(
                            lookup_protein_gumbel_params(&args.scoring, total_db_len).context(
                                "missing NCBI gumbel parameters for supported blastp scoring",
                            )?,
                        )
                    } else {
                        None
                    },
                    gapped: !args.ungapped,
                    show_linked_set_size: false,
                    num_descriptions: args.num_descriptions,
                    num_alignments: args.num_alignments,
                };
                // (A report without queries writes no alignment and no table: the line length
                // and the sorts are those of the run, `ResolvedBlastpArgs`.)
                let config = PairwiseConfig {
                    line_length: args.line_length,
                    show_frame: false,
                    program: args.app.ncbi_name().to_string(),
                    protein_matrix: args.scoring.matrix,
                    sort: args.sort,
                    view: view.unwrap_or_default(),
                    html: args
                        .html
                        .as_ref()
                        .map(|links| crate::report::html_ids::HtmlReport::new(links, [])),
                    unreadable_query: None,
                };
                // No query writes warnings; an `-html` prolog's references post on the
                // diagnostics (`pairwise::html_with_diagnostics`).
                let mut warnings = crate::report::query_warnings::QueryWarnings::of_posts(
                    &[],
                    None,
                    &mut *outputs.diagnostics,
                )
                .with_html_posts(
                    outputs
                        .html_posts
                        .as_deref_mut()
                        .map(|posts| posts as &mut dyn std::io::Write),
                );
                write_blastp_pairwise_report(
                    &[],
                    &mut writer,
                    &config,
                    &[],
                    &[],
                    &[],
                    &[],
                    &[],
                    &report,
                    None,
                    pairwise_database.as_ref(),
                    None,
                    Some(&mut warnings),
                    true,
                    true,
                    None,
                    &crate::report::db_subject::HitLoads::default(),
                )?;
            }
            OutputFormat::TabularWithComments => {
                writeln!(writer, "# BLAST processed 0 queries")?;
            }
            OutputFormat::Tabular => {}
        }
        writer.flush()?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:240-251
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_message.c:227-229
/// The error of a query batch's search that fails in its setup or preliminary stage (a database
/// without sequences, a sequence file that is missing or unreadable, RPS-BLAST's files, the
/// average subject length without the Gumbel block): NCBI has written the outfmt 0 prolog and
/// read the batch (its reader's messages) when the batch's search fails with `BLAST engine error:
/// <text>`, exit status 3. LOSAT searches all batches together, so the failure is the first
/// searched batch's: the call's first batch, or, when that batch holds invalid queries only,
/// which NCBI does not search (local_blast.cpp:177-180), a later one, after the reports of the
/// batches before it (`unsearched_batches_then_engine_error`). `first_batch` is an error of the
/// call's first batch whatever its queries' validity: a split query's chunk (its batch is
/// searched; its call has no search strategy requests, `search_split_query_batches` writes them),
/// or RPS-BLAST's files, which NCBI sets up before it checks the batch's queries
/// (`rps_hooks::RpsEngineError`). `database_name` is empty for a subject set.
#[allow(clippy::too_many_arguments)]
fn preliminary_search_error(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    pool: &crate::utils::threading::SearchPool<'_>,
    ranges: &BlastpRanges<'_>,
    reports: &BlastpReports<'_>,
    query_batches: &[seq_range::RangedBatch],
    query_batch_skipped: &[bool],
    first_batch: bool,
    database_name: &str,
    text: &str,
) -> anyhow::Error {
    if !(first_batch || query_batch_skipped.first() == Some(&false)) {
        return unsearched_batches_then_engine_error(
            args,
            query_records,
            subject_source,
            subject_label,
            outputs,
            pool,
            ranges,
            reports,
            query_batches,
            query_batch_skipped,
            database_name,
            text,
        );
    }
    let written = (|| -> Result<()> {
        // The failed batch's search strategy request: the call's first, whose search fails; it is
        // written before the prolog, and NCBI's posts of the bytes replaced in it after the
        // batch's reader messages.
        if let (Some(call), Some(batch)) = (
            reports.call_requests,
            query_batches
                .iter()
                .find(|batch| !batch.searched.is_empty()),
        ) {
            call.end_at(batch.searched.end)?;
        }
        // EXPERIMENT (LOSAT_X_BPQBLOCK): a block after the first of the block search follows
        // the prolog that the first block wrote.
        if reports.x_prolog {
            write_pairwise_prologs(args, subject_source, subject_label, outputs)?;
        }
        if let Some(first_batch_messages) = reports.reader_warnings.first() {
            outputs.diagnostics.write_all(first_batch_messages)?;
        }
        write_request_posts(args, &mut *outputs.diagnostics)
    })();
    match written {
        Ok(()) => crate::blastinput::app::engine_error(text),
        Err(error) => error,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_options_local_priv.hpp:1718-1725
/// The `EffectiveSearchSpace` option of the call's batch whose first query is `first` (BLASTP has
/// one context per query): its context's search space (`SplitQuery_SetEffectiveSearchSpace`), or
/// for a context that is not valid `-searchsp` (0 without it). The archive's
/// (`FirstContext::context0_space`) and the search strategy's (`CallRequests::set_searched`).
fn context0_search_space(
    query_valid: &[bool],
    searchsp: Option<i64>,
    exact_search_spaces: &[Option<i64>],
    search_spaces: &[SearchSpace],
    first: usize,
) -> i64 {
    if !query_valid[first] {
        return searchsp.unwrap_or(0);
    }
    exact_search_spaces[first].unwrap_or(search_spaces[first].effective_space.round() as i64)
}

/// The batch's context 0 as the `EffectiveSearchSpace` option of the archive and the search
/// strategy reads it (`archive::RequestLength::space`): its search space
/// (`context0_search_space`) and its query length, 0 for a context that is not valid.
fn context0_searched(
    contexts: &[QueryContext],
    query_valid: &[bool],
    searchsp: Option<i64>,
    exact_search_spaces: &[Option<i64>],
    search_spaces: &[SearchSpace],
    first: usize,
) -> crate::report::archive::SearchedContext0 {
    crate::report::archive::SearchedContext0 {
        space: context0_search_space(
            query_valid,
            searchsp,
            exact_search_spaces,
            search_spaces,
            first,
        ),
        query_length: if query_valid[first] {
            contexts[first].aa_len as i64
        } else {
            0
        },
    }
}

/// What the blocks and statistics of a batch's first query are made of (`FirstContexts::first_context`).
struct FirstContexts<'a> {
    contexts: &'a [QueryContext],
    query_valid: &'a [bool],
    searchsp: Option<i64>,
    exact_search_spaces: &'a [Option<i64>],
    search_spaces: &'a [SearchSpace],
    /// RPS-BLAST's `EffectiveSearchSpace` option of each batch, by its first query
    /// (`rps_request_spaces`); `None` otherwise.
    request_spaces: Option<&'a [Option<i64>]>,
    /// rpsblast's gapped K per query (`rps_search::engine::batch_gapped_k`); `None` otherwise.
    rps_gapped_k: Option<&'a [f64]>,
    gapped_params: crate::stats::KarlinParams,
    ungapped: bool,
}

impl FirstContexts<'_> {
    /// The blocks and statistics of the batch whose first query is `first` (BLASTP has one
    /// context per query): `-1` triples for a query without a valid context, and for the gapped
    /// block of an ungapped search.
    fn first_context(&self, first: usize) -> super::archive::FirstContext {
        let context0 = context0_searched(
            self.contexts,
            self.query_valid,
            self.searchsp,
            self.exact_search_spaces,
            self.search_spaces,
            first,
        );
        // rpsblast's option is the space its last search unit leaves in the options handle.
        let option_space = self
            .request_spaces
            .and_then(|spaces| spaces[first])
            .unwrap_or(context0.space);
        if !self.query_valid[first] {
            return super::archive::FirstContext {
                karlin: None,
                search_space: 0,
                length_adjustment: 0,
                context0_space: option_space,
                context0_length: context0.query_length,
            };
        }
        // rpsblast's first query of a batch has its profile's K (`RPS_K_MULT` times the `.aux`
        // value, `rps_search::engine::batch_gapped_k`), the others the table's.
        let mut gapped = self.gapped_params;
        if let Some(k) = self.rps_gapped_k {
            gapped.k = k[first];
        }
        super::archive::FirstContext {
            karlin: Some((
                self.contexts[first].karlin_params,
                (!self.ungapped).then_some(gapped),
            )),
            search_space: context0.space,
            length_adjustment: self.search_spaces[first].length_adjustment,
            context0_space: option_space,
            context0_length: context0.query_length,
        }
    }
}

/// rpsblast's `EffectiveSearchSpace` option of each batch (by its first query `first` of
/// `firsts`), the space its last search unit leaves in the options handle
/// (`rps_hooks::request_space`); `None` for another search.
fn rps_request_spaces(
    firsts: impl Iterator<Item = usize>,
    rps_units: Option<&crate::algorithm::rps_search::volumes::RpsUnits>,
    args: &ResolvedBlastpArgs,
    contexts: &[QueryContext],
    gapped_params: &KarlinParams,
) -> Result<Option<Vec<Option<i64>>>> {
    let Some(units) = rps_units else {
        return Ok(None);
    };
    let mut spaces = vec![None; contexts.len()];
    for first in firsts {
        spaces[first] = Some(rps_hooks::request_space(
            units,
            args,
            &contexts[first],
            gapped_params,
        )?);
    }
    Ok(Some(spaces))
}

/// The searched space of each batch of `call` (`CallRequests::set_searched`): its context 0
/// (`context0_searched`) where NCBI searches the batch (not `query_batch_skipped`), with
/// rpsblast's the space its last search unit leaves in the options handle
/// (`rps_hooks::request_space`).
fn set_call_searched(
    call: &CallRequests<'_>,
    args: &ResolvedBlastpArgs,
    rps_units: Option<&crate::algorithm::rps_search::volumes::RpsUnits>,
    contexts: &[QueryContext],
    gapped_params: &KarlinParams,
    (query_valid, query_batch_skipped): (&[bool], &[bool]),
    (exact_search_spaces, search_spaces): (&[Option<i64>], &[SearchSpace]),
) -> Result<()> {
    let rps_request_spaces = rps_request_spaces(
        call.batch_firsts(),
        rps_units,
        args,
        contexts,
        gapped_params,
    )?;
    call.set_searched(|first| {
        (!query_batch_skipped[first]).then(|| {
            let mut context0 = context0_searched(
                contexts,
                query_valid,
                args.searchsp,
                exact_search_spaces,
                search_spaces,
                first,
            );
            if let Some(space) = rps_request_spaces.as_ref().and_then(|spaces| spaces[first]) {
                context0.space = space;
            }
            context0
        })
    });
    Ok(())
}

/// The batch whose search strategy request the run writes right after a call's batches are
/// reported (`BlastpReports::strategy_next`): the batch that stops the run (every query without
/// letters) after the searched batches, or the batch that fails with a header file NCBI cannot read
/// after the batches before it, which are searched again.
#[derive(Clone, Copy)]
struct NextRequest<'a> {
    /// The run's index of the batch's first query.
    first_query: usize,
    records: &'a [FastaRecord],
    /// Where the first query's letters lie in its record (`-query_loc`).
    placement: Option<seq_range::RecordPlacement>,
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-262
/// The search strategy's requests (`-export_search_strategy`) of a call's batches. NCBI writes a
/// batch's request before its search, and BLASTP searches a call's batches together, so the call
/// writes the requests where it ends (`end_at`: every batch before the formats, or the batch an
/// NCBI error of the call belongs to), through the run's cursor (`StrategyExport::write_with`: a
/// request written already, by a search of the batches again, is not written twice), each
/// followed by its batch's searched space, and marks that place; `run_resolved_in_pool` releases
/// the call's requests in one place (`release`): on success every batch's (those that no mark
/// wrote) and the request of the batch that follows the call (`NextRequest`). A call that ends
/// with an NCBI error without a mark is LOSAT's internal error, so that no error path drops
/// requests silently. NCBI's posts of the bytes replaced in a request wait in the run's export
/// until the run writes the start of the request's batch (`StrategyExport::take_posts`): the
/// formats place those of the call's batches before their first queries' reports (`take_posts`),
/// and an error of a batch writes those of its batch after its reader messages
/// (`write_request_posts`), wherever the request was written (a search of the batches again
/// places those of an earlier attempt).
struct CallRequests<'a> {
    strategy: crate::report::archive::strategy::StrategyExport,
    requests: crate::report::archive::requests::RunRequests<'a>,
    /// The call's queries and their placements, from the run's query `base`.
    records: &'a [FastaRecord],
    placements: &'a seq_range::Placements,
    base: usize,
    /// The call's batches (call-local, not empty), and each one's searched space once the call's
    /// set-up has computed it (`set_searched`).
    spans: Vec<std::ops::Range<usize>>,
    searched: std::sync::Mutex<Vec<Option<crate::report::archive::SearchedContext0>>>,
    /// The call's query where the requests NCBI wrote end (`end_at`).
    through: std::sync::Mutex<Option<usize>>,
    next: Option<NextRequest<'a>>,
}

impl<'a> CallRequests<'a> {
    /// The requests of the call of `query_records` (`ranges`: their batches and the run's index of
    /// the first), with the batch that follows it (`next`). `None` when the batch size cannot be
    /// read: the call stops with that error before any batch.
    fn of_call(
        strategy: crate::report::archive::strategy::StrategyExport,
        args: &ResolvedBlastpArgs,
        query_records: &'a [FastaRecord],
        subject_source: BlastpSubjectSource<'a>,
        ranges: &'a BlastpRanges<'a>,
        next: Option<NextRequest<'a>>,
    ) -> Option<Self> {
        let batch_size = args.query_batch_size().ok()?;
        let spans: Vec<std::ops::Range<usize>> = ranges
            .query_input
            .batches(batch_size)
            .into_iter()
            .map(|batch| batch.searched)
            .filter(|searched| !searched.is_empty())
            .collect();
        Some(Self {
            strategy,
            requests: super::archive::requests(
                subject_source
                    .database()
                    .is_none()
                    .then(|| subject_source.fasta_records()),
                args.lcase_masking,
                args.psi_search.as_ref().and_then(|psi| psi.request_pssm()),
            ),
            records: query_records,
            placements: &ranges.query_input.placements,
            base: ranges.searched_base,
            searched: std::sync::Mutex::new(vec![None; spans.len()]),
            spans,
            through: std::sync::Mutex::new(None),
            next,
        })
    }

    /// The call-local index of each batch's first query, in order.
    fn batch_firsts(&self) -> impl Iterator<Item = usize> + '_ {
        self.spans.iter().map(|span| span.start)
    }

    /// The searched space of each batch (`searched(first query)`, `None` for a batch NCBI does not
    /// search), once the call's set-up has computed them.
    fn set_searched(
        &self,
        searched: impl Fn(usize) -> Option<crate::report::archive::SearchedContext0>,
    ) {
        let mut spaces = self
            .searched
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        for (space, span) in spaces.iter_mut().zip(&self.spans) {
            *space = searched(span.start);
        }
    }

    /// The call ends at its query `through`: NCBI has written the requests of the batches before
    /// it (every batch before the formats; the batch an NCBI error belongs to), so they are
    /// written now (`write_through`), before the output of the call's end, and the place is
    /// marked. NCBI's posts of the bytes that ASN.1 text replaced in them wait for the starts of
    /// their batches (`take_posts`, `write_request_posts`). A request whose stream throws (an
    /// element stated as `NULL` in a PSSM query, `encoder::NullElement`) ends the call there too:
    /// NCBI's error follows it, so the place is marked before the error returns (`release`).
    fn end_at(&self, through: usize) -> Result<()> {
        let written = self.write_through(through);
        if matches!(&written, Err(error) if !error.is::<crate::cli::NativeError>()) {
            return written;
        }
        *self
            .through
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner) = Some(through);
        written
    }

    /// Takes NCBI's posts of the bytes replaced in the requests written of the call's batches,
    /// each with the call's index of its batch's first query (`StrategyExport::take_posts`).
    fn take_posts(&self) -> Vec<(usize, crate::report::diag_log::Posts)> {
        self.strategy
            .take_posts(self.base + self.records.len())
            .into_iter()
            .map(|(first, posts)| (first.saturating_sub(self.base), posts))
            .collect()
    }

    /// The end of the call's first batch, where an NCBI error of that batch's reports ends the
    /// call (`end_at`).
    fn first_batch_end(&self) -> usize {
        self.spans.first().map_or(0, |span| span.end)
    }

    /// Writes the requests of the call's batches that start before its query `through` (those
    /// written already are not made again), each followed by its batch's searched space.
    fn write_through(&self, through: usize) -> Result<()> {
        let searched = self
            .searched
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .clone();
        for (span, &space) in self.spans.iter().zip(&searched) {
            if span.start >= through {
                break;
            }
            let queries = self.base + span.start..self.base + span.end;
            self.requests.write(
                &self.strategy,
                queries.start,
                &self.records[span.clone()],
                self.placements.get(span.start),
            )?;
            self.strategy.note_searched(queries, space)?;
        }
        Ok(())
    }

    /// Writes the requests of the call that ended with `result` (`end_at` wrote those up to its
    /// mark): on success every batch's and the next one's, whose posts the caller writes after
    /// that batch's reader messages (`write_request_posts`); after an error those up to the
    /// call's mark. A call that stops with a header file NCBI cannot read writes none (`search`
    /// searches the batches before it again), and so does a refusal of LOSAT's.
    fn release(&self, result: &Result<()>) -> Result<()> {
        let mark = *self
            .through
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let through = match (result, mark) {
            (Ok(()), _) => self.records.len(),
            (Err(_), Some(through)) => through,
            (Err(error), None) => {
                anyhow::ensure!(
                    !error.is::<crate::cli::NativeError>(),
                    "LOSAT internal error: a blastp call ended with NCBI's error without the batch where its search strategy requests end ({error})"
                );
                return Ok(());
            }
        };
        self.write_through(through)?;
        if let (Ok(()), Some(next)) = (result, self.next) {
            self.requests.write(
                &self.strategy,
                next.first_query,
                next.records,
                next.placement,
            )?;
        }
        Ok(())
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:259-276
/// Writes NCBI's posts of the bytes replaced in the search strategy requests written whose
/// batches' starts are not written (`StrategyExport::write_posts`): those of the batch whose
/// error follows, after the batch's reader messages and before the error.
fn write_request_posts(args: &ResolvedBlastpArgs, diagnostics: &mut dyn Write) -> Result<()> {
    if let Some(strategy) = args.strategy.as_ref() {
        strategy.write_posts(diagnostics, args.archive.as_ref())?;
    }
    Ok(())
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:177-180
/// The error `text` of a query batch's search (`BLAST engine error: <text>`, exit status 3) when
/// the call's first batch holds invalid queries only: NCBI does not search such a batch and
/// reports its queries, so the reports of the batches before the first searched one come first
/// (searched again here, alone and without the epilog: none of their queries is searched), then
/// that batch's reader messages and the error (NCBI 2.17.0 blastp `BATCH_SIZE=1`, an all-X query
/// before tk, with `-subject` of `MK`, ``, `` and `-searchsp 1000`, or `-db` of an empty
/// database: outfmt 0 and 7, exit 3; line PS DECISIONS 28 and 30). LOSAT refuses where this is
/// not ported: no searched batch, another program, or a format other than 0 to 4, 6, 7 and 10.
#[cold]
#[inline(never)]
#[allow(clippy::too_many_arguments)]
fn unsearched_batches_then_engine_error(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    pool: &crate::utils::threading::SearchPool<'_>,
    ranges: &BlastpRanges<'_>,
    reports: &BlastpReports<'_>,
    query_batches: &[seq_range::RangedBatch],
    query_batch_skipped: &[bool],
    database_name: &str,
    text: &str,
) -> anyhow::Error {
    let first_searched = query_batches
        .iter()
        .find(|batch| !batch.searched.is_empty() && !query_batch_skipped[batch.searched.start]);
    let ported_formats = outputs.formats.iter().all(|format| {
        crate::blastinput::app::parse_formatting_string(format.outfmt).is_ok_and(|choice| {
            matches!(choice.number, 6 | 7 | 10)
                || blastp_pairwise_view_of(args.app, choice.number).is_some()
        })
    });
    let Some(batch) = first_searched
        .filter(|_| ported_formats && args.app == BlastpApp::Blastp && args.archive.is_none())
    else {
        let searched = if database_name.is_empty() {
            "the subjects".to_string()
        } else {
            format!("the database {database_name}")
        };
        return anyhow::anyhow!(
            "a search of {searched} that fails with NCBI BLAST+'s error \"{text}\" after a first query batch that NCBI BLAST+ does not search is not supported by LOSAT's {}",
            args.app.name()
        );
    };
    let (searched_end, input_end) = (batch.searched.start, batch.input.start);
    let input = ranges.query_input;
    let prefix_input = seq_range::QueryInput {
        placements: seq_range::Placements(
            input
                .placements
                .0
                .iter()
                .take(searched_end)
                .copied()
                .collect(),
        ),
        ordinals: input.ordinals[..searched_end].to_vec(),
        input_lengths: input.input_lengths[..input_end].to_vec(),
        skipped: input.skipped[..input_end].to_vec(),
        range_length: input.range_length,
    };
    let prefix_ranges = BlastpRanges {
        input_records: &ranges.input_records[..input_end],
        query_input: &prefix_input,
        subject_placements: ranges.subject_placements,
        x_ordinal_base: ranges.x_ordinal_base,
        searched_base: ranges.searched_base,
        split_search_spaces: None,
    };
    // The first searched batch's search strategy request follows those of the batches before
    // (`CallRequests`), written before their reports; the reports of those batches place NCBI's
    // posts of the bytes replaced in their requests, and the failed batch's follow its reader
    // messages.
    if let Some(call) = reports.call_requests {
        if let Err(error) = call.end_at(batch.searched.end) {
            return error;
        }
    }
    let prefix_reports = BlastpReports {
        reader_warnings: &reports.reader_warnings
            [..searched_end.min(reports.reader_warnings.len())],
        identity_warned: reports.identity_warned,
        epilog: false,
        x_prolog: reports.x_prolog,
        x_queries_before: reports.x_queries_before,
        x_xml2_carry: reports.x_xml2_carry,
        // The engine error follows this report and ends the run: it is the run's last report
        // call (the XML formats, the only readers of `last_call`, keep the refusal before here).
        last_call: true,
        psi_round: None,
        // BLASTP only (`args.app` above): no psiblast PSSM stops here.
        stopped_pssms: None,
        query_scopes: reports.query_scopes,
        prefetch_scope: None,
        strategy_next: None,
        call_requests: None,
        x_plan_blocks: false,
    };
    if let Err(error) = run_resolved_in_pool(
        args.clone(),
        &query_records[..searched_end],
        subject_source,
        "",
        subject_label,
        outputs,
        pool,
        None,
        &prefix_ranges,
        &prefix_reports,
        x_lookup_ahead::XAheadBlock::none(),
    ) {
        return error;
    }
    if let Some(messages) = reports.reader_warnings.get(searched_end) {
        if let Err(error) = outputs.diagnostics.write_all(messages) {
            return error.into();
        }
    }
    if let Err(error) = write_request_posts(args, &mut *outputs.diagnostics) {
        return error;
    }
    crate::blastinput::app::engine_error(text)
}

/// What the reports of the searched batches need besides the search: the reader's
/// messages before each query's report (`read_query_batches`) and whether the reports end
/// with the epilog.
#[derive(Clone, Copy)]
struct BlastpReports<'a> {
    reader_warnings: &'a [Vec<u8>],
    /// The searched queries of the run before this one have IDENTITY's warning
    /// (`report::query_warnings::IDENTITY_CBS_RESET_MESSAGE`): those of the first query batch
    /// when IDENTITY reset the composition mode, none otherwise.
    identity_warned: usize,
    epilog: bool,
    /// EXPERIMENT (LOSAT_X_BPQBLOCK): whether the outfmt 0 report starts with the prolog
    /// (false for the blocks after the first of the block search, `x_query_blocks.rs`).
    x_prolog: bool,
    /// EXPERIMENT (LOSAT_X_BPQBLOCK): the queries of the blocks before this one, which outfmt
    /// 7's last line counts with this call's (0 outside the block search).
    x_queries_before: usize,
    /// EXPERIMENT (LOSAT_X_BPQBLOCK): the framing state of an outfmt 13 to 16 run that the
    /// blocks hand on (the count of reports, the start of the document, a `-subject` run's held
    /// report; `None` outside the block search).
    x_xml2_carry: Option<&'a std::sync::Mutex<crate::report::blast_xml2::Xml2Carry>>,
    /// Whether this call is the run's last report call, so that without the epilog the run stops
    /// after it (`blast_xml::BlastpXml::stops`). Not a block of the block search before the last
    /// (LOSAT_X_BPQBLOCK), a chunk's search (which reports nothing) or a psiblast round (whose
    /// XML document the rounds end, `psi_round::write_round_epilogs`).
    last_call: bool,
    /// psiblast's round of a FASTA query (`psi_round`): the round's number and previous
    /// inclusion set for its reports, and where its outcome goes; `None` outside the rounds.
    psi_round: Option<&'a crate::algorithm::psiblast::rounds::PsiRoundIo>,
    /// The PSSM files and read messages of the queries up to a psiblast PSSM's stop
    /// (`psiblast::pssm::PssmStop`), for the search again of those queries (`search`): the call
    /// takes them for its own (`psiblast_pssms` makes none), so the PSSMs' reads look their taxids
    /// up after each query's report as NCBI's did (`ListedTitles::pssm_load_texts`) and their
    /// messages come where a run that does not stop writes them. `None` (or empty) elsewhere.
    stopped_pssms:
        Option<&'a std::sync::Mutex<Option<crate::algorithm::psiblast::pssm::PssmFiles>>>,
    /// Under `-parse_deflines` against a database, the scopes of the call's queries
    /// (`query_batch::QueryScopes`, numbered from the call's first), whose ids the reported
    /// subjects must not have (`ListedTitles::refuse_ids_of_queries`); `None` elsewhere.
    query_scopes: Option<crate::blastinput::query_batch::QueryScopesFrom<'a>>,
    /// psiblast's round of a query (`psi_round`), a call of its own one record: the query records
    /// that NCBI's scope holds at the round's prefetch (the query's batch, the records that
    /// `-query_loc` skips before it included, or the PSSM's query of `-in_pssm`), made only where
    /// the prefetch stops (`common::database::PrefetchQueries`); `None` elsewhere: the call's
    /// batches are its input records (`report_prep::prefetch_reported_subjects`).
    prefetch_scope: Option<&'a (dyn Fn() -> crate::blastinput::query_batch::QueryScope + Sync)>,
    /// `-export_search_strategy`: the batch whose request the run writes right after this call's
    /// batches are reported (`NextRequest`), set by the caller.
    strategy_next: Option<NextRequest<'a>>,
    /// `-export_search_strategy`: the call's requests (`CallRequests`), set by
    /// `run_resolved_in_pool` for its own call.
    call_requests: Option<&'a CallRequests<'a>>,
    /// EXPERIMENT (LOSAT_X_BPQBLOCK): whether the call may search its preliminary stage in blocks
    /// of whole batches (`x_prelim_block_plan`): the search of all queries of `search` only, not a
    /// block of the block search (whose plan was made for all queries), a query chunk, a block of
    /// the preliminary stage, a psiblast round or the batches before a failing one.
    x_plan_blocks: bool,
}

impl<'a> BlastpReports<'a> {
    /// The reports of a call of the preliminary stage only (`PrelimOnly`): none (the call of the
    /// whole query or of all queries writes them, and the archive). `x_prolog`: an error of the
    /// call's search is written after the outfmt 0 prolog (`preliminary_search_error`).
    fn prelim_only(x_prolog: bool) -> Self {
        Self {
            reader_warnings: &[],
            identity_warned: 0,
            epilog: false,
            x_prolog,
            x_queries_before: 0,
            x_xml2_carry: None,
            last_call: false,
            psi_round: None,
            stopped_pssms: None,
            query_scopes: None,
            prefetch_scope: None,
            strategy_next: None,
            call_requests: None,
            x_plan_blocks: false,
        }
    }
}

/// The environment that NCBI's blastp reads around its search that LOSAT rejects:
/// `BL2SEQ_LEGACY` in a search of subject records (`fasta_subjects`). The values that NCBI
/// converts in each query batch stop its first batch: `CHUNK_SIZE` and `OVERLAP_CHUNK_SIZE` in its
/// query splitter (`read_query_batches`; not `OVERLAP_CHUNK_SIZE` in RPSBLAST, whose queries it
/// never splits) and `PRE_FETCH_SEQS_LIMIT` in its prefetch (`run_call`). Other batch
/// sizes change a blastp report only through the query batches (`read_query_batches`; the redo of
/// a query reads the heaps of its batch, `kappa_query_batches`) and the query split
/// (`search_split_query_batches`).
fn check_blastp_environment(app: BlastpApp, fasta_subjects: bool) -> Result<()> {
    crate::blastinput::app::check_search_environment(app.name(), fasta_subjects)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:145-188
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-283
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:215-219
// NCBI searches every query batch on its own (RP-4). A blastp query's preliminary hit list
// does not depend on the other queries of its batch, so the batches that NCBI does not
// split are the search of all queries above; the queries of a batch that NCBI splits
// (`blastp/protein_query_split.rs`) take the merged hit lists of the batch's query chunks,
// each searched as a query set of its parts. NCBI gives a chunk the options' search spaces of
// its batch's contexts (`SplitQuery_SetEffectiveSearchSpace`: one context per query, 0 for a
// query whose Karlin-Altschul block failed on the whole query; `option_space_of_query`,
// `BlastpRanges::split_search_spaces`), which the chunk's setup reads by the context's index in
// the chunk (a part of another query than the batch's query of that index, in general), and
// the chunk's own length adjustments; a context given 0 gets the space of its own part. With
// Spouge's FSC the preliminary cutoffs and E-values do not read them; psiblast's `-sum_stats
// true` links with them, and without the Gumbel block (`OLD_FSC`) the cutoffs read them. A
// chunk whose contexts all fail the Karlin-Altschul setup finds nothing, as NCBI ignores its
// error.
//
// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:219-221
// NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:173-181
#[allow(clippy::too_many_arguments)]
fn search_split_query_batches(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    input_records: &[FastaRecord],
    subject_source: BlastpSubjectSource<'_>,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
    query_lengths: &[usize],
    option_space_of_query: &dyn Fn(usize) -> Result<i64>,
    query_input: &seq_range::QueryInput,
    preliminary_hit_lists: &mut [Option<BlastpHitList>],
    x_compact_hit_lists: &mut [Option<BlastpHitList<XPrelimHsp>>],
    x_prolog: bool,
    call: Option<&CallRequests<'_>>,
) -> Result<()> {
    use crate::algorithm::blastp::protein_query_split::{
        protein_chunk_part_in_record, protein_chunk_would_be_split, split_protein_batch,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:86-93
    // An RPS-BLAST query batch is never split.
    if args.rps_search.is_some() {
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:61
    // An ungapped search splits no query batch (and so builds one compressed lookup table for
    // each batch, `ncbi_compressed_lookup_tables`).
    if args.ungapped {
        return Ok(());
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:122-129 (cited at
    // `psi_round::run_queries`: a query factory that is not `CObjMgr_QueryFactory` makes one
    // chunk)
    // A PSSM search splits only a later round of a FASTA query (`PsiSearch::split_setup`); a PSSM
    // query (`-in_pssm`, `CObjMgrFree_QueryFactory`) is never split.
    if args
        .psi_search
        .as_deref()
        .is_some_and(|psi| psi.split_setup.is_none())
    {
        return Ok(());
    }
    let sizes = protein_query_split_sizes(10000)?;
    let batch_size = args.query_batch_size()?;
    // The batches count whole records (`seq_range::QueryInput::batches`); a batch is
    // split by the lengths of its searched queries (their intervals with `-query_loc`).
    for range in query_input
        .batches(batch_size)
        .into_iter()
        .map(|batch| batch.searched)
        .filter(|searched| !searched.is_empty())
    {
        let Some(chunks) = split_protein_batch(&query_lengths[range.clone()], sizes) else {
            // A later round that the round driver split (`psi_round::run_queries`, the same
            // lengths and sizes) has no PSSM words to search: its chunks must run here.
            anyhow::ensure!(
                args.psi_search
                    .as_deref()
                    .is_none_or(|psi| psi.split_setup.is_none()),
                "internal error of LOSAT's {}: a later round set up as split has no query chunk",
                args.app.name()
            );
            continue;
        };
        // NCBI stops in the search of the batch (its exception) after it has printed the reports of
        // the batches before (`blastp_app.cpp:257-276`, cited in `search`): so do the refusals
        // below, as the batch's failure (`database::HeaderReadError`, the batch's first query;
        // TBLASTN's and TBLASTX's are `blastinput::app::ErrorAfterBatchReports`). psiblast searches
        // one query at a time and NCBI's stops in the query's round 1 (it splits that round too),
        // after the rounds of the queries before: every round of a query makes the same chunks,
        // so a later round never reaches these refusals, and the refusal of the query's search is
        // the plain error.
        let batch_refusal = |refusal: anyhow::Error| -> anyhow::Error {
            if args.psi_search.is_some() {
                refusal
            } else {
                // Before the batch's search: nothing of it read a header.
                database::HeaderReadError {
                    query: range.start,
                    error: refusal,
                    utf8: None,
                }
                .into()
            }
        };
        // NCBI's `CQuerySplitter::Split` makes the query factory of every chunk before it sets any
        // chunk up, and a chunk without a query has no query vector to make it from
        // (`CCoreException::eNullPtr`), so NCBI stops there before a chunk is split again (as with
        // a positive pair whose evened chunk size falls below the overlap, or a negative pair whose
        // wrapped size passes the batch's length: chunk 0 covers the batch); `split_protein_batch`
        // makes no part for such a chunk.
        //
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:234,876-879
        if chunks.iter().any(|chunk| chunk.queries.is_empty()) {
            return Err(batch_refusal(anyhow::anyhow!(
                "a query chunk without a query is not supported by LOSAT's {}",
                args.app.name()
            )));
        }
        // Every chunk holds a query: NCBI sets the chunks up one by one after the split
        // (prelim_stage.cpp:235-244) and stops at the missing lookup table of a chunk that it would
        // split again (prelim_search_runner.hpp:77).
        if protein_chunk_would_be_split(&chunks, sizes) {
            return Err(batch_refusal(anyhow::anyhow!(
                "an OVERLAP_CHUNK_SIZE so close to CHUNK_SIZE that NCBI BLAST+ would split a query chunk again (it then stops with a null reference) is not supported by LOSAT's {}",
                args.app.name()
            )));
        }
        // A negative CHUNK_SIZE splits a batch only with a negative OVERLAP_CHUNK_SIZE below it,
        // and NCBI's `size_t` chunk ranges then wrap (split_query_cxx.cpp:145-171, cited at
        // `split_protein_batch`): where the first range covers the batch, the chunks after it
        // hold no query and NCBI stops above. Otherwise the ranges leave gaps between chunks that
        // each hold a query and would not be split again, which NCBI searches on with the hits in
        // the gaps lost: the maintainer kept that setting refused (PD-LOSAT-NCBI-DEFECTS version
        // 1.1), and the refusal comes where the batch's search would start, after the reports of
        // the batches before (the form of 1.40 for a refusal in the middle of a run).
        if sizes.negative_chunk_size() {
            return Err(batch_refusal(anyhow::anyhow!(
                "a negative CHUNK_SIZE above a negative OVERLAP_CHUNK_SIZE splits this query batch with NCBI BLAST+'s size_t chunk ranges wrapped, which is not supported by LOSAT's {}",
                args.app.name()
            )));
        }
        let mut merged: Vec<Option<BlastpHitList>> = (0..range.len()).map(|_| None).collect();
        // The options' search spaces of the batch's contexts (one per query), which every chunk
        // reads by its own context index (`SplitQuery_SetEffectiveSearchSpace`, cited at
        // `stats::setup::preliminary_option_space`): computed once for the batch, 0 for a query
        // whose Karlin-Altschul block failed on the whole query (unless `-searchsp` gives it).
        let batch_option_spaces: Vec<i64> = range
            .clone()
            .map(option_space_of_query)
            .collect::<Result<_>>()?;
        for chunk in &chunks {
            let parts: Vec<FastaRecord> = chunk
                .queries
                .iter()
                .map(|part| {
                    let query = range.start + part.query;
                    let mut cut = query_records[query].cut(part.from, part.to);
                    // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:278-281
                    // `-lcase_masking`: the chunk's lower case is the runs of the query's whole
                    // record restricted to the part as NCBI restricts them (each run one letter
                    // longer, and with `-query_loc` in the part's coordinates without the
                    // interval's start).
                    if args.lcase_masking {
                        cut.sequence = protein_chunk_part_in_record(
                            input_records[query_input.ordinal(query)].seq(),
                            query_input.placements.offset(query),
                            part,
                        );
                    }
                    cut
                })
                .collect();
            let mut chunk_lists = Vec::new();
            let part_input = seq_range::QueryInput::whole(&parts);
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_aux_priv.cpp:191-196
            // A chunk is set up without the PSSM: a split later round's chunks are searched as
            // round 1's, with the round's options (its composition-based statistics mode and
            // search spaces).
            let mut chunk_args = args.clone();
            chunk_args.psi_search = None;
            let searched = run_resolved_in_pool(
                chunk_args,
                &parts,
                subject_source,
                "",
                subject_label,
                outputs,
                blastp_parallel_pool,
                Some(PrelimOnly::Chunk(&mut chunk_lists)),
                &BlastpRanges {
                    input_records: &parts,
                    query_input: &part_input,
                    subject_placements: &seq_range::Placements::default(),
                    x_ordinal_base: 0,
                    searched_base: 0,
                    split_search_spaces: Some(&batch_option_spaces),
                },
                // A chunk's search ends after its preliminary stage, before any report.
                // EXPERIMENT (LOSAT_X_BPQBLOCK): a chunk of a block after the first follows the
                // first block's outfmt 0 prolog (`preliminary_search_error`).
                &BlastpReports::prelim_only(x_prolog),
                x_lookup_ahead::XAheadBlock::none(),
            );
            // A chunk's search that fails is its batch's (`preliminary_search_error`): the
            // requests of the call's batches up to it come before the error (NCBI searched the
            // batches before it, one at a time, and wrote their requests). No NCBI error reaches
            // here today: a chunk reads the database only for a valid context of its parts, whose
            // query is valid in the call, so the call's own search above has read the same
            // database and made the same checks first (`first_oid_read_fails`, the reads,
            // `OLD_FSC`'s average length; RG K4b review M1b, probe K/k4b/agent6/m1b: an all-X query
            // longer than a chunk, a query whose second chunk is all X, an empty database and one
            // without its sequence file, 28 cases equal to NCBI). Such an error would write no
            // reader messages of its batch, and the reports of the batches before it would be
            // missing.
            if let Err(error) = searched {
                if let Some(call) = call {
                    call.end_at(range.end)?;
                    write_request_posts(args, &mut *outputs.diagnostics)?;
                }
                return Err(error);
            }
            let chunk_parts: Vec<super::query_split::BlastpChunkPart> = chunk
                .queries
                .iter()
                .zip(&chunk.context_offsets)
                .map(|(part, &offset)| super::query_split::BlastpChunkPart {
                    batch_query: part.query,
                    q_idx: (range.start + part.query) as u32,
                    offset,
                    query_length: query_lengths[range.start + part.query],
                })
                .collect();
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/split_query_cxx.cpp:174-178
            super::query_split::merge_query_chunk(
                &mut merged,
                chunk_lists,
                &chunk_parts,
                sizes.overlap as u32 as i32,
            );
        }
        for (local, list) in merged.into_iter().enumerate() {
            // The chunks' list replaces the compact one too.
            let q_idx = range.start + local;
            x_compact_hit_lists[q_idx] = None;
            preliminary_hit_lists[q_idx] = list;
        }
    }
    Ok(())
}

// EXPERIMENT (LOSAT_X_BPLKPAHEAD): the next query block's query frames, lookup table and contexts
// built during the block before it (`x_query_blocks.rs` block search).
pub(super) mod x_lookup_ahead;

// `-lcase_masking`: the word scan of a `-subject` with lower-case letters over its scan ranges.
mod x_soft_subject;
// BLASTP's sequences for the query-anchored reports (outfmt 1 to 4).
mod anchored;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_filter.c:1381-1382
/// The query frames of a search: each record encoded with the search's SEG parameters
/// (`encode_protein_query_frame_with_seg`; `encode_protein_query_frame_soft_masked` with
/// `-soft_masking true`; `encode_protein_query_frame_lcase` with `-lcase_masking`).
/// `run_call` encodes its queries with it, and so
/// does the build of the next block's lookup ahead of that block (LOSAT_X_BPLKPAHEAD,
/// `x_lookup_ahead.rs`).
fn blastp_query_frames(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
) -> Vec<Vec<QueryFrame>> {
    blastp_query_frames_marked(args, query_records, false)
}

/// `blastp_query_frames`; with `marked` and LOSAT_TIMING the encoding and the SEG go to the
/// common line (`encode_protein_query_frame_with_seg_timed`). Only the search's own call marks.
fn blastp_query_frames_marked(
    args: &ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    marked: bool,
) -> Vec<Vec<QueryFrame>> {
    let seg_params = args.seg.params();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:621-624
    query_records
        .iter()
        .map(|record| {
            vec![if args.lcase_masking {
                encode_protein_query_frame_lcase(
                    record.seq(),
                    seg_params.as_ref(),
                    args.soft_masking,
                    marked,
                )
            } else if args.soft_masking {
                encode_protein_query_frame_soft_masked(record.seq(), seg_params.as_ref(), marked)
            } else if marked {
                encode_protein_query_frame_with_seg_timed(record.seq(), seg_params.as_ref())
            } else {
                encode_protein_query_frame_with_seg(record.seq(), seg_params.as_ref())
            }]
        })
        .collect()
}

// No NCBI counterpart: reads the LOSAT_X_BPLK16 switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPLK16): strict set (on unless `0`). Exact blastp's block tables hold 16-bit hits in 8-byte cells where the lookup query is at most 65,535
/// letters (`BlastAaLookupTable<BackboneCell16>`, lookup/backbone.rs), halving the overflow array
/// and the cells; the hits are the same in the same order.
fn x_bplk16() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPLK16"))
}

// No NCBI counterpart: chooses LOSAT's table (LOSAT_X_BPLK16); it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPLK16): whether a table of `args`' search may have 16-bit hits: a search
/// with a window (the fused two-hit loop scans it; the offset-pair scan of a search without one
/// reads the 32-bit table) that is not a query chunk of a split batch (`x_chunk`).
fn x_bplk16_for(args: &ResolvedBlastpArgs, x_chunk: bool) -> bool {
    args.window_size != 0 && !x_chunk
}

// No NCBI counterpart: LOSAT_X_BPLK16's table in the runtime lookup; it does not change any value
// NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLK16): the runtime lookup of a table built with the switch on, counted for
/// `LOSAT_X_REACH` (taken: the table has 16-bit hits).
fn x_bp_runtime_lookup(built: XBpLookup) -> BlastpRuntimeLookup {
    crate::utils::xreach::x_reach!("LOSAT_X_BPLK16", matches!(built, XBpLookup::Narrow(_)));
    match built {
        XBpLookup::Wide(table) => BlastpRuntimeLookup::Aa(table),
        XBpLookup::Narrow(table) => BlastpRuntimeLookup::Aa16(table),
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:237,245
/// The lookup table of exact blastp (`AaLookupTable`) and the query contexts of `query_frames`, as
/// `call_lookup` (`blast_engine/call_setup.rs`) builds them, of the search's word size (3, 2:
/// blastp-short and `-word_size 2`, item SC3, or 4, item SC8b). The threshold takes the x86-64
/// conversion: `INT_MIN` for a threshold beyond `Int4` (every word). The table is filled in place,
/// without the reference builder's chains (lookup/x_lut_direct.rs); a table of 4-letter words by
/// the joint builder on `pool`'s threads (`x_bp_build_wide_lookup`). The build of the next block's table ahead of that block
/// (LOSAT_X_BPLKPAHEAD, `x_lookup_ahead.rs`, `blastp_aa_lookup_ahead`) takes the same threshold and
/// query frames and gives the same table; the Karlin parameters it is passed are not used.
fn blastp_aa_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    x_chunk: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    // `-word_size 4` (item SC8b): the joint builder of 32-bit indices on the pool's threads.
    if let Some((lookup, contexts)) = crate::algorithm::tblastx::lookup::x_bp_build_wide_lookup(
        query_frames,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        args.scoring.matrix,
        args.word_size,
        pool,
    ) {
        return (BlastpRuntimeLookup::Aa(lookup), contexts);
    }
    // EXPERIMENT (LOSAT_X_BPLK16): the table of 16-bit hits where the switch is on.
    if x_bplk16() {
        let (built, contexts) = x_bp_build_ncbi_lookup_as(
            query_frames,
            crate::core::blast_util::ncbi_int4_from_double(args.threshold),
            false,
            x_bplk16_for(args, x_chunk),
            args.scoring.matrix,
            args.word_size,
        );
        return (x_bp_runtime_lookup(built), contexts);
    }
    // Another matrix than BLOSUM62, or `-word_size 2` (item SC3): the same table of the matrix's
    // words of that size (`x_bp_build_ncbi_lookup` builds BLOSUM62's words of 3).
    if args.scoring.matrix != ScoringMatrix::Blosum62
        || args.word_size != crate::algorithm::tblastx::lookup::LOOKUP_WORD_LENGTH
    {
        return match x_bp_build_ncbi_lookup_as(
            query_frames,
            crate::core::blast_util::ncbi_int4_from_double(args.threshold),
            false,
            false,
            args.scoring.matrix,
            args.word_size,
        ) {
            (XBpLookup::Wide(lookup), contexts) => (BlastpRuntimeLookup::Aa(lookup), contexts),
            (XBpLookup::Narrow(_), _) => unreachable!("a build without `narrow` has Int4 hits"),
        };
    }
    let (lookup, contexts) = x_bp_build_ncbi_lookup(
        query_frames,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        false,
    );
    (BlastpRuntimeLookup::Aa(lookup), contexts)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:245
/// EXPERIMENT (LOSAT_X_BPLKPAHEAD): `blastp_aa_lookup`'s table and contexts (the same arguments),
/// built on the vectors and with the neighbour lists that `spares` keeps
/// (`XAaSpares::x_build_ahead`, where `XAaSpares::x_applies`); `None` where no table was made (the
/// block then builds its own).
fn blastp_aa_lookup_ahead(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    spares: &crate::algorithm::tblastx::lookup::XAaSpares,
) -> Option<(BlastpRuntimeLookup, Vec<QueryContext>)> {
    // `-word_size 4` (item SC8b): no table ahead; the block builds its own on the pool's threads
    // (`blastp_aa_lookup`).
    if !crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(args.word_size) {
        return None;
    }
    // EXPERIMENT (LOSAT_X_BPLK16): the table of 16-bit hits where the switch is on (a block of the
    // block loop is never a query chunk).
    if x_bplk16() {
        let threshold = crate::core::blast_util::ncbi_int4_from_double(args.threshold);
        if x_bplk16_for(args, false) {
            let (built, contexts) = spares.x_build_ahead16(
                query_frames,
                threshold,
                false,
                args.scoring.matrix,
                args.word_size,
            )?;
            return Some((x_bp_runtime_lookup(built), contexts));
        }
        let (lookup, contexts) = spares.x_build_ahead(
            query_frames,
            threshold,
            &lookup_protein_params_ungapped(args.scoring.matrix),
            false,
            args.scoring.matrix,
            args.word_size,
        )?;
        return Some((x_bp_runtime_lookup(XBpLookup::Wide(lookup)), contexts));
    }
    let (lookup, contexts) = spares.x_build_ahead(
        query_frames,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        &lookup_protein_params_ungapped(args.scoring.matrix),
        false,
        args.scoring.matrix,
        args.word_size,
    )?;
    Some((BlastpRuntimeLookup::Aa(lookup), contexts))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:245
/// EXPERIMENT (LOSAT_X_BPLKREUSE): `blastp_aa_lookup`'s table and contexts (the same arguments),
/// built on the vectors of the block before that `spare` keeps
/// (`XAaSpares::x_build_on_block_spare`); a table of 4-letter words (item SC8b) is
/// `blastp_aa_lookup`'s, on `pool`'s threads.
fn blastp_aa_lookup_on_block_spare(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    spare: &crate::algorithm::tblastx::lookup::XAaSpares,
    x_chunk: bool,
    pool: &crate::utils::threading::SearchPool<'_>,
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    if !crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(args.word_size) {
        return blastp_aa_lookup(args, query_frames, x_chunk, pool);
    }
    // EXPERIMENT (LOSAT_X_BPLK16): the table of 16-bit hits where the switch is on.
    if x_bplk16() {
        let threshold = crate::core::blast_util::ncbi_int4_from_double(args.threshold);
        if x_bplk16_for(args, x_chunk) {
            let (built, contexts) = spare.x_build_on_block_spare16(
                query_frames,
                threshold,
                false,
                args.scoring.matrix,
                args.word_size,
            );
            return (x_bp_runtime_lookup(built), contexts);
        }
        let (lookup, contexts) = spare.x_build_on_block_spare(
            query_frames,
            threshold,
            &lookup_protein_params_ungapped(args.scoring.matrix),
            false,
            args.scoring.matrix,
            args.word_size,
        );
        return (x_bp_runtime_lookup(XBpLookup::Wide(lookup)), contexts);
    }
    let (lookup, contexts) = spare.x_build_on_block_spare(
        query_frames,
        crate::core::blast_util::ncbi_int4_from_double(args.threshold),
        &lookup_protein_params_ungapped(args.scoring.matrix),
        false,
        args.scoring.matrix,
        args.word_size,
    );
    (BlastpRuntimeLookup::Aa(lookup), contexts)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1349-1354
/// The compressed lookup table of blastp-fast (`CompressedAaLookupTable`) and the query contexts of
/// `query_frames`, as `call_lookup` (`blast_engine/call_setup.rs`) builds them: `query_batches` and
/// `query_lengths` are the search's (its batches of the whole input, its records' lengths), and
/// `x_for` says which build of a block search this is (the spare vectors it takes,
/// `x_query_blocks::XLookupFor`; a query chunk's NCBI tables are its one table). The build of the
/// next block's table ahead of that block (LOSAT_X_BPLKPAHEAD, `x_lookup_ahead.rs`) calls it too,
/// with `XLookupFor::Ahead`. `x_appends` is the most appends of the build where the block search
/// counted them (`x_query_blocks::x_block_appends`), else 0.
#[allow(clippy::too_many_arguments)]
fn blastp_compressed_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    query_batches: &[seq_range::RangedBatch],
    query_lengths: &[usize],
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
    x_for: super::x_query_blocks::XLookupFor,
    x_appends: usize,
) -> Result<(BlastpRuntimeLookup, Vec<QueryContext>)> {
    // Another matrix than BLOSUM62: `blastp_matrix_compressed_lookup`.
    if args.scoring.matrix != ScoringMatrix::Blosum62 {
        return blastp_matrix_compressed_lookup(
            args,
            query_frames,
            query_batches,
            query_lengths,
            blastp_parallel_pool,
            x_for,
            x_appends,
        );
    }
    let query_chunk = x_for == super::x_query_blocks::XLookupFor::QueryChunk;
    // EXPERIMENT (LOSAT_X_BPQBLOCK): the contexts' Karlin blocks on the search's pool
    // (`x_query_blocks::x_bp_lookup_query`); without the switch
    // `prepare_blosum62_lookup_query_for_word_size`. EXPERIMENT (LOSAT_X_BPLKPFLAT): they
    // can be left to the build below and set after it (`x_query_blocks::XLookupFlat`); a
    // query chunk takes no spare table vectors.
    let (x_flat, query, locations, mut contexts) = super::x_query_blocks::XLookupFlat::lookup_query(
        query_frames,
        args.word_size,
        ScoringMatrix::Blosum62,
        blastp_parallel_pool,
        x_for,
    );
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1292
    // One table for all query batches; NCBI's tables (`ncbi_compressed_lookup_tables`)
    // only bound the overflow cells. With LOSAT_X_BPQBLOCK the queries are a block of
    // whole batches, and `query_batches` are the block's batches of the whole input, as
    // ranges of the block's lookup segments.
    let ncbi_tables = ncbi_compressed_lookup_tables(
        query_batches,
        query_lengths,
        &contexts,
        &locations,
        query_chunk,
        !args.ungapped,
    )?;
    // EXPERIMENT (LOSAT_X_BPCLUT; BPCLUT is a strict-set switch): the same table built
    // flat, each cell's hits in one slice in NCBI's scan order, with the same refusal
    // (lookup/x_clut_direct.rs), by the search's pool when it has more than one thread
    // (lookup/x_clut_par.rs); without the switch `try_build_blosum62_compressed_lookup`.
    let lookup = x_bp_build_compressed_lookup(
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(100.0 * args.threshold),
        &query,
        &locations,
        &ncbi_tables,
        blastp_parallel_pool,
        x_flat.build_args(x_appends),
    )
    .context("failed to build NCBI compressed amino-acid lookup table")?
    .map_err(|()| {
        anyhow::anyhow!(
            "-threshold {} with the compressed lookup table of -word_size {} is not supported by LOSAT's BLASTP: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
            args.threshold,
            args.word_size
        )
    })?;
    let lookup = match lookup {
        XCompressedLookup::Reference(lookup) => {
            // The flat build declined (a backbone past its stripes): the reference builder ran.
            if x_bpclut_on() {
                crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", false);
            }
            BlastpRuntimeLookup::Compressed(lookup)
        }
        XCompressedLookup::Flat(lookup) => {
            crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", true);
            BlastpRuntimeLookup::CompressedFlat(lookup)
        }
    };
    x_flat.finish(&mut contexts);
    Ok((lookup, contexts))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_aalookup.c:1303-1310
/// [`blastp_compressed_lookup`] of a matrix other than BLOSUM62 (items SC2b, SC2d-2): the contexts
/// with the Karlin blocks of the matrix's scores (`XLookupFlat::lookup_query`:
/// `prepare_matrix_lookup_query_for_word_size`, on the search's pool where it has more than one
/// thread), NCBI's tables of the query batches (`ncbi_compressed_lookup_tables`) and the compressed
/// table of the matrix's words (`x_bp_build_compressed_lookup_of`: the flat table, by the search's
/// pool where it has more than one thread, on the spare vectors of `x_for` as BLOSUM62's; the
/// reference one where the flat builders decline or with `LOSAT_X_BPCLUT=0`), with the refusals of
/// the BLOSUM62 table. `x_appends` goes to the build as BLOSUM62's (only BLOSUM62's serial builder
/// reserves it, `x_build_flat_in_pool`).
#[inline(never)]
fn blastp_matrix_compressed_lookup(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    query_batches: &[seq_range::RangedBatch],
    query_lengths: &[usize],
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
    x_for: super::x_query_blocks::XLookupFor,
    x_appends: usize,
) -> Result<(BlastpRuntimeLookup, Vec<QueryContext>)> {
    let query_chunk = x_for == super::x_query_blocks::XLookupFor::QueryChunk;
    let (x_flat, query, locations, mut contexts) = super::x_query_blocks::XLookupFlat::lookup_query(
        query_frames,
        args.word_size,
        args.scoring.matrix,
        blastp_parallel_pool,
        x_for,
    );
    let ncbi_tables = ncbi_compressed_lookup_tables(
        query_batches,
        query_lengths,
        &contexts,
        &locations,
        query_chunk,
        !args.ungapped,
    )?;
    let lookup = x_bp_build_compressed_lookup_of(
        args.scoring.matrix,
        args.word_size,
        crate::core::blast_util::ncbi_int4_from_double(100.0 * args.threshold),
        &query,
        &locations,
        &ncbi_tables,
        blastp_parallel_pool,
        x_flat.build_args(x_appends),
    )
    .ok_or_else(|| {
        anyhow::anyhow!(
            "-word_size {} (the compressed lookup table) is not supported by LOSAT's BLASTP",
            args.word_size
        )
    })?
    .map_err(|()| {
        anyhow::anyhow!(
            "-threshold {} with the compressed lookup table of -word_size {} is not supported by LOSAT's BLASTP: the neighboring words of a query batch fill NCBI BLAST+'s 1024 overflow banks, past which NCBI BLAST+ 2.17.0 writes and corrupts its heap",
            args.threshold,
            args.word_size
        )
    })?;
    let lookup = match lookup {
        XCompressedLookup::Reference(lookup) => {
            // The flat builders declined: the reference builder ran.
            if x_bpclut_on() {
                crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", false);
            }
            BlastpRuntimeLookup::Compressed(lookup)
        }
        XCompressedLookup::Flat(lookup) => {
            crate::utils::xreach::x_reach!("LOSAT_X_BPCLUT", true);
            BlastpRuntimeLookup::CompressedFlat(lookup)
        }
    };
    x_flat.finish(&mut contexts);
    Ok((lookup, contexts))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:2771-2803 (cited at
// `x_context_karlin` in lookup/backbone.rs)
/// EXPERIMENT (LOSAT_X_BPQBLOCK): the contexts of the search of all queries whose preliminary stage
/// runs in blocks (`run_call`'s `x_blocks`), as each table's builder prepares them
/// (`prepare_matrix_lookup_query_for_word_size` of the search's matrix and word size; their Karlin
/// blocks on the pool's threads), and an empty table: each block builds its own (unit test
/// `x_blocks_contexts_are_the_builders_contexts`).
fn x_blocks_contexts(
    args: &ResolvedBlastpArgs,
    query_frames: &[Vec<QueryFrame>],
    pool: &crate::utils::threading::SearchPool<'_>,
) -> (BlastpRuntimeLookup, Vec<QueryContext>) {
    let (_, _, contexts) =
        crate::algorithm::tblastx::lookup::prepare_matrix_lookup_query_for_word_size(
            query_frames,
            args.word_size,
            false,
            args.scoring.matrix,
            Some(pool),
        );
    (empty_aa_lookup(args.word_size), contexts)
}

// No NCBI counterpart: a table that no scan reads; it changes no value NCBI computes.
/// An empty table of exact blastp's words of `word_size`, for a call that scans with another
/// table: rpsblast's of the database (`rps_hooks::rps_contexts`), or each block's of the
/// preliminary stage in blocks (`x_blocks_contexts`).
fn empty_aa_lookup(word_size: usize) -> BlastpRuntimeLookup {
    BlastpRuntimeLookup::Aa(crate::algorithm::tblastx::lookup::BlastAaLookupTable {
        backbone: Vec::new(),
        overflow: Vec::new(),
        pv: Vec::new(),
        frame_bases: Vec::new(),
        num_contexts: 0,
        query_length: 0,
        word_length: word_size as i32,
        alphabet_size: crate::core::blast_util::BLASTAA_SIZE as i32,
        charsize: 5,
        mask: 0,
        longest_chain: 0,
        threshold: 0,
        row_max: Vec::new(),
    })
}

// No NCBI counterpart: memory only (which vectors a later build may take).
/// Whether a later block's build may take the vectors of exact blastp's `table`: those of the
/// direct builders' word sizes (2 and 3); the joint builder of a table of 4-letter words (item
/// SC8b, `blastp_aa_lookup`) takes none, so those are freed with the table.
fn x_bp_aa_vectors_reused(table: &crate::algorithm::tblastx::lookup::BlastAaLookupTable) -> bool {
    crate::algorithm::tblastx::lookup::x_pair_lut::x_direct_covers(table.word_length as usize)
}

// No NCBI counterpart: memory only; it changes no value NCBI computes.
/// EXPERIMENT (LOSAT_X_BPLKDROP): a strict-set switch, on unless it is `0`. The search's lookup
/// table after its preliminary stage, which nothing reads any more: kept (`Some`) only where the
/// search gives its vectors to a spare at its end that a later block's build may take
/// (`x_query_blocks::x_bp_lookup_spare` / `x_bp_aa_spares` / `x_bp_aa_block_spare` and
/// `XAheadBlock::x_table_taken_later`), or with the switch off; otherwise dropped here, before the
/// redo, instead of at the search's end. Item SC2d: first, where the next block's flat table is built
/// ahead, the table's vectors go to that build, which goes on beside the rest of this block's search
/// (`XAheadBlock::after_prelim`; `None`).
fn x_lookup_after_prelim(
    lookup: BlastpRuntimeLookup,
    pool: &crate::utils::threading::SearchPool<'_>,
    ahead: &x_lookup_ahead::XAheadBlock<'_>,
) -> Option<BlastpRuntimeLookup> {
    let lookup = ahead.after_prelim(lookup)?;
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    if !*ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPLKDROP")) {
        return Some(lookup);
    }
    let taken_later = ahead.x_table_taken_later()
        && match &lookup {
            BlastpRuntimeLookup::CompressedFlat(_) => {
                super::x_query_blocks::x_bp_lookup_spare(pool).is_some()
            }
            BlastpRuntimeLookup::Aa(table) if !x_bp_aa_vectors_reused(table) => false,
            BlastpRuntimeLookup::Aa(_) | BlastpRuntimeLookup::Aa16(_) => {
                super::x_query_blocks::x_bp_aa_spares(pool).is_some()
                    || super::x_query_blocks::x_bp_aa_block_spare(pool).is_some()
            }
            BlastpRuntimeLookup::Compressed(_) => false,
        };
    crate::utils::xreach::x_reach!("LOSAT_X_BPLKDROP", !taken_later);
    taken_later.then_some(lookup)
}

/// A call of the preliminary stage only, which writes no report and hands its hit lists to its
/// caller: a query chunk of a split batch (`search_split_query_batches`), whose lists leave as
/// lists of full HSPs, or (EXPERIMENT LOSAT_X_BPQBLOCK) a block of whole batches of the search of
/// all queries (`x_prelim_block_plan`), whose lists leave compact, numbered from the block's first
/// query.
enum PrelimOnly<'a> {
    Chunk(&'a mut Vec<Option<BlastpHitList>>),
    /// The block's lists, and what the call of all queries made once for all its blocks: its
    /// subjects, encoded, with `-lcase_masking` their scan ranges (`x_soft_subject`), and the
    /// block's query frames (SEG made once); `first_query`: the index among the call's queries of
    /// the block's first (the block counts its queries from it; the diagnostics name the call's).
    Block {
        lists: &'a mut Vec<Option<BlastpHitList<XPrelimHsp>>>,
        subjects: &'a [EncodedProtein],
        soft_subjects: Option<&'a [Vec<(i32, i32)>]>,
        frames: &'a [Vec<QueryFrame>],
        first_query: usize,
    },
}

impl PrelimOnly<'_> {
    /// Whether the call is a query chunk's (its table is of one chunk: no 16-bit hits, no spare
    /// vectors, `XLookupFor::QueryChunk`; its batch is searched).
    fn is_chunk(&self) -> bool {
        matches!(self, PrelimOnly::Chunk(_))
    }

    /// Hands the call's hit lists to its caller where its preliminary stage ends (`run_call`).
    fn hand_off(
        self,
        mut preliminary_hit_lists: Vec<Option<BlastpHitList>>,
        mut x_compact_hit_lists: Vec<Option<BlastpHitList<XPrelimHsp>>>,
    ) {
        match self {
            PrelimOnly::Chunk(chunk_lists) => {
                // A chunk's lists leave as `BlastpHitList`s of full HSPs, so the chunk merge
                // (`query_split.rs`) reads what it reads from a search without chunks.
                x_finish_compact_lists(&mut preliminary_hit_lists, &mut x_compact_hit_lists);
                *chunk_lists = preliminary_hit_lists;
            }
            // EXPERIMENT (LOSAT_X_BPQBLOCK): a block's lists leave compact, to the search of all
            // queries, which moves them to their queries (`x_blocks`); the block's table and
            // scratch are freed when `run_call` returns right after this.
            PrelimOnly::Block { lists, .. } => *lists = x_compact_hit_lists,
        }
    }
}

/// A call of the batches of `query_records` (`run_call`), and, with `-export_search_strategy`, the
/// search strategy's requests of its batches, released where the call ends (`CallRequests`): one
/// place for every way the call ends. A call of the preliminary stage only (`prelim`) has none.
#[allow(clippy::too_many_arguments)]
fn run_resolved_in_pool(
    args: ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    subject_source: BlastpSubjectSource<'_>,
    query_label: &str,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
    prelim: Option<PrelimOnly<'_>>,
    ranges: &BlastpRanges<'_>,
    reports: &BlastpReports<'_>,
    x_ahead: x_lookup_ahead::XAheadBlock<'_>,
) -> Result<()> {
    let call = match (args.strategy.clone(), prelim.is_none()) {
        (Some(strategy), true) => CallRequests::of_call(
            strategy,
            &args,
            query_records,
            subject_source,
            ranges,
            reports.strategy_next,
        ),
        _ => None,
    };
    let reports = BlastpReports {
        call_requests: call.as_ref(),
        ..*reports
    };
    let result = run_call(
        args,
        query_records,
        subject_source,
        query_label,
        subject_label,
        outputs,
        blastp_parallel_pool,
        prelim,
        ranges,
        &reports,
        x_ahead,
    );
    match &call {
        Some(call) => call.release(&result).and(result),
        None => result,
    }
}

/// An ungapped search's failure before the search of its first query batch
/// (`ungapped_setup_failure`).
enum UngappedSetupFailure {
    /// NCBI's error of the first batch (`preliminary_search_error` writes it after the outfmt 0
    /// prolog), or LOSAT's refusal of the matrix file that NCBI reads, in its place.
    Engine(String),
    /// LOSAT's refusal before any output.
    Refusal(anyhow::Error),
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_stat.c:1679-1680 (the score block of
/// the first query batch: a `-matrix` name that the matrix loader does not know fails with -1,
/// posted for all queries, blast_setup.c:498-501; `blastx::query_setup::unknown_matrix_message`)
/// An ungapped search takes any `-matrix` name (a gapped one stops at NCBI's matrix error with
/// the options): its first query batch fails to set up, after the outfmt 0 prolog, with NCBI's
/// error of each query of the batch (`BLAST engine error: Error: Unknown error code -1 `, exit 3),
/// or LOSAT's refusal of the file that NCBI reads for the name (`matrix_file`).
///
/// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
/// NCBI searches a longer subject in overlapping parts and merges the parts' HSP lists
/// (`Blast_HSPListsMerge`, blast_engine.c:578-584); an ungapped search does not port the merge
/// (BLASTX refuses the same subject, `blastx::database::too_long`). The refusal belongs to the
/// first batch with a valid context, as the reader's errors do.
fn ungapped_setup_failure(
    call_in: &call_setup::CallSetup<'_, '_>,
    query_batches: &[seq_range::RangedBatch],
    query_valid: &[bool],
    subject_source: BlastpSubjectSource<'_>,
) -> Option<UngappedSetupFailure> {
    let args = call_in.args;
    if args.unknown_matrix {
        if let Some(first) = query_batches
            .iter()
            .find(|batch| !batch.searched.is_empty())
        {
            return Some(UngappedSetupFailure::Engine(match &args.matrix_file {
                Some(file) => crate::blastinput::ncbi_environment::matrix_file_refusal(
                    &args.matrix_name,
                    file,
                    "BLASTP",
                ),
                None => unknown_matrix_message(&call_in.query_records[first.searched.clone()]),
            }));
        }
    }
    if args.ungapped && query_valid.iter().any(|valid| *valid) {
        let longest = match subject_source {
            BlastpSubjectSource::Fasta(records) => records
                .iter()
                .enumerate()
                .find(|(_, record)| {
                    record.sequence.len() > crate::core::subject_split::MAX_DBSEQ_LEN
                })
                .map(|(index, record)| {
                    (
                        format!("subject sequence {}", index + 1),
                        record.sequence.len(),
                    )
                }),
            BlastpSubjectSource::Database(database) => {
                let longest = database.seqdb().max_length() as usize;
                (longest > crate::core::subject_split::MAX_DBSEQ_LEN)
                    .then(|| ("a database sequence".to_string(), longest))
            }
        };
        if let Some((what, length)) = longest {
            return Some(UngappedSetupFailure::Refusal(too_long_subject_refusal(
                &what, length,
            )));
        }
    }
    None
}

/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:82-89
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:123-125
/// The statistics of each query: what the XML reports print and what the observer receives.
/// A PSSM's search (psiblast's later rounds, `-in_pssm`) reads the PSSM's gapped block there,
/// as its E-values and the pairwise report's `Gapped` line do (`PsiSearch::gapped_karlin`).
fn query_xml_stats(
    args: &ResolvedBlastpArgs,
    gapped_params: KarlinParams,
    contexts: &[QueryContext],
    query_valid: &[bool],
    search_spaces: &[SearchSpace],
    exact_search_spaces: &[Option<i64>],
    rps_gapped_k: Option<&[f64]>,
    queries: usize,
) -> Vec<crate::report::blast_xml::XmlStats> {
    let gapped = args
        .psi_search
        .as_deref()
        .map_or(gapped_params, |psi| psi.gapped_karlin(&gapped_params));
    (0..queries)
        .map(|q_idx| {
            if !query_valid[q_idx] {
                return crate::report::blast_xml::XmlStats::NONE;
            }
            // NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_results.cpp:97-103
            // An ungapped search has no `kbp_gap`: the XML statistics are the query's own
            // `kbp_std` (probe 1: kappa 0.12756, lambda 0.314758 of query 1).
            let block = call_setup::context_block(args, &gapped, &contexts[q_idx]);
            crate::report::blast_xml::XmlStats {
                hsp_len: search_spaces[q_idx].length_adjustment as i32,
                eff_space: exact_search_spaces[q_idx]
                    .unwrap_or(search_spaces[q_idx].effective_space.round() as i64),
                // rpsblast: the footer's K, 1.2 times the profile's of the batch's first query
                // (`rps_hooks::rps_search`); the block's for every other.
                kappa: rps_gapped_k.map_or(block.k, |k| k[q_idx]),
                lambda: block.lambda,
                entropy: block.h,
            }
        })
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn run_call(
    args: ResolvedBlastpArgs,
    query_records: &[FastaRecord],
    subject_source: BlastpSubjectSource<'_>,
    _query_label: &str,
    subject_label: &str,
    outputs: &mut ReportOutputs<'_>,
    blastp_parallel_pool: &crate::utils::threading::SearchPool<'_>,
    prelim: Option<PrelimOnly<'_>>,
    ranges: &BlastpRanges<'_>,
    reports: &BlastpReports<'_>,
    mut x_ahead: x_lookup_ahead::XAheadBlock<'_>,
) -> Result<()> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1681
    // The records of a FASTA subject set (none for a database, whose sequences are read
    // in blocks below and, for the hit lists, after the preliminary stage).
    let subject_records = subject_source.fasta_records();
    let database = subject_source.database();
    // A query chunk's search (`PrelimOnly::Chunk`).
    let query_chunk = prelim.as_ref().is_some_and(PrelimOnly::is_chunk);
    // EXPERIMENT (LOSAT_X_BPQBLOCK): a block of the preliminary stage of the search of all queries
    // (`PrelimOnly::Block`): the call's subjects, their scan ranges and the block's query frames,
    // from that call.
    let x_block_inputs = match &prelim {
        Some(PrelimOnly::Block {
            subjects,
            soft_subjects,
            frames,
            ..
        }) => Some((*subjects, *frames, *soft_subjects)),
        _ => None,
    };
    // Diagnostic (LOSAT_DEBUG_BLASTP_TRACE_HSP): the index among the call's queries of this call's
    // first query, which the preliminary stage's trace (`BlastpPrelimKernel::gapped_subject` in
    // `blast_engine/prelim_kernel.rs`) adds to its query index: a block of the preliminary stage
    // (`PrelimOnly::Block`) counts its queries from its first, and the trace names the queries of
    // the search of all queries, as its redo's does. No NCBI counterpart.
    let x_trace_query_base = match &prelim {
        Some(PrelimOnly::Block { first_query, .. }) => *first_query as u32,
        _ => 0,
    };
    let timing_enabled = blastp_timing_env_enabled();
    let kappa_range_counters_enabled = blastp_kappa_range_counters_env_enabled();
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1705
    reset_blastp_kappa_range_counters(kappa_range_counters_enabled);
    let t_total = blastp_timing_start(timing_enabled);
    let timing = if timing_enabled {
        Some(Arc::new(BlastpTiming::new()))
    } else {
        None
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
    //
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1744-1748
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let blastp_num_threads = blastp_parallel_pool.threads();
    let (requested_formats, needs) = call_setup::requested_reports(&args, outputs)?;

    let call_setup::BlastpQueryMeta {
        query_ids,
        query_titles,
        query_lengths,
    } = call_setup::query_meta(query_records);
    let call_in = call_setup::CallSetup {
        args: &args,
        query_records,
        subject_records,
        prelim: &prelim,
        query_chunk,
        x_block_inputs,
        blastp_parallel_pool,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
        ranges,
        reports,
    };

    let query_encoding_start = blastp_timing_start(timing_enabled);
    // EXPERIMENT (LOSAT_X_BPLKPAHEAD): a block's query frames, lookup table and contexts built ahead
    // of it, during the block before it (`blast_engine/x_lookup_ahead.rs`); else encoded by
    // `query_frames_of` and built by `call_lookup` (`blast_engine/call_setup.rs`).
    let mut x_prebuilt = x_ahead.take_prebuilt(query_records);
    let mut x_owned_frames: Vec<Vec<QueryFrame>> = Vec::new();
    let query_frames: &[Vec<QueryFrame>] =
        call_setup::query_frames_of(&call_in, &mut x_prebuilt, &mut x_owned_frames);
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.query_encoding_ns, query_encoding_start);
    }

    let subject_encoding_start = blastp_timing_start(timing_enabled);
    let mut x_owned_subjects: Vec<EncodedProtein> = Vec::new();
    let call_setup::BlastpCallSubjects {
        subjects,
        x_meta,
        subject_titles,
    } = call_setup::call_subjects(&call_in, &mut x_owned_subjects)?;
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.subject_encoding_ns, subject_encoding_start);
    }
    // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs): the
    // subjects' encoding (and their titles). No NCBI counterpart.
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Encode, subject_encoding_start);

    let lookup_start = blastp_timing_start(timing_enabled);
    let call_setup::BlastpScoring {
        ideal_ungapped_lambda,
        gapped_params,
        x_drop_gapped,
        x_drop_gapped_final,
    } = call_setup::call_scoring(&args)?;
    let call_setup::BlastpCallLookup {
        query_batches,
        x_blocks,
        lookup,
        contexts,
        context_lookup_bounds,
        concatenated_query,
        x_context_table,
        query_valid,
        query_batch_skipped,
    } = call_setup::call_lookup(
        &call_in,
        query_frames,
        &query_lengths,
        &mut x_prebuilt,
        &x_ahead,
    )?;
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_duration(&timing.lookup_ns, lookup_start);
    }
    // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs): the
    // lookup table and the query contexts. No NCBI counterpart.
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::BuildLookup, lookup_start);
    // Diagnostic (LOSAT_TIMING): the search, from here to the end of the final traceback, for the
    // common line (utils/xtiming.rs); a query chunk's search (`search_split_query_batches`)
    // marks its own. No NCBI counterpart.
    let t_search = crate::utils::xtiming::clock();

    // The setup of the first searched batch fails: its error follows the outfmt 0 prolog and the
    // reports of the batches before it (`preliminary_search_error`).
    let setup_error =
        |outputs: &mut ReportOutputs<'_>, first_batch: bool, database_name: &str, text: &str| {
            preliminary_search_error(
                &args,
                query_records,
                subject_source,
                subject_label,
                outputs,
                blastp_parallel_pool,
                ranges,
                reports,
                &query_batches,
                &query_batch_skipped,
                first_batch,
                database_name,
                text,
            )
        };
    // An ungapped search's failures before the search of its first query batch
    // (`ungapped_setup_failure`): NCBI's error of a `-matrix` that the score block cannot load,
    // and LOSAT's refusal of a subject over NCBI's part length.
    let failure = ungapped_setup_failure(&call_in, &query_batches, &query_valid, subject_source);
    if let Some(failure) = failure {
        return Err(match failure {
            UngappedSetupFailure::Engine(text) => setup_error(outputs, true, "", &text),
            UngappedSetupFailure::Refusal(error) => error,
        });
    }

    // A database without sequences (makeblastdb builds none; a volume whose index holds
    // no sequence) and without an alias file's LENGTH or STATS_TOTLEN
    // (`SeqDb::first_oid_read_fails`): NCBI opens it, writes the outfmt 0 prolog, reads the
    // first query batch and fails in that batch's preliminary search, where its reader cannot
    // map the sequence file (fixture `LOSAT/tests/seqdb/db/empty`).
    // A first batch of invalid queries only is not searched: the batches up to the first one that
    // is are reported before the error (`unsearched_batches_then_engine_error`); a call of invalid
    // queries only searches no batch, and NCBI reports them (NCBI 2.17.0 blastp, an all-X query
    // and the empty database: exit 0).
    if let Some(database) = database.filter(|database| {
        database.first_oid_read_fails()
            && (query_chunk || query_batch_skipped.iter().any(|skipped| !skipped))
    }) {
        return Err(setup_error(
            outputs,
            query_chunk,
            database.name(),
            "Database memory map file error",
        ));
    }

    // The database's extent (`call_setup::db_extent`) and NCBI's error of its setup's average
    // subject length (`call_setup::setup_fails`, after LOSAT's refusal of a count of 0).
    let extent = call_setup::db_extent(&args, database, subjects)?;
    // The setup of the first searched batch fails: NCBI has written the outfmt 0 prolog and the
    // batch's reader messages (`preliminary_search_error`; NCBI 2.17.0 blastp and psiblast
    // `-subject` of the letters `M`, ``, `` with `-searchsp` or `OLD_FSC`: exit 3). After a first
    // batch of invalid queries only, which NCBI does not search, the error follows the reports of
    // the batches before the first searched one (`unsearched_batches_then_engine_error`).
    if call_setup::setup_fails(&args, &extent, &contexts)? {
        return Err(setup_error(
            outputs,
            query_chunk,
            "",
            crate::stats::setup::SUBJECT_LENGTH_TOO_SHORT,
        ));
    }
    let call_setup::BlastpDbExtent {
        total_db_len,
        sorts_results,
        gapped_gumbel,
        num_subjects,
        eff_db_length,
        cutoff_subject_length,
        rps_units,
        ..
    } = extent;
    let call_statistics = call_setup::CallStatistics {
        args: &args,
        database,
        ranges,
        contexts: &contexts,
        gapped_params: &gapped_params,
        gapped_gumbel: &gapped_gumbel,
        total_db_len,
        num_subjects,
        eff_db_length,
        cutoff_subject_length,
    };
    let call_setup::BlastpSearchSpaces {
        search_spaces,
        exact_search_spaces,
        per_subject,
        option_space_of,
        preliminary_option_spaces,
        preliminary_spaces,
    } = call_setup::call_search_spaces(&call_statistics)?;
    let (preliminary_search_spaces, preliminary_exact_spaces) = match &preliminary_spaces {
        Some((spaces, exact)) => (spaces, exact),
        None => (&search_spaces, &exact_search_spaces),
    };
    let call_setup::BlastpCallCutoffs {
        hit_cutoff_scores,
        word_cutoff_scores,
        sum_stats,
        hit_cutoff_score_min,
        prelim_hitlist_size,
        prelim_evalue,
    } = call_setup::prelim_cutoffs(&call_statistics, preliminary_exact_spaces)?;
    // The statistics of each subject (`per_subject`), which `prelim_run::run_preliminary_with` reads.
    let subject_parameters = call_setup::subject_parameters(
        &call_statistics,
        &preliminary_option_spaces,
        preliminary_search_spaces,
        preliminary_exact_spaces,
    );
    let call_setup::BlastpPrelimSinks {
        mut prelim_filter,
        subject_space_overflow,
        mut x_prelim_freed,
        mut redone_matches,
        mut preliminary_hit_lists,
        x_expanded_at_start,
        x_pooled_at_start,
        mut x_compact_hit_lists,
    } = call_setup::prelim_sinks(
        &args,
        x_blocks.as_deref(),
        prelim_hitlist_size,
        query_ids.len(),
    );

    // The scan's sizes (`call_setup::scan_shape`; the owner of the scan ranges it makes).
    let mut x_owned_soft_subjects: Option<Vec<Vec<(i32, i32)>>> = None;
    let call_setup::BlastpScanShape {
        x_dropoff_per_context,
        window,
        word_size,
        x_ungap_sent,
        offset_array_size,
        diag_array_size,
        diag_mask,
        total_query_span,
        x_soft_subjects,
    } = call_setup::scan_shape(
        &call_in,
        &lookup,
        &contexts,
        &concatenated_query,
        x_blocks.as_deref(),
        &mut x_owned_soft_subjects,
    )?;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1498
    // The kernel takes the subject from its caller: a record of the encoded FASTA set,
    // or a sequence of a database block (`BlastpSubjectRef`); `s_idx` is the subject's
    // index in the set, or its ordinal id in the database.
    // The kernel runs no rayon work: with `-db` the window of `database/x_db_keep.rs` makes the
    // workers wait for the one that runs it (header: why the waits end).
    let batch_values = BlastpPrelimValues {
        word_cutoff_scores: &word_cutoff_scores,
        hit_cutoff_scores: &hit_cutoff_scores,
        search_spaces: preliminary_search_spaces,
        prelim_evalue,
        sum_stats: sum_stats.as_ref(),
    };
    // The kernel of one subject (`prelim_kernel::BlastpPrelimKernel::search_subject`: its scan,
    // ungapped extensions and gapped half) with what it reads of the call; the preliminary stage
    // calls it with the subject's preliminary statistics `values` (`prelim_run::
    // run_preliminary_with`), alone or inside the masked kernels (`db_masked.rs`).
    let kernel = prelim_kernel::BlastpPrelimKernel {
        args: &args,
        query_ids: &query_ids,
        contexts: &contexts,
        context_lookup_bounds,
        concatenated_query: &concatenated_query,
        lookup: &lookup,
        x_context_table: &x_context_table,
        x_dropoff_per_context: &x_dropoff_per_context,
        window,
        word_size,
        x_ungap_sent: &x_ungap_sent,
        diag_mask,
        total_query_span,
        x_soft_subjects,
        x_drop_gapped,
        gapped_params: &gapped_params,
        gapped_gumbel: &gapped_gumbel,
        x_trace_query_base,
        timing_enabled,
        timing: &timing,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1498
    // The preliminary stage (`blast_engine/prelim_run.rs`): the choice of the pool and the call's
    // branch, into the hit lists and the culling writer; RPS-BLAST's gapped K and its volumes'
    // merge are read below.
    // `-db_soft_mask` and `-db_hard_mask` (`db_masked.rs`): a database search with either searches
    // its sequences with the masked kernels below, chosen once here.
    let db_masking = database
        .zip(SearchMask::of(args.db_mask_resolved.as_ref()))
        .map(|(database, mask)| db_masked::BlastpDbMasking {
            db: database.seqdb(),
            mask,
            program: args.app.name(),
            error: std::sync::Mutex::new(None),
        });
    let prelim_stage = prelim_run::PrelimRun {
        args: &args,
        query_records,
        subject_source,
        subject_label,
        blastp_parallel_pool,
        ranges,
        reports,
        query_batches: &query_batches,
        query_batch_skipped: &query_batch_skipped,
        query_chunk,
        rps_recount_identities: recounts_rps_identities(&requested_formats),
        database,
        subjects,
        x_soft_subjects,
        query_frames,
        x_blocks: &x_blocks,
        x_block_call: x_block_inputs.is_some(),
        lookup: &lookup,
        contexts: &contexts,
        gapped_params: &gapped_params,
        gapped_gumbel: &gapped_gumbel,
        exact_search_spaces: &exact_search_spaces,
        preliminary_exact_spaces,
        preliminary_search_spaces,
        preliminary_option_spaces: &preliminary_option_spaces,
        per_subject,
        num_subjects,
        cutoff_subject_length,
        rps_units,
        prelim_hitlist_size,
        prelim_evalue,
        subject_space_overflow: &subject_space_overflow,
        db_masking: db_masking.as_ref(),
        window,
        diag_mask,
        offset_array_size,
        diag_array_size,
        total_query_span,
        timing_enabled,
        timing: &timing,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
    };
    let prelim_run::PrelimOutcome {
        rps_gapped_k,
        rps_multi_volume,
        rps_unlisted,
    } = match db_masking.as_ref() {
        None => prelim_run::run_preliminary_with(
            &prelim_stage,
            &|scratch: &mut BlastpSubjectScratch,
              s_idx: usize,
              subject: BlastpSubjectRef<'_>,
              subject_result: &mut BlastpSubjectResult,
              values: &BlastpPrelimValues<'_>| {
                kernel.search_subject(scratch, s_idx, subject, subject_result, values)
            },
            &batch_values,
            &subject_parameters,
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            &mut prelim_filter,
            outputs,
        )?,
        Some(masking) => prelim_run::run_preliminary_with(
            &prelim_stage,
            &|scratch: &mut BlastpSubjectScratch,
              s_idx: usize,
              subject: BlastpSubjectRef<'_>,
              subject_result: &mut BlastpSubjectResult,
              values: &BlastpPrelimValues<'_>| {
                masking.search_subject(&kernel, scratch, s_idx, subject, subject_result, values)
            },
            &batch_values,
            &subject_parameters,
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            &mut prelim_filter,
            outputs,
        )?,
    };
    let no_alignment_set = no_alignment_sets(&query_batch_skipped, rps_multi_volume);

    // `-culling_limit`: `s_BlastHSPCullingFinal` makes the hit lists of the queries from the
    // trees (hspfilter_culling.c:516-581; `culling.rs`), where the collector's lists were
    // updated as the subjects came; `-best_hit_*`: `s_BlastHSPBestHitFinal` from the best hit
    // lists (hspfilter_besthit.c:194-231; `besthit.rs`). A query chunk's lists leave below as the
    // collector's do.
    if let Some(filter) = prelim_filter.take() {
        filter.finish(&mut x_compact_hit_lists);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:268-271
    // A query chunk's preliminary stage ends here, before its stream is closed: its hit
    // lists merge into the batch's (`search_split_query_batches`).
    if let Some(prelim) = prelim {
        prelim.hand_off(preliminary_hit_lists, x_compact_hit_lists);
        crate::utils::xtiming::add(crate::utils::xtiming::Phase::Search, t_search);
        return Ok(());
    }
    // Diagnostic (LOSAT_TIMING): the search stops here and goes on after the chunks' searches.
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Search, t_search);
    // The search strategy's requests carry each batch's `EffectiveSearchSpace` option, the
    // archive's (`FirstContext::context0_space`), now that the call's set-up has computed them.
    if let Some(call) = reports.call_requests {
        set_call_searched(
            call,
            &args,
            rps_units,
            &contexts,
            &gapped_params,
            (&query_valid, &query_batch_skipped),
            (&exact_search_spaces, &search_spaces),
        )?;
    }
    search_split_query_batches(
        &args,
        query_records,
        ranges.input_records,
        subject_source,
        subject_label,
        outputs,
        blastp_parallel_pool,
        &query_lengths,
        &|query| Ok(option_space_of(&contexts[query])?),
        ranges.query_input,
        &mut preliminary_hit_lists,
        &mut x_compact_hit_lists,
        reports.x_prolog,
        reports.call_requests,
    )?;
    prelim_run::close_stream(&args, &mut preliminary_hit_lists, &mut x_compact_hit_lists);
    let t_search = crate::utils::xtiming::clock();
    prelim_run::sort_preliminary_lists(
        &args,
        &mut preliminary_hit_lists,
        &mut x_compact_hit_lists,
        timing.as_deref(),
    );

    // A database: the sequences of the preliminary hit lists are read now, by ordinal id,
    // and stand for the subject set in the redo; the hit lists name them by their rank
    // among the listed ordinal ids from here on (`database::ListedSubjects`). Their
    // definition lines are read after the redo (`ListedSubjects::read_titles`).
    // The listed subjects have no copy of their residues (`ListedSubjects::x_list`); the redo
    // reads each match's subject from the database.
    let mut listed_subjects = match database {
        Some(database) => Some(database::ListedSubjects::x_list(
            database,
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            !args.ungapped,
        )?),
        None => None,
    };
    // The subject set of a `-subject` search, which the reports read after the redo.
    let subject_set: &[EncodedProtein] = subjects;
    let subjects = match (listed_subjects.as_ref(), database) {
        // `-db_hard_mask`: the redo and the traceback read the subjects with `X` over their masks.
        (Some(listed), Some(database)) => listed.redo_subjects(
            database,
            SearchMask::of(args.db_mask_resolved.as_ref())
                .filter(|mask| mask.hard())
                .map(|mask| mask.algo_id),
        ),
        _ => super::x_redo_subjects::XRedoSubjects::Encoded(subject_set),
    };
    // The redo reads lists of full HSPs. With more than one thread its parallel paths read the
    // lists of all queries, so all are expanded here (but for the windows of the redo below); the
    // one-thread loop (`redo_run::one_thread_redo`) expands each query's list when it reaches it.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_compact_expand_all = blastp_num_threads > 1;
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let x_compact_expand_all = false;
    // Where the redo in windows (LOSAT_X_BPREDOPAR) may run,
    // its windows expand the compact lists as they reach them (`x_redo_par.rs`); where it does not
    // run, every list is expanded in `redo_run::plan_redo_and_expand_lists` before the other paths.
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    let x_win_expand = x_compact_expand_all && super::x_redo_par::x_bp_redo_par();
    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    let x_win_expand = false;
    // `-comp_based_stats 0`: the ordinary traceback below reads the compact lists as they are.
    // psiblast `-in_pssm`: the PSSM's traceback (`pssm_hooks::pssm_traceback`, mode 0 or the
    // position-based redo) takes the place of both, and reads them as they are too.
    let ordinary_traceback = args.comp_based_stats.mode
        == BlastpCompositionMode::NoCompositionBasedStats
        || args.psi_search.is_some()
        || args.rps_search.is_some();
    if x_compact_expand_all && !x_win_expand && !ordinary_traceback {
        x_finish_compact_lists(&mut preliminary_hit_lists, &mut x_compact_hit_lists);
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3329-3334
    // EXPERIMENT (LOSAT_X_BPIHPOOL): no subject scan runs any more; the pool's arrays of initial
    // hits are freed, as the subject loop's scratches freed theirs (item SC2d: before the next
    // block's build goes on below).
    x_ih_pool::release();
    // EXPERIMENT (LOSAT_X_BPLKDROP): nothing reads the lookup table after the preliminary stage;
    // it is dropped here unless a later block's build may take its vectors (at the end below).
    // Item SC2d: a flat table's vectors go to the next block's build here.
    let lookup = x_lookup_after_prelim(lookup, blastp_parallel_pool, &x_ahead);
    let kappa_redo_start = blastp_timing_start(timing_enabled);
    let fetch_notes = database::FetchNotes::new(listed_subjects.as_ref(), database);
    // `-comp_based_stats 0`: the ordinary traceback (blast_traceback.c:1500-1723) takes every
    // preliminary list here, so the composition-based redo below finds none; its hit lists
    // replace the redo's after it.
    let ordinary_hit_lists = if args.rps_search.is_some() {
        // rpsblast: the RPS engine's final hit lists (`rps_hooks::rps_search`).
        Some(preliminary_hit_lists.iter_mut().map(Option::take).collect())
    } else if ordinary_traceback && args.psi_search.is_some() {
        // psiblast `-in_pssm`: the PSSM's traceback (`pssm_hooks::pssm_traceback`: mode 0's, or
        // the position-based redo of modes 1 to 3, `FetchNotes::pssm_redo`).
        Some(pssm_hooks::pssm_traceback(
            &args,
            &contexts,
            total_query_span,
            subjects,
            &gapped_params,
            gapped_gumbel.as_ref(),
            &exact_search_spaces,
            &search_spaces,
            (cutoff_subject_length, prelim_evalue),
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            blastp_parallel_pool,
            &|heaps: &[_]| fetch_notes.pssm_redo(heaps),
            sorts_results,
        )?)
    } else if args.ungapped {
        // `-ungapped`: NCBI skips the traceback (`perform_traceback` is the gapped mode,
        // blast_traceback.c:1526); the lists stay as the preliminary stage made them.
        redo_run::sw_ungapped_check(args.use_sw_tback, &query_batches, &query_batch_skipped)?;
        let finish = traceback::UngappedFinish {
            contexts: &contexts,
            subjects,
            matrix: args.scoring.matrix,
            hitlist_size: args.traceback_hitlist_size(),
            sorts_results,
        };
        Some(finish.run(&mut preliminary_hit_lists, &mut x_compact_hit_lists)?)
    } else if ordinary_traceback {
        let ordinary = traceback::OrdinaryTraceback {
            contexts: &contexts,
            total_query_span,
            subjects,
            matrix: args.scoring.matrix,
            gap_open: args.scoring.gap_open,
            gap_extend: args.scoring.gap_extend,
            x_drop_final: x_drop_gapped_final,
            karlin: &gapped_params,
            gumbel: gapped_gumbel.as_ref(),
            search_spaces: &search_spaces,
            expect_value: args.evalue,
            hitlist_size: args.traceback_hitlist_size(),
            sorts_results,
            link: args
                .sum_stats
                .then(crate::core::link_hsps::EvenGapLinkParams::traceback),
        };
        Some(ordinary.run(
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            blastp_parallel_pool,
        )?)
    } else {
        None
    };
    // The parameters of each query's redo (`redo_run::RedoParamsBuilder`), which the redo's paths
    // below read.
    let redo_params = redo_run::RedoParamsBuilder {
        args: &args,
        contexts: &contexts,
        gapped_params: &gapped_params,
        hit_cutoff_score_min,
        ideal_ungapped_lambda,
        x_drop_gapped_final,
        initial_matrix_info: std::sync::OnceLock::new(),
    };
    let build_redo_align_params = |q_idx: usize| redo_params.build(q_idx);

    // Each query's batch (`redo_run::kappa_query_batches`): the heaps that its redo's early
    // termination test reads.
    let kappa_query_batches = redo_run::kappa_query_batches(query_ids.len(), &query_batches);
    let evalue_adjust_contexts = redo_run::evalue_adjust_contexts(
        &args,
        &contexts,
        &search_spaces,
        &kappa_query_batches,
        &query_ids,
    )?;
    // What the redo's plan and its paths read (`redo_run::RedoInputs`).
    let redo = redo_run::RedoInputs {
        args: &args,
        blastp_parallel_pool,
        needs,
        subjects,
        contexts: &contexts,
        search_spaces: &search_spaces,
        gapped_params: &gapped_params,
        gapped_gumbel: &gapped_gumbel,
        kappa_query_batches: &kappa_query_batches,
        evalue_adjust_contexts: &evalue_adjust_contexts,
        build_redo_align_params: &build_redo_align_params,
        fetch_excluded: listed_subjects
            .as_ref()
            .map_or(&[][..], |listed| listed.fetch_excluded.as_slice()),
        fetch_notes: &fetch_notes,
        x_win_expand,
        ordinary_traceback,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        blastp_num_threads,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        timing_enabled,
    };
    // The plan chooses the redo's path (`redo_run::RedoPath`) and expands the compact lists that
    // the path reads expanded; the path redoes the matches.
    match redo_run::plan_redo_and_expand_lists(
        &redo,
        &mut preliminary_hit_lists,
        &mut x_compact_hit_lists,
    ) {
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        redo_run::RedoPath::Windows => redo_run::x_windows_redo(
            &redo,
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            &mut redone_matches,
            &mut x_prelim_freed,
        )?,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        redo_run::RedoPath::MatchParallel(q_idx) => redo_run::match_parallel_redo(
            &redo,
            q_idx,
            &preliminary_hit_lists,
            &mut redone_matches,
        )?,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        redo_run::RedoPath::QueryParallel(kappa_parallel_query_indices) => {
            redo_run::query_parallel_redo(
                &redo,
                &kappa_parallel_query_indices,
                &preliminary_hit_lists,
                &mut redone_matches,
            )?
        }
        redo_run::RedoPath::OneThread => redo_run::one_thread_redo(
            &redo,
            &mut preliminary_hit_lists,
            &mut x_compact_hit_lists,
            &mut redone_matches,
            &mut x_prelim_freed,
        )?,
    }

    if timing_enabled {
        eprintln!("[TIMING] x_bppfree: freed_lists={x_prelim_freed}");
        eprintln!(
            "[TIMING] x_bppcompact: expanded_lists={}",
            x_expanded_lists() - x_expanded_at_start
        );
        eprintln!(
            "[TIMING] x_bpihpool: pooled={}",
            x_ih_pool::x_pooled() - x_pooled_at_start
        );
    }
    fetch_notes.after(&redone_matches, &ordinary_hit_lists, &kappa_query_batches)?;
    let mut hit_lists =
        fill_results_from_compo_heaps(args.traceback_hitlist_size(), &mut redone_matches);
    if let Some(ordinary_hit_lists) = ordinary_hit_lists {
        hit_lists = ordinary_hit_lists;
    }
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.kappa_redo_ns,
            &timing.kappa_redo_calls,
            kappa_redo_start,
        );
    }

    let final_traceback_start = blastp_timing_start(timing_enabled);
    traceback::finish_hit_lists(
        &args,
        &mut hit_lists,
        prelim_hitlist_size,
        rps_multi_volume,
        sorts_results,
    );
    // NCBI's Seq-aligns of the hit lists (`LocalBlastResults2SeqAlign` below): none for a subject
    // that the seqid list keeps without a definition line (line DBID item lineless; `-ungapped`:
    // also under a negative list). psiblast's lists are cut to the preliminary size here; the
    // formatter's cut to `-max_target_seqs` below takes the Seq-aligns: no place for it there.
    if let (Some(listed), Some(database)) = (listed_subjects.as_ref(), database) {
        listed.drop_unaligned(database, &mut hit_lists, !args.ungapped);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/traceback_stage.cpp:267-275
    // A database: the definition lines of the final hit lists are read first, where NCBI
    // reads the Seq-ids of its Seq-aligns (`database::ListedSubjects::read_titles`), before
    // the HSPs of score 0 are dropped.
    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the titles of the reported subjects only, and the
    // letters from the stored residues (`database/x_listed_lean.rs`).
    let x_listed_lean = database::x_bp_listed_lean();
    if listed_subjects.is_some() && database.is_some() && x_listed_lean {
        crate::utils::xreach::x_reach!("LOSAT_X_BPLISTEDLEAN", true);
    }
    // What the steps from here to the labels of the reports read (`report_prep::ReportSetup`).
    let setup = report_prep::ReportSetup {
        args: &args,
        reports,
        ranges,
        query_records,
        query_batches: &query_batches,
        database,
        requested_formats: &requested_formats,
        list_cut: outputs
            .formats
            .first()
            .and_then(|format| crate::blastinput::app::parse_formatting_string(format.outfmt).ok())
            .map_or(crate::algorithm::blastp::args::ListCut::Uncut, |choice| {
                args.app.list_cut(choice.number)
            }),
        contexts: &contexts,
        gapped_params: &gapped_params,
        subject_source,
        subject_label,
        total_db_len,
        unlisted: &rps_unlisted,
    };
    report_prep::read_shown_residues(
        &setup,
        outputs,
        &mut listed_subjects,
        &hit_lists,
        x_listed_lean,
    )?;
    let subjects: &[EncodedProtein] = match listed_subjects.as_ref() {
        Some(listed) => &listed.encoded,
        None => subject_set,
    };
    let reads = database::CallTitles::new(database, &query_batches);
    let listed_titles = reads.titles(listed_subjects.as_ref(), &hit_lists, x_listed_lean)?;
    report_prep::refuse_prefetch_environment(&setup, outputs, reads.utf8())?;
    if let (Some(titles), Some(scopes)) = (listed_titles.as_ref(), reports.query_scopes) {
        titles.refuse_ids_of_queries(&hit_lists, scopes)?;
    }
    let (subject_records, subject_titles, shown_subjects) =
        report_prep::shown_of(&listed_titles, subjects, subject_records, &subject_titles);
    let report_prep::ReportValues {
        batches_prefetch,
        prefetch_taxids,
        prefetch_only,
        psiblast_prefetched,
        batch_ranges,
        batch_prefetch,
        pssm_reads,
        mut psiblast_pssms,
        psi_found_before,
        pssm_made,
        psiblast_hits_found,
        engine_counts,
    } = report_prep::report_values(
        &setup,
        &listed_subjects,
        &listed_titles,
        shown_subjects,
        &mut hit_lists,
    )?;

    let report_prep::ReportLabels {
        subject_input_label,
        shown_num_sequences,
        shown_total_letters,
        shown_database,
        context,
    } = report_prep::report_labels(&setup, subject_records);

    // The requested output formats; only the reporting depends on the format, the search above
    // does not.
    let format_plans: &[BlastpReport] = &requested_formats;
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/tabular.cpp:973-986
    // NCBI reference (598d8ae6): c++/src/objtools/align_format/showalign.cpp:1970-1973
    // The final hits in record coordinates and lengths, which the reports print
    // (`shift_hsp_to_records`); the rows are made from them (`rows_of`).
    for hit_list in hit_lists.iter_mut().flatten() {
        for hsp_list in &mut hit_list.hsplist_array {
            for hsp in &mut hsp_list.hsps {
                shift_hsp_to_records(hsp, ranges);
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hspstream.c:271-327
    // The rendered rows of the HSPs of the queries `window` (`PairwiseHit`, in the order of the
    // lists): each row is read from the searched letters, so its hit moves back from the records
    // by the intervals' starts, and the row then prints record coordinates and lengths
    // (`shift_to_records`). Every report that shows the alignment takes its rows from here. The
    // rows' errors are LOSAT's internal invariants (`build_pairwise_hits`: an HSP whose alignment
    // cannot be rendered), so making them while the reports are written changes no run's output,
    // as BLASTN's (`blastn::blast_engine::run`).
    let x_rows_pool = crate::report::x_pairwise_par::pool(blastp_parallel_pool);
    // psiblast's later rounds: the reports read the record's letters
    // (`psi_round::report_contexts`).
    let report_contexts = psi_round::report_contexts(&args, &contexts, query_records);
    let rows_of = |window: std::ops::Range<usize>| -> Result<Vec<PairwiseHit>> {
        let end = window.end.min(hit_lists.len());
        let first = window.start.min(end);
        let mut hits =
            collect_hits_and_methods_from_hit_lists(&hit_lists[first..end], first, &engine_counts);
        for (hit, _, _, _) in hits.iter_mut() {
            let q_shift = ranges.query_input.placements.offset(hit.q_idx as usize);
            let s_shift = ranges.subject_placements.offset(hit.s_idx as usize);
            hit.q_start -= q_shift;
            hit.q_end -= q_shift;
            hit.s_start -= s_shift;
            hit.s_end -= s_shift;
        }
        let mut rows = build_pairwise_hits(
            hits,
            &report_contexts,
            subjects,
            args.scoring.matrix,
            true,
            ShownQueries::of(query_records, &args),
            shown_subjects,
            x_rows_pool,
        )?;
        for row in rows.iter_mut() {
            let s_idx = row.hit.s_idx as usize;
            shift_to_records(&mut row.hit, ranges);
            row.subject_length = row
                .subject_length
                .map(|searched| ranges.subject_placements.length(s_idx, searched));
        }
        Ok(rows)
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411-1585
    // NCBI makes and prints the alignments of one query's result set at a time: the reports
    // make the rows of a window of queries and drop them before the next window's
    // (`query_rows`: outfmt 0 and the tabular fields that show the alignment,
    // `TabularNeeds::rows`, as BLASTN's and TBLASTX's outfmt 0; BLAST XML one query's). The rows of every query are
    // made once, before any format is written, only where a consumer reads them all: the
    // caller's hit records, and a format observer of outfmt 0, whose HSP numbers count over the
    // whole report.
    let rows_up_front = outputs.hits.is_some()
        || (outputs.observer.is_some()
            && format_plans
                .iter()
                .any(|report| matches!(report, BlastpReport::Pairwise(_))));
    let pairwise_hits = if rows_up_front {
        Some(rows_of(0..hit_lists.len())?)
    } else {
        None
    };

    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.final_traceback_ns,
            &timing.final_traceback_calls,
            final_traceback_start,
        );
    }
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Search, t_search);

    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1411
    // The caller receives the same final result that every formatter prints.
    let pairwise_hits = pairwise_hits.map(std::sync::Arc::new);
    if let (Some(hits_sink), Some(pairwise_hits)) = (outputs.hits.as_mut(), pairwise_hits.as_ref())
    {
        hits_sink(std::sync::Arc::clone(pairwise_hits));
    }
    // The library's records name a database's sequences through this event.
    if let (Some(titles), Some(database), true, Some(observer)) = (
        listed_titles.as_ref(),
        database,
        outputs.hits.is_some(),
        outputs.observer.as_deref_mut(),
    ) {
        observer.database_subjects(titles.library_subjects(database));
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
    let output_format_start = blastp_timing_start(timing_enabled);
    // NCBI reference (598d8ae6): c++/src/app/blast/blast_formatter.cpp:429-467
    // NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:1450-1452
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/local_blast.cpp:197-201,221
    // The warnings of each query are written before the query's report, once, with the
    // first format (`QueryWarnings`): the reader's messages of the query's batch (before
    // its first query, `read_query_batches`), then the query's own: the O characters read
    // as X, and the message of its set-up (a query without letters) or of its
    // Karlin-Altschul parameters (an invalid query with letters).
    let filtered_empty = database.is_some_and(|database| database.filtered_empty());
    // The warning lines of each query, the checks of the reported subjects and each batch's
    // prefetch of them, and the refusals of what LOSAT does not report as NCBI does, before any
    // format is written (`report_prep`).
    let prep = report_prep::ReportPrep {
        args: &args,
        outputs,
        reports,
        ranges,
        query_records,
        query_valid: &query_valid,
        query_batch_skipped: &query_batch_skipped,
        query_batches: &query_batches,
        database,
        listed_titles: &listed_titles,
        subject_records,
        subject_titles: &subject_titles,
        format_plans,
        hit_lists: &hit_lists,
        pairwise_hits: &pairwise_hits,
        engine_counts: &engine_counts,
        psi_found_before: &psi_found_before,
        psiblast_hits_found: &psiblast_hits_found,
        filtered_empty,
        rps_multi_volume,
        batches_prefetch,
        prefetch_taxids,
        prefetch_only,
        psiblast_prefetched: &psiblast_prefetched,
        batch_ranges: &batch_ranges,
        batch_prefetch: &batch_prefetch,
        pssm_made,
        pssm_reads: &pssm_reads,
    };
    let (mut query_warning_lines, kept) = report_prep::query_warning_lines(&prep);
    let prefetched = report_prep::prefetch_reported_subjects(
        &prep,
        &mut query_warning_lines,
        &mut psiblast_pssms,
    )
    .map_err(|error| database::HeaderReadError::after_reads(error, reads.utf8()))?;
    // An ungapped HSP of length 0 at the first letter of a `-subject`'s searched part ends the call
    // with NCBI's exception after the formats' prologs (`empty_alignments`).
    let subject_start_exception = empty_alignments::subject_start_exception(
        &hit_lists,
        ranges.subject_placements,
        database.is_none(),
    );
    let unreportable = report_prep::refuse_unreportable(&prep, &rows_of);
    // It was the first of `refuse_unreportable`'s refusals and its error still comes first; it is
    // made after them because they read `prep`, which holds the outputs (line WQA's sweep, D7,
    // 2026-11-19).
    refuse_shared_record_masks(
        &args,
        format_plans,
        outputs,
        query_records,
        query_frames,
        &ranges.query_input.placements,
    )?;
    unreportable?;
    // The subject IDs of the tabular formats, checked for the subjects with rows before any
    // format is written (`tabular_subject_ids`). EXPERIMENT (LOSAT_X_BPQBLOCK): the block
    // search's IDs made once for all blocks (`x_subject_meta`).
    // `sacc` and `saccver` of the database sequences with parsed Seq-ids, by rank (empty
    // without them: the others' are their `sseqid`).
    let mut subject_accessions: Vec<Option<[Arc<[u8]>; 2]>> = Vec::new();
    // The records of a subject set whose ids `-parse_deflines` read, by subject index.
    let records_parsed = if args.parse_deflines && listed_titles.is_none() {
        crate::report::db_subject::parsed_records(
            subject_records
                .iter()
                .map(|record| (record.believed(), record.title.as_slice())),
            args.show_gis,
        )
    } else {
        Vec::new()
    };
    // LOSAT's refusals of these records in outfmt 0 to 4, each such format's (`FormatRefusal`).
    if !records_parsed.is_empty() {
        for (format, outfmt) in outputs.formats.iter_mut().zip(format_plans.iter()) {
            if matches!(outfmt, BlastpReport::Pairwise(_)) {
                format
                    .refusal
                    .check(crate::report::db_subject::check_pairwise_records(
                        &records_parsed,
                        true,
                        args.app.name(),
                    ))?;
            }
        }
    }
    let subject_ids: std::borrow::Cow<'_, [crate::algorithm::common::tabular::SubjectId]> =
        if format_plans.iter().any(|outfmt| outfmt.is_tabular()) {
            match (listed_titles.as_ref(), x_meta) {
                (None, _) if !records_parsed.is_empty() => {
                    let (ids, accessions) =
                        record_tabular_ids(subject_records, &records_parsed, args.app.name())?;
                    subject_accessions = accessions;
                    ids.into()
                }
                (Some(titles), _) => {
                    let (ids, accessions) = titles.tabular_subject_ids(
                        &hit_lists,
                        format_plans
                            .iter()
                            .any(|outfmt| outfmt.is_tabular() && outfmt.reads_subjects()),
                    )?;
                    subject_accessions = accessions;
                    ids.into()
                }
                // EXPERIMENT (LOSAT_X_BPQBLOCK): the block search's ids made once for all
                // blocks (`x_subject_meta`), `tabular_subject_ids` of the same records.
                (None, Some(meta)) => meta.shown().into(),
                (None, None) => tabular_subject_ids(subject_records).into(),
            }
        } else {
            Vec::new().into()
        };
    // The rows of subjects whose lookups fail where NCBI has no result, before any format is
    // written (`check_tabular_lookups`).
    let parsed_for_checks: &[Option<Arc<crate::report::db_subject::ParsedSubject>>] = listed_titles
        .as_ref()
        .map_or(records_parsed.as_slice(), |titles| &titles.parsed);
    check_tabular_lookups(
        &mut outputs.formats,
        format_plans,
        &hit_lists,
        &subject_ids,
        parsed_for_checks,
        &prefetched.hit_loads,
        listed_titles.as_ref(),
        subject_records,
        args.app.name(),
    )?;
    // The `stitle` text of the custom tabular reports by subject: the stored titles of a
    // database (the first definition line's of a sequence with parsed Seq-ids); none for the
    // records of a subject set, which have no BLAST definition line (`N/A`).
    // EXPERIMENT (LOSAT_X_BPLISTEDLEAN): the titles that outfmt 0 borrows too.
    let subject_stitles: std::borrow::Cow<'_, [Arc<[u8]>]> = match listed_titles.as_ref() {
        Some(titles)
            if format_plans
                .iter()
                .any(|outfmt| matches!(outfmt, BlastpReport::CustomTabular(_))) =>
        {
            match &titles.shown {
                database::ListedShown::Lean(titles) => std::borrow::Cow::Borrowed(&titles[..]),
                database::ListedShown::Records(_) => titles.stitles().into(),
            }
        }
        _ => Vec::new().into(),
    };
    // `-db_soft_mask` and `-db_hard_mask`: the `Masked using:` lines of each name.
    let pairwise_database =
        database.map(|database| database.pairwise_database_masked(args.db_mask_resolved.as_ref()));
    // The reported database sequences with parsed Seq-ids by rank (empty without them).
    let listed_parsed: &[Option<Arc<crate::report::db_subject::ParsedSubject>>] = listed_titles
        .as_ref()
        .map_or(records_parsed.as_slice(), |titles| &titles.parsed);
    // The sources from which NCBI makes the empty stored titles of the others, by rank (empty
    // without them; `report::pairwise::stored_source`).
    let listed_sources: &[crate::report::defline::SourceTaxname] = listed_titles
        .as_ref()
        .map_or(&[][..], |titles| &titles.sources);
    // The warnings are written once, with the first format (`QueryWarnings`).
    // NCBI makes each query's PSSM after the query's report and writes its file then
    // (`psiblast_pssms`): the warnings of its reads and the posts of the file's write
    // (`PssmFiles::take_read_warnings`) come before the next query's warnings, and the last
    // query's after its report, before the epilog (psiblast_app.cpp:304-308, 319-321 and 562,
    // `PrintEpilog`; review 2 of line DBID pssmread, H4), where the first format's writer puts them
    // (`QueryWarnings::after_last_query`; the others have no epilog: after the reports).
    let mut pssm_posts_after = crate::report::diag_log::Posts::new();
    // Those posts at the start of each query's posts, before the reader's messages of its batch
    // (`query_posts` below).
    let mut pssm_posts_before: Vec<crate::report::diag_log::Posts> = Vec::new();
    if let Some(files) = psiblast_pssms.as_mut() {
        let queries = query_warning_lines.len();
        for (query, posts) in files.take_read_warnings() {
            let before = if query + 1 < queries {
                query + 1
            } else if format_plans
                .first()
                .is_some_and(|outfmt| outfmt.writes_after_last_query())
            {
                query_warning_lines.resize(queries + 1, Vec::new());
                queries
            } else {
                pssm_posts_after.append(&posts);
                continue;
            };
            if pssm_posts_before.len() <= before {
                pssm_posts_before.resize(before + 1, crate::report::diag_log::Posts::new());
            }
            pssm_posts_before[before].append(&posts);
        }
    }
    // Each requested format prints the same final result without searching again.
    let mut output_result: Result<()> = Ok(());
    // The statistics of each query: what the XML reports print and what the observer receives
    // (`query_xml_stats`).
    let query_stats = || {
        query_xml_stats(
            &args,
            gapped_params,
            &contexts,
            &query_valid,
            &search_spaces,
            &exact_search_spaces,
            rps_gapped_k.as_deref(),
            query_records.len(),
        )
    };
    // Where the run's XML formats fail for the matrix name (`blastp_xml_failing_query`).
    let scopes = reports.query_scopes.map(|scopes| scopes.run());
    let xml_matrix_failure = blastp_xml_failing_query(&args, outputs, database.is_some(), scopes);
    let observer = &mut outputs.observer;
    if let Some(observer) = observer.as_deref_mut() {
        for (q_idx, stats) in query_stats().into_iter().enumerate() {
            observer.query_stats(q_idx, stats);
        }
    }
    let rps_request_spaces = rps_request_spaces(
        query_batches.iter().map(|batch| batch.searched.start),
        rps_units,
        &args,
        &contexts,
        &gapped_params,
    )?;
    let first_contexts = FirstContexts {
        contexts: &contexts,
        query_valid: &query_valid,
        searchsp: args.searchsp,
        exact_search_spaces: &exact_search_spaces,
        search_spaces: &search_spaces,
        request_spaces: rps_request_spaces.as_deref(),
        rps_gapped_k: rps_gapped_k.as_deref(),
        gapped_params,
        ungapped: args.ungapped,
    };
    let first_context = |first: usize| first_contexts.first_context(first);
    // The queries of the call's batches, as the archive and the search strategy take them.
    let batch_spans = || -> Vec<std::ops::Range<usize>> {
        query_batches
            .iter()
            .map(|batch| batch.searched.clone())
            .collect()
    };
    // (`masks` and `spans` are each format's.)
    let call_batches = super::archive::CallBatches {
        queries: query_records,
        query_placements: &ranges.query_input.placements,
        subject_records: database.is_none().then_some(subject_records),
        masks: &[],
        lcase_masking: args.lcase_masking,
        unsearched: &query_batch_skipped,
        spans: &[],
    };
    // The search strategy's requests of the call's batches (`-export_search_strategy`): NCBI has
    // written every one of them before it formats them, so they are written here (`CallRequests`;
    // a matrix name that the XML formats do not know stops them after the batch whose formatting
    // fails, `search`).
    if let Some(call) = reports.call_requests {
        call.end_at(query_records.len())?;
    }
    // The posts before each query's report: NCBI's posts of the bytes replaced in the search
    // strategy request of each of the call's batches follow the reader's messages of the batch
    // (after psiblast's warnings of the PSSM reads of the query before), in its first query
    // (`query_warnings::query_posts`).
    let mut query_posts = crate::report::query_warnings::query_posts(
        query_warning_lines,
        reports
            .call_requests
            .map_or_else(Vec::new, CallRequests::take_posts),
        prefetched.posts,
        |query| reports.reader_warnings.get(query).map_or(0, Vec::len),
        reads.header_read(args.app.ncbi_name()),
    );
    for (query, mut posts) in pssm_posts_before.into_iter().enumerate() {
        if let Some(after) = query_posts.get_mut(query).filter(|_| !posts.is_empty()) {
            posts.append(after);
            *after = posts;
        }
    }
    // NCBI reference (598d8ae6): blast_format.cpp:1450-1452 (at `QueryWarnings::with_later`).
    let query_posts = Posts::with_later(query_posts, kept.as_deref().unwrap_or_default());
    let mut warning_lines: Option<&[crate::report::diag_log::Posts]> = Some(&query_posts);
    // Every structured format of the call is XML2 or JSON: none writes the report that a
    // `-subject` document holds back (`blast_xml::write_blastp_xml`).
    let xml_only_xml2 = !format_plans.contains(&BlastpReport::Xml);
    // What the reports of the loop over the formats read of the call (`report_write`).
    let report_view = report_write::ReportView {
        args: &args,
        reports,
        ranges,
        query_records,
        query_frames,
        query_ids: &query_ids,
        query_titles: &query_titles,
        query_lengths: &query_lengths,
        query_valid: &query_valid,
        query_batch_skipped: &query_batch_skipped,
        contexts: &contexts,
        report_contexts: &report_contexts,
        gapped_params,
        gapped_gumbel,
        search_spaces: &search_spaces,
        exact_search_spaces: &exact_search_spaces,
        database,
        listed_subjects: &listed_subjects,
        listed_titles: &listed_titles,
        listed_parsed,
        listed_sources,
        pairwise_database: &pairwise_database,
        shown_database: &shown_database,
        subject_input_label: &subject_input_label,
        shown_num_sequences,
        shown_total_letters,
        subjects,
        subject_records,
        subject_titles: &subject_titles,
        shown_subjects,
        subject_ids: &subject_ids,
        subject_accessions: &subject_accessions,
        subject_stitles: &subject_stitles,
        hit_lists: &hit_lists,
        hit_loads: &prefetched.hit_loads,
        engine_counts: &engine_counts,
        pairwise_hits: &pairwise_hits,
        psiblast_hits_found: &psiblast_hits_found,
        psi_found_before: &psi_found_before,
        rps_gapped_k: &rps_gapped_k,
        no_alignment_set: &no_alignment_set,
        context: &context,
        call_batches: &call_batches,
        filtered_empty,
        rps_multi_volume,
        xml_matrix_failure,
        xml_only_xml2,
        x_rows_pool,
        blastp_parallel_pool,
        prolog_only: subject_start_exception.is_some(),
    };
    for (format_index, (format, &outfmt)) in outputs
        .formats
        .iter_mut()
        .zip(format_plans.iter())
        .enumerate()
    {
        // A format that failed alone writes nothing more; the next takes the warnings.
        if format.refusal.refused() {
            continue;
        }
        // NCBI makes the titles of the subjects that the pairwise report shows, those with hits
        // (blast_format.cpp:1540); LOSAT rejects the titles that it does not write as NCBI
        // (`check_shown_subject_titles`), before the report's first byte and before the format
        // takes the warnings that the formats share, so that a format that fails alone
        // (`FormatRefusal::check`) leaves them to the next. The shown subjects are those of the
        // final lists (the rows' subjects).
        let pairwise_shown: Option<std::collections::BTreeSet<usize>> = match outfmt {
            BlastpReport::Pairwise(_) => {
                let shown: std::collections::BTreeSet<usize> = hit_lists
                    .iter()
                    .flatten()
                    .flat_map(|hit_list| hit_list.hsplist_array.iter().take(hit_list.hsplist_count))
                    .flat_map(|hsp_list| &hsp_list.hsps)
                    .map(|hsp| hsp.s_idx as usize)
                    .collect();
                let checked = check_pairwise_shown_titles(
                    listed_titles.as_ref(),
                    database,
                    &shown,
                    subject_records,
                    &prefetched.hit_loads,
                );
                if !format.refusal.check(checked)? {
                    continue;
                }
                Some(shown)
            }
            _ => None,
        };
        // The warnings that a format other than the archive writes are not the archive's
        // messages (a call with several formats, `run_local`).
        if let Some(handle) = args.archive.as_ref() {
            handle.set_muted(outfmt != BlastpReport::Archive);
        }
        let mut probe = observer
            .as_deref_mut()
            .map(|observer| FormatProbe::new(observer, format_index));
        // The caller's receiver of the side files of outfmt 13 and 14, if it gives one.
        let side_sink = format.sink.side_files().cloned();
        // The taxonomy database's messages of XML2 and JSON, and of the organism report's
        // lookups, go to the diagnostics at every format, and the archive saves the warnings of
        // its batches' queries; so do the Seq-annots' posts of the bytes that NCBI's stream
        // replaces (`seqalign::write_query_to`).
        // The first format's writer also writes the texts of the subjects that the reports
        // read one by one (`ReportedPrefetch::subject_reads`).
        let first_lines = warning_lines.take();
        let own = format.later_warnings(kept.is_some(), &mut outputs.later_warnings, format_index);
        let mut warnings = (first_lines.is_some()
            || own.is_some()
            || crate::report::query_warnings::report_posts(
                format.report_format(),
                args.html.is_some(),
            )
            || matches!(
                outfmt,
                BlastpReport::CustomTabular(spec)
                    if spec.kind == crate::report::tabular::TabularKind::TaxReport
            ))
        .then(|| {
            let warnings = crate::report::query_warnings::QueryWarnings::of_posts(
                first_lines.unwrap_or(&[]),
                args.archive.as_ref(),
                &mut *outputs.diagnostics,
            )
            .with_later(&query_posts, own)
            .for_report(format_index)
            .with_html_posts(
                outputs
                    .html_posts
                    .as_deref_mut()
                    .map(|posts| posts as &mut dyn std::io::Write),
            );
            match first_lines {
                Some(_) => warnings.with_subject_reads(prefetched.subject_reads.queries()),
                None => warnings,
            }
        });
        let format_result: Result<()> = match outfmt {
            _ if report_view.prolog_only => empty_alignments::write_exception_prolog(
                &report_view,
                outfmt,
                format,
                &mut warnings,
                &mut probe,
                &pairwise_shown,
                &rows_of,
            ),
            BlastpReport::Seqalign(serial) => {
                report_write::write_seqalign(&report_view, serial, format, &mut warnings)
            }
            BlastpReport::Archive => report_write::write_archive(
                &report_view,
                format,
                &mut warnings,
                &first_context,
                &batch_spans,
            ),
            BlastpReport::Xml | BlastpReport::Xml2(_) => report_write::write_xml(
                &report_view,
                outfmt,
                format,
                &mut warnings,
                &side_sink,
                &query_stats,
                &rows_of,
            ),
            BlastpReport::Pairwise(view) => report_write::write_pairwise(
                &report_view,
                view,
                format,
                &mut warnings,
                &mut probe,
                &pairwise_shown,
                &rows_of,
            ),
            BlastpReport::CustomTabular(spec) => report_write::write_custom_tabular(
                &report_view,
                spec,
                format,
                &mut warnings,
                &mut probe,
                &rows_of,
            ),
            BlastpReport::Text(
                outfmt @ (OutputFormat::Tabular | OutputFormat::TabularWithComments),
            ) => report_write::write_text(&report_view, outfmt, format, &mut warnings, &mut probe),
            // Outfmt 0 is `BlastpReport::Pairwise` (`requested_formats`).
            BlastpReport::Text(OutputFormat::Pairwise) => Err(anyhow::anyhow!(
                "LOSAT internal error: outfmt 0 as a BLASTP tabular report"
            )),
        };
        // An error of the format's report fails the call, or this format alone (`FormatRefusal`):
        // the formats after it are written then.
        if let Err(error) = format.refusal.check(format_result) {
            output_result = Err(error);
            break;
        }
    }
    if let Some(handle) = args.archive.as_ref() {
        handle.set_muted(false);
    }
    if let Some(timing) = timing.as_ref() {
        BlastpTiming::record_call(
            &timing.output_format_ns,
            &timing.output_format_calls,
            output_format_start,
        );
    }
    // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs): the
    // formatting and writing of every report, after the search. No NCBI counterpart.
    crate::utils::xtiming::add(crate::utils::xtiming::Phase::Output, output_format_start);
    if let Some(exception) = subject_start_exception {
        output_result =
            output_result.and(Err(crate::cli::NativeError::exception(&exception).into()));
    }
    output_result?;
    pssm_posts_after.write(&mut *outputs.diagnostics, args.archive.as_ref())?;
    // NCBI reference (598d8ae6): c++/src/app/blast/psiblast_app.cpp:319-321
    // The files were staged before the reports (`psiblast_pssms`); an error before this point
    // has dropped them, and no PSSM file has changed (`PssmFiles`).
    if let Some(files) = psiblast_pssms {
        files.commit()?;
    }
    // EXPERIMENT (LOSAT_X_BPLKPFLAT): the flat table's vectors kept for the next block's build
    // (`x_query_blocks::x_bp_lookup_spare`) instead of being freed. EXPERIMENT
    // (LOSAT_X_BPLKPAHEAD): exact blastp's, for the build of a later block's table ahead of it
    // (`x_query_blocks::x_bp_aa_spares`); EXPERIMENT (LOSAT_X_BPLKREUSE): or for the next block's own
    // build (`x_query_blocks::x_bp_aa_block_spare`). EXPERIMENT (LOSAT_X_BPLKDROP): `None` when the
    // table was dropped before the redo.
    match lookup {
        Some(BlastpRuntimeLookup::CompressedFlat(table)) => {
            if let Some(spare) = super::x_query_blocks::x_bp_lookup_spare(blastp_parallel_pool) {
                spare.give(table);
            }
        }
        Some(BlastpRuntimeLookup::Aa(table)) if x_bp_aa_vectors_reused(&table) => {
            if let Some(spares) = super::x_query_blocks::x_bp_aa_spares(blastp_parallel_pool)
                .or_else(|| super::x_query_blocks::x_bp_aa_block_spare(blastp_parallel_pool))
            {
                spares.give(table);
            }
        }
        // EXPERIMENT (LOSAT_X_BPLK16): the vectors of a table of 16-bit hits, for the next block's
        // build of that width.
        Some(BlastpRuntimeLookup::Aa16(table)) => {
            if let Some(spares) = super::x_query_blocks::x_bp_aa_spares(blastp_parallel_pool)
                .or_else(|| super::x_query_blocks::x_bp_aa_block_spare(blastp_parallel_pool))
            {
                spares.give16(table);
            }
        }
        Some(BlastpRuntimeLookup::Aa(_) | BlastpRuntimeLookup::Compressed(_)) | None => {}
    }

    if let Some(timing) = timing.as_ref() {
        print_blastp_timing(timing, t_total);
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1633-1705
    if kappa_range_counters_enabled {
        print_blastp_kappa_range_counters();
    }

    Ok(())
}

/// NCBI's error of a query batch whose score block cannot load an unknown `-matrix` name
/// (`blastx::query_setup::unknown_matrix_message`): the error of each query, then the messages of
/// its reader, as NCBI's `TSearchMessages::ToString` writes them: the O warning of a protein
/// query and `Sequence contains no data` of a record without letters (`query_messages`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/setup_factory.cpp:158-182 (cited at
/// `blastx::query_setup::unknown_matrix_message`)
fn unknown_matrix_message(records: &[FastaRecord]) -> String {
    crate::algorithm::common::unknown_matrix::unknown_matrix_message_of(records.iter().map(
        |record| crate::report::query_warnings::query_messages(record.seq(), true, false, false),
    ))
}

// NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_gapalign.h:53-54
/// The refusal of an ungapped search of a subject longer than `MAX_DBSEQ_LEN` (`what` names it):
/// NCBI searches it in overlapping parts and merges their HSP lists (`Blast_HSPListsMerge`,
/// blast_engine.c:578-584), which `-ungapped` does not port (BLASTX refuses the same subject,
/// `blastx::database::too_long`).
fn too_long_subject_refusal(what: &str, length: usize) -> anyhow::Error {
    anyhow::anyhow!(
        "{what} has {length} residues; a subject longer than {} residues, which NCBI BLAST+ searches in overlapping parts, is not supported by LOSAT's BLASTP",
        crate::core::subject_split::MAX_DBSEQ_LEN
    )
}

#[cfg(test)]
#[path = "blast_engine_tests.rs"]
mod tests;
