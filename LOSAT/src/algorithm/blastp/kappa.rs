//! `blast_kappa.c`-style post-processing helpers for the default `blastp`
//! parity path.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c
//! NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c

use crate::utils::glibc_fma_libm::GlibcLibm;
use anyhow::Result;
use std::collections::HashMap;
use std::ffi::c_void;
use std::ptr::NonNull;
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering as AtomicOrdering};
use std::sync::OnceLock;

use crate::algorithm::tblastx::lookup::QueryContext;
use crate::common::{GapEditOp, Hit};
use crate::config::ScoringMatrix;
use crate::core::blast_seg::{SegMasker, SegParams};
use crate::core::composition_adjustment::adjust_scores::{
    blast_overall_p_value, read_aa_composition, x_newton_input, AdjustedProteinMatrix,
    BlastCompositionWorkspace,
};
use crate::core::composition_adjustment::redo_alignment::{
    blast_compo_alignment_new, blast_redo_one_match_with_workspace_queries_and_matrix,
    build_query_word_hashes, preliminary_test_near_identical, redo_calc_lambda,
    test_near_identical, windows_from_protein_aligns, x_whole_subject_composition,
    BlastCompoAlignment, BlastCompoAlignmentContext, BlastCompoMatchingSequence,
    BlastCompoQueryInfo, BlastCompoSequenceData, BlastCompoSequenceRange, BlastRedoAlignCallbacks,
    BlastRedoRangeResult, EMatrixAdjustRule, XSubjectComposition,
};
use crate::stats::search_space::SearchSpace;
use crate::stats::{blast_spouge_stoe, e_to_p, p_to_e, BlastGumbelBlk, KarlinParams};
use crate::utils::matrix::{blosum62_score_ncbistdaa_direct, ncbistdaa, protein_score};

use super::gapalign::{
    blast_gapped_alignment_with_traceback_with_scratch, BlastpPreliminaryHsp, GapAlignScratch,
};
use super::hsp::{
    blast_compo_early_termination, reap_hsplist_by_evalue, sort_hsplist_by_score,
    update_best_evalue, BlastCompoHeap, BlastpHsp, BlastpHspList,
};

pub(crate) use crate::core::composition_adjustment::redo_alignment::BlastRedoAlignParams;

// NCBI's Smith-Waterman redo of a match with a square matrix (`-use_sw_tback`; item X4).
pub(crate) mod smith_waterman;

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1614-1646
//
// LOSAT diagnostics count the corresponding query/subject range copies and
// SEG range scans without changing the redo-alignment control flow.
static BLASTP_KAPPA_RANGE_COUNTERS_ENABLED: AtomicBool = AtomicBool::new(false);
static BLASTP_KAPPA_QUERY_RANGE_BYTES_COPIED: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_RANGE_BYTES_COPIED: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_SEG_RANGE_BYTES: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_SEG_CALLS: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_SEG_MASKED_RESIDUES: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_SEG_CACHE_HITS: AtomicU64 = AtomicU64::new(0);
static BLASTP_KAPPA_SUBJECT_SEG_CACHE_MISSES: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Copy, Debug)]
pub(crate) struct BlastpTraceHspTarget {
    pub(crate) q_idx: Option<u32>,
    pub(crate) s_idx: Option<u32>,
    pub(crate) q_start: usize,
    pub(crate) q_end: usize,
    pub(crate) s_start: usize,
    pub(crate) s_end: usize,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3658-3716
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:746-824
pub(crate) fn blastp_trace_hsp_target() -> Option<BlastpTraceHspTarget> {
    if !crate::utils::xdiag::ON {
        return None;
    }
    static TARGET: OnceLock<Option<BlastpTraceHspTarget>> = OnceLock::new();
    *TARGET.get_or_init(|| {
        let value = crate::utils::xdiag::var("LOSAT_DEBUG_BLASTP_TRACE_HSP").ok()?;
        let fields: Vec<&str> = value.split(',').map(str::trim).collect();
        match fields.as_slice() {
            [q_idx, s_idx, q_start, q_end, s_start, s_end] => Some(BlastpTraceHspTarget {
                q_idx: q_idx.parse::<u32>().ok(),
                s_idx: s_idx.parse::<u32>().ok(),
                q_start: q_start.parse().ok()?,
                q_end: q_end.parse().ok()?,
                s_start: s_start.parse().ok()?,
                s_end: s_end.parse().ok()?,
            }),
            [q_start, q_end, s_start, s_end] => Some(BlastpTraceHspTarget {
                q_idx: None,
                s_idx: None,
                q_start: q_start.parse().ok()?,
                q_end: q_end.parse().ok()?,
                s_start: s_start.parse().ok()?,
                s_end: s_end.parse().ok()?,
            }),
            _ => None,
        }
    })
}

#[inline]
pub(crate) fn blastp_trace_matches_pair(
    target: BlastpTraceHspTarget,
    q_idx: u32,
    s_idx: u32,
) -> bool {
    target.q_idx.is_none_or(|expected| expected == q_idx)
        && target.s_idx.is_none_or(|expected| expected == s_idx)
}

#[inline]
fn blastp_trace_matches_hsp(target: BlastpTraceHspTarget, hsp: &BlastpHsp) -> bool {
    blastp_trace_matches_pair(target, hsp.q_idx, hsp.s_idx)
        && hsp.q_start == target.q_start
        && hsp.q_end == target.q_end
        && hsp.s_start.min(hsp.s_end) == target.s_start
        && hsp.s_start.max(hsp.s_end) == target.s_end
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_gapalign.c:364-366
fn blastp_trace_gap_summary(gap_info: Option<&[GapEditOp]>) -> String {
    let Some(gap_info) = gap_info else {
        return "ungapped".to_string();
    };
    let mut sub_ops = 0usize;
    let mut del_ops = 0usize;
    let mut ins_ops = 0usize;
    let mut total_len = 0usize;
    for op in gap_info {
        total_len += op.num() as usize;
        match op {
            GapEditOp::Sub(_) => sub_ops += 1,
            GapEditOp::Del(_) => del_ops += 1,
            GapEditOp::Ins(_) => ins_ops += 1,
        }
    }
    format!(
        "ops={} sub_ops={} del_ops={} ins_ops={} total_len={}",
        gap_info.len(),
        sub_ops,
        del_ops,
        ins_ops,
        total_len
    )
}

#[inline]
pub(crate) fn blastp_kappa_range_counters_env_enabled() -> bool {
    std::env::var_os("LOSAT_TIMING").is_some()
        || std::env::var("LOSAT_DIAGNOSTICS")
            .map(|value| value == "1" || value.eq_ignore_ascii_case("true"))
            .unwrap_or(false)
}

pub(crate) fn reset_blastp_kappa_range_counters(enabled: bool) {
    BLASTP_KAPPA_QUERY_RANGE_BYTES_COPIED.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_RANGE_BYTES_COPIED.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_SEG_RANGE_BYTES.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_SEG_CALLS.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_SEG_MASKED_RESIDUES.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_SEG_CACHE_HITS.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_SUBJECT_SEG_CACHE_MISSES.store(0, AtomicOrdering::Relaxed);
    BLASTP_KAPPA_RANGE_COUNTERS_ENABLED.store(enabled, AtomicOrdering::Relaxed);
}

#[inline]
fn record_blastp_kappa_range_counter(counter: &AtomicU64, value: usize) {
    if BLASTP_KAPPA_RANGE_COUNTERS_ENABLED.load(AtomicOrdering::Relaxed) {
        counter.fetch_add(value.min(u64::MAX as usize) as u64, AtomicOrdering::Relaxed);
    }
}

#[inline]
fn record_blastp_kappa_range_event(counter: &AtomicU64) {
    if BLASTP_KAPPA_RANGE_COUNTERS_ENABLED.load(AtomicOrdering::Relaxed) {
        counter.fetch_add(1, AtomicOrdering::Relaxed);
    }
}

pub(crate) fn print_blastp_kappa_range_counters() {
    eprintln!(
        "[TIMING] blastp_kappa_query_range_bytes_copied: {}",
        BLASTP_KAPPA_QUERY_RANGE_BYTES_COPIED.load(AtomicOrdering::Relaxed)
    );
    eprintln!(
        "[TIMING] blastp_kappa_subject_range_bytes_copied: {}",
        BLASTP_KAPPA_SUBJECT_RANGE_BYTES_COPIED.load(AtomicOrdering::Relaxed)
    );
    eprintln!(
        "[TIMING] blastp_kappa_subject_seg_range: {} bytes (calls={}, masked_residues={})",
        BLASTP_KAPPA_SUBJECT_SEG_RANGE_BYTES.load(AtomicOrdering::Relaxed),
        BLASTP_KAPPA_SUBJECT_SEG_CALLS.load(AtomicOrdering::Relaxed),
        BLASTP_KAPPA_SUBJECT_SEG_MASKED_RESIDUES.load(AtomicOrdering::Relaxed)
    );
    eprintln!(
        "[TIMING] blastp_kappa_subject_seg_cache: hits={} misses={}",
        BLASTP_KAPPA_SUBJECT_SEG_CACHE_HITS.load(AtomicOrdering::Relaxed),
        BLASTP_KAPPA_SUBJECT_SEG_CACHE_MISSES.load(AtomicOrdering::Relaxed)
    );
}

#[inline(always)]
fn blastp_positive_score<const BLOSUM62: bool>(
    matrix: ScoringMatrix,
    query_residue: u8,
    subject_residue: u8,
) -> bool {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-818
    if BLOSUM62 {
        blosum62_score_ncbistdaa_direct(query_residue, subject_residue) > 0
    } else {
        protein_score(matrix, query_residue, subject_residue) > 0
    }
}

#[derive(Debug)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3675-3716
pub(crate) struct BlastpPostprocessResult {
    pub hsp_list: Option<BlastpHspList>,
    pub best_score: i32,
    pub best_evalue: f64,
}

#[derive(Debug, Clone)]
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2336
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:432-435
pub(crate) struct BlastpKappaQueryWorkspace {
    pub redo_query_info: BlastCompoQueryInfo,
    pub query_length: i32,
    pub length_adjustment: i32,
    /// What the adjustment of the query's E-values with a composition P-value reads
    /// ([`evalue_adjust_contexts`]): not the query's own values, in general.
    pub evalue_adjust_context: EvalueAdjustContext,
}

/// The three values of `queryInfo->contexts[context_index]` that NCBI's redo reads when it
/// adjusts a query's E-values with a composition P-value (the unified P-values of
/// `-comp_based_stats 1u`, `2u` and `3u`; `adjust_evalues_for_composition`).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct EvalueAdjustContext {
    pub query_length: i32,
    pub length_adjustment: i32,
    pub eff_searchsp: f64,
}

/// What LOSAT's refusal of unified P-values names ([`evalue_adjust_contexts`]): the program
/// ("... are not supported by LOSAT's BLASTP") and each query's id as the reports show it.
#[derive(Clone, Copy)]
pub(crate) struct UnifiedPRefusal<'a> {
    pub program: &'a str,
    pub query_ids: &'a [std::sync::Arc<[u8]>],
}

/// LOSAT's refusal of unified P-values for query `query`, whose E-values NCBI adjusts with the
/// search space 0 of the query after it ([`evalue_adjust_contexts`]).
fn unified_p_without_search_space(refusal: UnifiedPRefusal<'_>, query: usize) -> anyhow::Error {
    let id = |q: usize| String::from_utf8_lossy(&refusal.query_ids[q]).into_owned();
    anyhow::anyhow!(
        "unified P-values (a 'u' after the -comp_based_stats mode) are not supported by LOSAT's {} for query {}, which is followed in its query batch by a query without a search space ({}: no residue that can be scored, and -searchsp is not given or is 0): NCBI BLAST+ 2.17.0 adjusts a query's e-values with the search space of the next query",
        refusal.program,
        id(query),
        id(query + 1)
    )
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3542-3543,3590-3592,3675-3683
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_setup.c:773-784,846-847
/// Each query's [`EvalueAdjustContext`], in query order: what NCBI's
/// `&queryInfo->contexts[context_index]` holds when the redo adjusts the query's E-values with a
/// composition P-value.
///
/// NCBI passes the context index after the loop over the query's frames has advanced it, so for
/// BLASTP (one frame) the context is the NEXT query's of the query batch (`query_batches`, each
/// query's batch as a range: `kappa_query_batches`): its length, its length adjustment and its
/// effective search space. For a query that is not valid the length adjustment is 0 and the
/// search space is the `-searchsp` value (0 when not given). Reproduced (PD-LOSAT-NCBI-DEFECTS
/// rule 2).
///
/// For the last query of a batch NCBI reads past the batch's contexts (heap bytes: its reports
/// differ from run to run). LOSAT reads the query's own context there, which is what NCBI
/// reads for that query when a copy of it follows in the same batch and the batch is searched
/// in the same chunks (the exception of PD-LOSAT-NCBI-DEFECTS for unified P-values; the redo
/// also ends no match early then, `ResolvedBlastpArgs::redone_match_heap`).
///
/// A next query whose effective search space is 0 makes NCBI's scale infinite: NCBI reports
/// every E-value of the query before it as 0, an alignment whose E-value was 0 without an
/// E-value (the product of 0 and infinity is not a number), and drops a match of such alignments
/// only. With unified P-values that batch is refused (`refuse` is `Some`; the error is a
/// `HeaderReadError` of the batch's first query, so the reports of the batches before are
/// written first, and its text names the two queries).
pub(crate) fn evalue_adjust_contexts(
    contexts: &[QueryContext],
    search_spaces: &[SearchSpace],
    query_batches: &[std::ops::Range<usize>],
    options_searchsp: i64,
    refuse: Option<UnifiedPRefusal<'_>>,
) -> Result<Vec<EvalueAdjustContext>> {
    let context_of = |q_idx: usize| {
        let ctx = &contexts[q_idx];
        if ctx.is_valid && ctx.aa_len > 0 {
            EvalueAdjustContext {
                query_length: ctx.aa_len as i32,
                length_adjustment: search_spaces[q_idx].length_adjustment as i32,
                eff_searchsp: search_spaces[q_idx].effective_space,
            }
        } else {
            EvalueAdjustContext {
                query_length: ctx.aa_len as i32,
                length_adjustment: 0,
                eff_searchsp: options_searchsp as f64,
            }
        }
    };
    let mut adjust_contexts = Vec::with_capacity(contexts.len());
    for (q_idx, ctx) in contexts.iter().enumerate() {
        let batch = &query_batches[q_idx];
        if q_idx + 1 < batch.end {
            let next = context_of(q_idx + 1);
            if let Some(refusal) = refuse {
                if ctx.is_valid && ctx.aa_len > 0 && next.eff_searchsp == 0.0 {
                    // Before the batch's redo: none of its fetches read a header.
                    return Err(crate::algorithm::common::database::HeaderReadError {
                        query: batch.start,
                        error: unified_p_without_search_space(refusal, q_idx),
                        utf8: None,
                        title_reads: Vec::new(),
                    }
                    .into());
                }
            }
            adjust_contexts.push(next);
        } else {
            adjust_contexts.push(context_of(q_idx));
        }
    }
    Ok(adjust_contexts)
}

#[derive(Clone)]
struct CachedSubjectRangeData {
    subject: BlastCompoSequenceData,
    subject_maybe_biased: bool,
    // `read_aa_composition(subject.data())`, read when the
    // entry was made.
    x_composition: Option<std::sync::Arc<XSubjectComposition>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1614-1646
//
// Protein redo windows use the full subject range for each query. SEG masking
// depends only on the subject range and SEG parameters after the near-identical
// branch has decided masking is required, so Rust reuses that masked range for
// later query-subject pairs while preserving the same returned sequence bytes.
pub(crate) struct BlastpKappaSubjectRangeCache {
    masked_ranges: HashMap<(i32, i32, i32), CachedSubjectRangeData>,
    // Graduated 2026-10-22 (was `LOSAT_X_SEGSHARE`): when set, used instead of `masked_ranges`.
    x_shared: Option<std::sync::Arc<XSharedSubjectRanges>>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1629-1637
// SEG of the subject range depends on the range and the SEG parameters only, so one masked copy
// serves every query that reaches the same range. Placement of an identical result; no value is
// computed differently.
/// Graduated 2026-10-22 (was `LOSAT_X_SEGSHARE`): one store of masked subject ranges for all
/// queries of a search.
///
/// The serial redo keeps a single `BlastpKappaSubjectRangeCache` for the whole
/// search, so a subject is SEG-masked once however many queries hit it. The
/// query-parallel redo creates one cache per query, in which a subject range
/// is never looked up twice, so every (query, subject) pair runs SEG again.
/// Those per-query caches read and write this store. The one-thread redo of
/// the block search keeps the blocks' store from block to block where its
/// subjects are those encoded for all blocks (`x_query_blocks::x_block_seg_store`). The
/// entry for a key is the value the serial cache holds for it: the masked
/// range and its bias flag depend on the subject range only. Two threads may
/// compute the same entry at the same time; they then store equal values.
pub(crate) struct XSharedSubjectRanges {
    shards: Vec<std::sync::Mutex<HashMap<(i32, i32, i32), CachedSubjectRangeData>>>,
}

impl XSharedSubjectRanges {
    pub(crate) fn new() -> Self {
        Self {
            shards: (0..64)
                .map(|_| std::sync::Mutex::new(HashMap::new()))
                .collect(),
        }
    }

    fn shard(
        &self,
        key: &(i32, i32, i32),
    ) -> std::sync::MutexGuard<'_, HashMap<(i32, i32, i32), CachedSubjectRangeData>> {
        let shard = &self.shards[key.0 as u32 as usize % self.shards.len()];
        // Nothing but the map operation runs under the lock.
        shard
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner())
    }

    // No NCBI counterpart: memory only (LOSAT_X_BPSEGWIN, `x_bp_seg_win`); it changes
    // no value NCBI computes.
    /// Drops every entry and the maps' room. A range looked up afterwards is masked again, into the
    /// value it had: an entry depends on its key only.
    pub(crate) fn x_clear(&self) {
        for shard in &self.shards {
            *shard
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner()) = HashMap::new();
        }
    }

    // No NCBI counterpart: memory only (LOSAT_X_BPSEGWIN's size rule); it changes no value NCBI
    // computes.
    /// The number of entries (LOSAT_X_BPSEGWIN's rule, after a window: the sum of the shards').
    pub(crate) fn x_len(&self) -> usize {
        self.shards
            .iter()
            .map(|shard| {
                shard
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner())
                    .len()
            })
            .sum()
    }
}

/// EXPERIMENT (LOSAT_X_BPSEGWIN): the built-in number of entries (about 320 bytes each) from
/// which the store is emptied after a window.
const X_SEG_WIN_ENTRIES: usize = 4_096;

// No NCBI counterpart: reads the LOSAT_X_BPSEGWIN switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPSEGWIN): the number of entries from which the search's store of
/// SEG-masked subject ranges (`XSharedSubjectRanges`) is emptied after a window, or `None` (off).
/// A strict-set value switch: the integer given, 2 or more, else 4,096; `None` for `0`. Where the
/// redo runs in windows, a store that holds that many entries once a window is computed is
/// emptied, so it holds about the ranges of one window, not those of the whole search; a smaller
/// store is kept for the next windows (masking its ranges again costs more time than it saves
/// memory: bp_wssv at 8 threads, 86 entries). A later window that reaches a range again masks it
/// again: an entry is a function of its key (the subject, the range and the SEG parameters), so
/// it gets the same bytes and bias flag.
pub(crate) fn x_bp_seg_win() -> Option<usize> {
    static MIN: OnceLock<Option<usize>> = OnceLock::new();
    *MIN.get_or_init(|| {
        crate::utils::xswitch::value("LOSAT_X_BPSEGWIN", 2..=usize::MAX, X_SEG_WIN_ENTRIES)
    })
}

impl BlastpKappaSubjectRangeCache {
    // No NCBI counterpart: constructor that selects the search-wide store (see
    // XSharedSubjectRanges); it does not change any value NCBI computes.
    /// A cache backed by `shared`, a store of more than one cache (`None`: the cache keeps its own
    /// map).
    pub(crate) fn new(shared: Option<std::sync::Arc<XSharedSubjectRanges>>) -> Self {
        Self {
            masked_ranges: HashMap::new(),
            x_shared: shared,
        }
    }

    fn get(
        &self,
        matching_seq: &BlastCompoMatchingSequence<'_>,
        subject_range: &BlastCompoSequenceRange,
    ) -> Option<CachedSubjectRangeData> {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
        // Dispatch point: with a shared store, a cached SEG result is read from it
        // instead of the per-query map. The stored value is what s_DoSegSequenceData left in the
        // buffer for this range.
        if let Some(shared) = &self.x_shared {
            let key = (matching_seq.index, subject_range.begin, subject_range.end);
            return shared.shard(&key).get(&key).cloned();
        }
        self.masked_ranges
            .get(&(matching_seq.index, subject_range.begin, subject_range.end))
            .cloned()
    }

    fn insert(
        &mut self,
        matching_seq: &BlastCompoMatchingSequence<'_>,
        subject_range: &BlastCompoSequenceRange,
        subject: &BlastCompoSequenceData,
        subject_maybe_biased: bool,
        x_composition: Option<std::sync::Arc<XSubjectComposition>>,
    ) {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
        // Dispatch point: with the switch on, the SEG-masked range is stored in the shared store
        // instead of the per-query map (same key, same value).
        if let Some(shared) = &self.x_shared {
            let key = (matching_seq.index, subject_range.begin, subject_range.end);
            let value = CachedSubjectRangeData {
                subject: subject.clone(),
                subject_maybe_biased,
                x_composition,
            };
            shared.shard(&key).insert(key, value);
            return;
        }
        self.masked_ranges.insert(
            (matching_seq.index, subject_range.begin, subject_range.end),
            CachedSubjectRangeData {
                subject: subject.clone(),
                subject_maybe_biased,
                x_composition,
            },
        );
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1910-1916
struct BlastpRedoAlignmentContext<'a> {
    preliminary_hits: &'a [BlastpPreliminaryHsp],
    matrix: ScoringMatrix,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1719-1720
    gap_scratch: NonNull<GapAlignScratch>,
    subject_range_cache: NonNull<BlastpKappaSubjectRangeCache>,
    _scratch_lifetime: std::marker::PhantomData<&'a mut ()>,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1413-1416
fn blastp_subject_seg_params() -> SegParams {
    SegParams::new(10, 1.8, 2.1)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803
fn build_incoming_alignment_list(
    preliminary_hits: &[BlastpPreliminaryHsp],
    score_divisor: f64,
) -> Option<Box<BlastCompoAlignment>> {
    let mut head: Option<Box<BlastCompoAlignment>> = None;
    let mut tail = &mut head;
    for (prelim_index, preliminary_hsp) in preliminary_hits.iter().enumerate() {
        let node = blast_compo_alignment_new(
            ((preliminary_hsp.raw_score as f64) * score_divisor) as i32,
            EMatrixAdjustRule::DontAdjustMatrix,
            preliminary_hsp.query_start,
            preliminary_hsp.query_end,
            preliminary_hsp.query_context,
            preliminary_hsp.subject_start,
            preliminary_hsp.subject_end,
            preliminary_hsp.subject_frame,
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:770-803
            //
            // Rust stores the stable preliminary-hit index instead of a raw
            // `BlastHSP*`, preserving the same lookup without borrowing the
            // caller-owned slice across the redo-alignment list.
            Some(BlastCompoAlignmentContext::PreliminaryHspIndex(
                prelim_index,
            )),
        );
        *tail = Some(node);
        tail = &mut tail.as_mut().expect("incoming align inserted").next;
    }
    head
}

#[derive(Debug, Clone)]
struct RedoneBlastpHit {
    hit: Hit,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1755-1770
// Kappa consumes raw signed state before any report Hit conversion.
struct RedoneBlastpTraceback {
    score: i32,
    edit_script: Vec<crate::common::GapEditOp>,
    query_start: i32,
    query_end: i32,
    subject_start: i32,
    subject_end: i32,
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1748-1770
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1902-1945
fn redo_preliminary_blastp_hit(
    preliminary_hsp: &BlastpPreliminaryHsp,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    matrix: ScoringMatrix,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    gap_open: i32,
    gap_extend: i32,
    x_drop: i32,
    gap_scratch: &mut GapAlignScratch,
) -> Option<RedoneBlastpTraceback> {
    let q_start = preliminary_hsp.gapped_query_start - query_range.begin;
    let s_start = preliminary_hsp.gapped_subject_start - subject_range.begin;
    if q_start < 0 || s_start < 0 {
        return None;
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1914-1944
    let aligned = blast_gapped_alignment_with_traceback_with_scratch(
        query_data.data(),
        subject_data.data(),
        q_start as usize,
        s_start as usize,
        matrix,
        adjusted_matrix,
        gap_open,
        gap_extend,
        x_drop,
        gap_scratch,
        None,
    )?;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1760-1764,1943-1946
    // Neither an empty edit script nor fence=true overrides successful status.
    let query_start = aligned.query_start + query_range.begin;
    let query_end = aligned.query_stop + query_range.begin;
    let subject_start = aligned.subject_start + subject_range.begin;
    let subject_end = aligned.subject_stop + subject_range.begin;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1757-1764
    if crate::utils::xdiag::var_os("LOSAT_DEBUG_BLASTP_WSSV").is_some()
        && query_start == 78
        && query_end == 789
        && subject_start == 0
        && subject_end == 665
    {
        eprintln!(
            "[BLASTP_WSSV] prelim_gapped_start=({}, {}) score={} template=({}, {}) qidx={} sidx={}",
            preliminary_hsp.gapped_query_start,
            preliminary_hsp.gapped_subject_start,
            aligned.score,
            preliminary_hsp.query_start,
            preliminary_hsp.subject_start,
            preliminary_hsp.q_idx,
            preliminary_hsp.s_idx,
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1767-1770
    Some(RedoneBlastpTraceback {
        score: aligned.score,
        edit_script: aligned.edit_script,
        query_start,
        query_end,
        subject_start,
        subject_end,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1910-1916
fn blastp_redo_alignment_context<'a>(
    params: &'a BlastRedoAlignParams,
) -> &'a BlastpRedoAlignmentContext<'a> {
    let context = params
        .gapping_params
        .context
        .get()
        .expect("NCBI blast_kappa.c stores redo callback context in gapping_params->context");
    // Safety: `postprocess_preliminary_hits` stores a pointer to a
    // `BlastpRedoAlignmentContext` in `gapping_params.context` for the full
    // duration of the `blast_redo_one_match` call that reaches this callback.
    unsafe { &*context.cast::<BlastpRedoAlignmentContext<'a>>().as_ptr() }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2389-2394
fn blastp_get_range_callback<'seq>(
    matching_seq: &BlastCompoMatchingSequence<'seq>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    subject_maybe_biased: bool,
    params: &BlastRedoAlignParams,
) -> Result<BlastRedoRangeResult> {
    let redo_context = blastp_redo_alignment_context(params);
    // Safety: `postprocess_preliminary_hits` stores a pointer to its
    // thread-local subject range cache for the synchronous duration of
    // `Blast_RedoOneMatch`. NCBI passes the same state through
    // `gapping_params->context` and invokes get_range/redo_one_alignment
    // serially for one match.
    let subject_range_cache = unsafe { redo_context.subject_range_cache.as_ptr().as_mut() }
        .expect("blastp redo subject range cache pointer is valid during callback");
    Ok(build_subject_range_data(
        matching_seq,
        subject_range,
        orig_query,
        query_range,
        query_words,
        align,
        should_test_identical,
        subject_maybe_biased,
        params,
        Some(subject_range_cache),
    ))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1896-1945
fn blastp_redo_one_alignment_callback(
    incoming_align: &BlastCompoAlignment,
    matrix_adjust_rule: EMatrixAdjustRule,
    adjusted_matrix: Option<&AdjustedProteinMatrix>,
    query_data: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    _ccat_query_length: i32,
    subject_data: &BlastCompoSequenceData,
    subject_range: &BlastCompoSequenceRange,
    _full_subject_length: i32,
    params: &BlastRedoAlignParams,
) -> Result<Option<Box<BlastCompoAlignment>>> {
    let Some(BlastCompoAlignmentContext::PreliminaryHspIndex(prelim_index)) =
        incoming_align.context.as_ref()
    else {
        return Ok(None);
    };
    let redo_context = blastp_redo_alignment_context(params);
    let preliminary_hsp = &redo_context.preliminary_hits[*prelim_index];
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1912-1914
    // Safety: `postprocess_preliminary_hits` stores a pointer to its
    // thread-local `GapAlignScratch` for the synchronous duration of
    // `Blast_RedoOneMatch`. The callback sequence is serial for a single
    // match, mirroring NCBI's `context->gap_align` reuse.
    let gap_scratch = unsafe { redo_context.gap_scratch.as_ptr().as_mut() }
        .expect("blastp redo gap scratch pointer is valid during callback");
    let Some(mut redone_hit) = redo_preliminary_blastp_hit(
        preliminary_hsp,
        query_data,
        query_range,
        subject_data,
        subject_range,
        redo_context.matrix,
        adjusted_matrix,
        params.gapping_params.gap_open,
        params.gapping_params.gap_extend,
        params.gapping_params.x_dropoff,
        gap_scratch,
    ) else {
        return Ok(None);
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1767-1772
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1767-1772
    let edit_script_context = Some(BlastCompoAlignmentContext::EditScript(std::mem::take(
        &mut redone_hit.edit_script,
    )));
    Ok(Some(blast_compo_alignment_new(
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1767
        // Wraps the redone score, matrix adjust rule, and alignment range via `BlastCompo_AlignmentNew`.
        redone_hit.score,
        matrix_adjust_rule,
        redone_hit.query_start,
        redone_hit.query_end,
        query_range.context,
        redone_hit.subject_start,
        redone_hit.subject_end,
        subject_range.context,
        edit_script_context,
    )))
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2389-2394
const BLASTP_REDO_ALIGN_CALLBACKS: BlastRedoAlignCallbacks = BlastRedoAlignCallbacks {
    calc_lambda: Some(redo_calc_lambda),
    get_range: blastp_get_range_callback,
    redo_one_alignment: blastp_redo_one_alignment_callback,
};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:323-331
fn redone_hit_from_alignment_owned(
    mut align: BlastCompoAlignment,
    template_hit: &BlastpPreliminaryHsp,
) -> Option<RedoneBlastpHit> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:326-331,1598-1643;
    // c++/src/algo/blast/core/blast_gapalign.c:4620-4653
    // The actual BLASTP protein-range getter copies NCBISTDAA 0..27,
    // adds NULLB=0, and SEG writes X=21; it never constructs FENCE_SENTRY=201.
    // Valid seeded starts therefore reach `found_end` as true. Report conversion
    // enforces that domain; it must never silently drop a signed raw payload.
    let query_start = usize::try_from(align.query_start)
        .expect("NCBI full-protein Kappa query start must be non-negative");
    let query_end = usize::try_from(align.query_end)
        .expect("NCBI full-protein Kappa query end must be non-negative");
    let subject_start = usize::try_from(align.match_start)
        .expect("NCBI full-protein Kappa subject start must be non-negative");
    let subject_end = usize::try_from(align.match_end)
        .expect("NCBI full-protein Kappa subject end must be non-negative");
    let gap_info = match align.context.take() {
        Some(BlastCompoAlignmentContext::EditScript(edit_script)) => Some(edit_script),
        None => None,
        Some(BlastCompoAlignmentContext::PreliminaryHspIndex(_)) => return None,
    };

    Some(RedoneBlastpHit {
        hit: Hit {
            identity: 0.0,
            length: 0,
            mismatch: 0,
            gapopen: 0,
            q_start: query_start + 1,
            q_end: query_end,
            s_start: subject_start + 1,
            s_end: subject_end,
            e_value: 0.0,
            bit_score: 0.0,
            num_ident: 0,
            query_frame: template_hit.query_frame,
            query_length: template_hit.query_length,
            q_idx: template_hit.q_idx,
            s_idx: template_hit.s_idx,
            raw_score: align.score,
            // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:125-148
            sort_query_offset: 0,
            sort_query_end: 0,
            sort_subject_offset: 0,
            sort_subject_end: 0,
            has_sort_offsets: false,
            gap_info,
            num_positives: 0,
        },
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:305-354
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3656-3661
fn hsplist_from_distinct_alignments(
    mut alignments: Option<Box<BlastCompoAlignment>>,
    template_hit: &BlastpPreliminaryHsp,
    oid: u32,
    query_index: u32,
) -> BlastpHspList {
    // Room for exactly the alignments of the list (`x_redo_fit.rs`).
    let capacity = super::x_redo_fit::x_alignment_count(alignments.as_deref());
    let mut hsp_list = BlastpHspList {
        oid,
        query_index,
        hsps: Vec::with_capacity(capacity),
        best_evalue: i32::MAX as f64,
    };
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:321-356
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:332-342
    let mut current = alignments.take();
    while let Some(mut align) = current {
        current = align.next.take();
        let comp_adjustment_method = match align.matrix_adjust_rule {
            EMatrixAdjustRule::DontAdjustMatrix => 0,
            EMatrixAdjustRule::CompoScaleOldMatrix => 1,
            _ => 2,
        };
        if let Some(redone_hit) = redone_hit_from_alignment_owned(*align, template_hit) {
            let mut hsp = BlastpHsp::from_hit(redone_hit.hit);
            hsp.comp_adjustment_method = comp_adjustment_method;
            hsp_list.hsps.push(hsp);
        }
    }
    sort_hsplist_by_score(&mut hsp_list);
    hsp_list
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1510-1647
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1682-1704
fn build_query_range_data(
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
) -> BlastCompoSequenceData {
    let query =
        BlastCompoSequenceData::copy_query_range_with_selenocysteine_fix(orig_query, query_range);
    record_blastp_kappa_range_counter(
        &BLASTP_KAPPA_QUERY_RANGE_BYTES_COPIED,
        query.length.max(0) as usize,
    );
    query
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1614-1646
fn build_subject_range_data_direct(
    matching_seq: &BlastCompoMatchingSequence<'_>,
    subject_range: &BlastCompoSequenceRange,
) -> BlastCompoSequenceData {
    let begin = usize::try_from(subject_range.begin)
        .expect("NCBI blast_kappa.c requires non-negative subject_range.begin");
    let end = usize::try_from(subject_range.end)
        .expect("NCBI blast_kappa.c requires non-negative subject_range.end");
    let subject_len = end.saturating_sub(begin);
    let mut buffer = vec![0u8; subject_len + 2];
    buffer[1..1 + subject_len].copy_from_slice(&matching_seq.data[begin..end]);
    record_blastp_kappa_range_counter(&BLASTP_KAPPA_SUBJECT_RANGE_BYTES_COPIED, subject_len);
    BlastCompoSequenceData {
        buffer,
        data_offset: 1,
        length: i32::try_from(subject_len)
            .expect("NCBI blast_kappa.c subject range length must fit in Int4"),
    }
}

fn build_subject_range_data(
    matching_seq: &BlastCompoMatchingSequence<'_>,
    subject_range: &BlastCompoSequenceRange,
    orig_query: &BlastCompoSequenceData,
    query_range: &BlastCompoSequenceRange,
    query_words: Option<&[u64]>,
    align: &BlastCompoAlignment,
    should_test_identical: bool,
    subject_maybe_biased_in: bool,
    params: &BlastRedoAlignParams,
    mut subject_range_cache: Option<&mut BlastpKappaSubjectRangeCache>,
) -> BlastRedoRangeResult {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1677-1709
    let query = build_query_range_data(orig_query, query_range);
    let mut subject_for_near_identical: Option<BlastCompoSequenceData> = None;
    let mut subject_maybe_biased = subject_maybe_biased_in;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1626-1636
    let should_apply_seg = params.uses_composition_based_stats()
        && subject_maybe_biased
        && (!should_test_identical || {
            let subject = subject_for_near_identical.get_or_insert_with(|| {
                build_subject_range_data_direct(matching_seq, subject_range)
            });
            !test_near_identical(
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1267-1277
                //
                // The Rust subject buffer already begins at `subject_range.begin`,
                // so the NCBI `seqOffset` argument remains the original range
                // origin.
                subject,
                subject_range.begin,
                &query,
                query_range.begin,
                query_words,
                align,
            )
        });

    if should_apply_seg {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1635-1636
        // A cached SEG result of this range (this cache or the shared store) replaces the call.
        if let Some(cache) = subject_range_cache.as_mut() {
            if let Some(cached) = cache.get(matching_seq, subject_range) {
                record_blastp_kappa_range_event(&BLASTP_KAPPA_SUBJECT_SEG_CACHE_HITS);
                record_blastp_kappa_range_counter(
                    &BLASTP_KAPPA_SUBJECT_RANGE_BYTES_COPIED,
                    cached.subject.length.max(0) as usize,
                );
                return BlastRedoRangeResult {
                    query,
                    subject: cached.subject,
                    subject_maybe_biased: cached.subject_maybe_biased,
                    x_subject_composition: cached.x_composition,
                };
            }
            record_blastp_kappa_range_event(&BLASTP_KAPPA_SUBJECT_SEG_CACHE_MISSES);
        }
        let mut subject = subject_for_near_identical
            .take()
            .unwrap_or_else(|| build_subject_range_data_direct(matching_seq, subject_range));
        let seg_params = blastp_subject_seg_params();
        let masker = SegMasker::with_params(&seg_params);
        record_blastp_kappa_range_event(&BLASTP_KAPPA_SUBJECT_SEG_CALLS);
        record_blastp_kappa_range_counter(
            &BLASTP_KAPPA_SUBJECT_SEG_RANGE_BYTES,
            subject.length.max(0) as usize,
        );
        let intervals = masker.mask_sequence(subject.data());
        subject_maybe_biased = !intervals.is_empty();
        for interval in &intervals {
            record_blastp_kappa_range_counter(
                &BLASTP_KAPPA_SUBJECT_SEG_MASKED_RESIDUES,
                interval.end.saturating_sub(interval.start),
            );
            for pos in interval.start..interval.end {
                subject.buffer[subject.data_offset + pos] = ncbistdaa::X;
            }
        }
        // The composition of the masked bytes, read once here
        // and kept with them.
        // NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:934-936
        let x_composition = Some(std::sync::Arc::new(XSubjectComposition(
            read_aa_composition(subject.data()),
        )));
        if let Some(cache) = subject_range_cache.as_mut() {
            cache.insert(
                matching_seq,
                subject_range,
                &subject,
                subject_maybe_biased,
                x_composition.clone(),
            );
        }
        return BlastRedoRangeResult {
            query,
            subject,
            subject_maybe_biased,
            x_subject_composition: x_composition,
        };
    }

    let subject = subject_for_near_identical
        .take()
        .unwrap_or_else(|| build_subject_range_data_direct(matching_seq, subject_range));
    BlastRedoRangeResult {
        query,
        subject,
        subject_maybe_biased,
        x_subject_composition: None,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:223-260
fn reap_contained_hits(hsp_list: &mut BlastpHspList) {
    if hsp_list.hsps.len() < 2 {
        return;
    }

    let mut keep = vec![true; hsp_list.hsps.len()];
    for iread in 1..hsp_list.hsps.len() {
        if !keep[iread] {
            continue;
        }
        let hsp1 = &hsp_list.hsps[iread];
        let hsp1_query_start = hsp1.q_start.saturating_sub(1) as i32;
        let hsp1_query_end = hsp1.q_end as i32;
        let hsp1_subject_start = hsp1.s_start.saturating_sub(1) as i32;
        let hsp1_subject_end = hsp1.s_end as i32;
        for iread_back in 0..iread {
            if !keep[iread_back] {
                continue;
            }
            let hsp2 = &hsp_list.hsps[iread_back];
            let hsp2_query_start = hsp2.q_start.saturating_sub(1) as i32;
            let hsp2_query_end = hsp2.q_end as i32;
            let hsp2_subject_start = hsp2.s_start.saturating_sub(1) as i32;
            let hsp2_subject_end = hsp2.s_end as i32;
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:248-259
            if hsp_list.hsps[iread_back].query_frame != hsp_list.hsps[iread].query_frame
                || hsp_list.hsps[iread_back].subject_frame != hsp_list.hsps[iread].subject_frame
            {
                continue;
            }
            let start_contained = hsp2_query_start <= hsp1_query_start
                && hsp2_query_end >= hsp1_query_start
                && hsp2_subject_start <= hsp1_subject_start
                && hsp2_subject_end >= hsp1_subject_start;
            let end_contained = hsp2_query_start <= hsp1_query_end
                && hsp2_query_end >= hsp1_query_end
                && hsp2_subject_start <= hsp1_subject_end
                && hsp2_subject_end >= hsp1_subject_end;
            if start_contained
                && end_contained
                && hsp_list.hsps[iread].raw_score <= hsp_list.hsps[iread_back].raw_score
            {
                keep[iread] = false;
                break;
            }
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:266-273
    let mut read_index = 0usize;
    hsp_list.hsps.retain(|_| {
        let keep_hit = keep[read_index];
        read_index += 1;
        keep_hit
    });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:101-116
fn normalize_scores(hsp_list: &mut BlastpHspList, lambda: f64, log_k: f64, score_divisor: f64) {
    for hsp in &mut hsp_list.hsps {
        hsp.raw_score = blast_nint((hsp.raw_score as f64) / score_divisor);
        hsp.bit_score =
            (hsp.raw_score as f64 * lambda * score_divisor - log_k) / std::f64::consts::LN_2;
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/ncbi_math.c:437-441
fn blast_nint(x: f64) -> i32 {
    (x + if x >= 0.0 { 0.5 } else { -0.5 }) as i32
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2120-2133
fn scaled_karlin_params(params: &KarlinParams, score_divisor: f64) -> KarlinParams {
    if score_divisor == 1.0 {
        return *params;
    }
    KarlinParams {
        lambda: params.lambda / score_divisor,
        k: params.k,
        h: params.h,
        alpha: params.alpha,
        beta: params.beta,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2308-2336
pub(crate) fn build_query_workspace(
    query: &[u8],
    query_length: i32,
    length_adjustment: i32,
    effective_space: f64,
    evalue_adjust_context: EvalueAdjustContext,
) -> BlastpKappaQueryWorkspace {
    let redo_query_info = BlastCompoQueryInfo {
        origin: 0,
        seq: BlastCompoSequenceData::from_ncbistdaa(query),
        composition: read_aa_composition(query),
        eff_search_space: effective_space,
        words: Some(build_query_word_hashes(query)),
    };
    BlastpKappaQueryWorkspace {
        redo_query_info,
        query_length,
        length_adjustment,
        evalue_adjust_context,
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-824
//
// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:518-529
fn recompute_num_identities(
    query_nomask: &[u8],
    subject: &[u8],
    hsps: &mut [BlastpHsp],
    matrix: ScoringMatrix,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-818
    if matrix == ScoringMatrix::Blosum62 {
        recompute_num_identities_impl::<true>(query_nomask, subject, hsps, matrix);
    } else {
        recompute_num_identities_impl::<false>(query_nomask, subject, hsps, matrix);
    }
}

/// NCBI's engine positives of the HSPs of one match (`BlastHSP::num_positives`, which the
/// Seq-aligns show): `s_ComputeNumIdentities` counts them after the redo with the score block
/// whose matrix the match's composition adjustments have overwritten (`Blast_AdjustScores`
/// writes its `matrix` argument, the thread-local `sbp->matrix->data`), so a pair is positive
/// by `adjusted`, the matrix of the match's last adjustment. `recompute_num_identities` counts
/// them with the search matrix for the other reports; this replaces only `num_positives`, and
/// leaves an HSP whose letters it cannot read as it was. With the Seq-aligns only
/// (`postprocess_preliminary_hits`' `engine_positives`).
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3500-3502,3625-3635
/// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1240-1241
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3691-3701
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:770-833
fn recount_engine_positives(
    query_nomask: &[u8],
    subject: &[u8],
    hsps: &mut [BlastpHsp],
    adjusted: &AdjustedProteinMatrix,
) {
    for hsp in hsps {
        let Some(q_offset) = hsp.q_start.checked_sub(1) else {
            continue;
        };
        let Some(s_offset) = hsp.s_start.min(hsp.s_end).checked_sub(1) else {
            continue;
        };
        let implicit = [GapEditOp::Sub(hsp.q_end.saturating_sub(q_offset) as u32)];
        let ops = hsp.gap_info.as_deref().unwrap_or(&implicit);
        let (mut q_index, mut s_index) = (q_offset, s_offset);
        let mut num_positives = 0usize;
        let mut valid = true;
        for op in ops {
            match *op {
                GapEditOp::Sub(n) => {
                    for _ in 0..n as usize {
                        let (Some(&q), Some(&s)) =
                            (query_nomask.get(q_index), subject.get(s_index))
                        else {
                            valid = false;
                            break;
                        };
                        if q == s || adjusted.is_positive(q, s) {
                            num_positives += 1;
                        }
                        q_index += 1;
                        s_index += 1;
                    }
                }
                GapEditOp::Del(n) => s_index += n as usize,
                GapEditOp::Ins(n) => q_index += n as usize,
            }
            if !valid {
                break;
            }
        }
        if valid {
            hsp.num_positives = num_positives;
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:745-824
fn recompute_num_identities_impl<const BLOSUM62: bool>(
    query_nomask: &[u8],
    subject: &[u8],
    hsps: &mut [BlastpHsp],
    matrix: ScoringMatrix,
) {
    for hsp in hsps {
        let Some(q_offset) = hsp.q_start.checked_sub(1) else {
            continue;
        };
        let Some(s_offset) = hsp.s_start.min(hsp.s_end).checked_sub(1) else {
            continue;
        };

        let mut q_index = q_offset;
        let mut s_index = s_offset;
        let mut num_ident = 0usize;
        let mut num_positives = 0usize;
        let mut gap_letters = 0usize;
        let mut gap_opens = 0usize;
        let mut alignment_len = 0usize;
        let mut valid = true;

        if let Some(gap_info) = hsp.gap_info.as_ref() {
            for op in gap_info {
                match *op {
                    GapEditOp::Sub(n) => {
                        for _ in 0..n as usize {
                            let (Some(&q), Some(&s)) =
                                (query_nomask.get(q_index), subject.get(s_index))
                            else {
                                valid = false;
                                break;
                            };
                            if q == s {
                                num_ident += 1;
                                num_positives += 1;
                            } else if blastp_positive_score::<BLOSUM62>(matrix, q, s) {
                                num_positives += 1;
                            }
                            q_index += 1;
                            s_index += 1;
                            alignment_len += 1;
                        }
                    }
                    GapEditOp::Del(n) => {
                        gap_opens += 1;
                        for _ in 0..n as usize {
                            if subject.get(s_index).is_none() {
                                valid = false;
                                break;
                            }
                            s_index += 1;
                            gap_letters += 1;
                            alignment_len += 1;
                        }
                    }
                    GapEditOp::Ins(n) => {
                        gap_opens += 1;
                        for _ in 0..n as usize {
                            if query_nomask.get(q_index).is_none() {
                                valid = false;
                                break;
                            }
                            q_index += 1;
                            gap_letters += 1;
                            alignment_len += 1;
                        }
                    }
                }
                if !valid {
                    break;
                }
            }

            if !valid {
                continue;
            }

            hsp.length = alignment_len;
            hsp.num_ident = num_ident;
            hsp.num_positives = num_positives;
            hsp.gapopen = gap_opens;
            hsp.mismatch = alignment_len.saturating_sub(num_ident + gap_letters) as i64;
            hsp.identity = if alignment_len > 0 {
                (num_ident as f64 * 100.0) / (alignment_len as f64)
            } else {
                0.0
            };
            continue;
        }

        let q_end = hsp.q_end;
        let s_start = hsp.s_start.min(hsp.s_end) - 1;
        let s_end = hsp.s_start.max(hsp.s_end);
        if q_end < q_offset || s_end < s_start {
            continue;
        }
        let q_span = q_end.saturating_sub(q_offset);
        let s_span = s_end.saturating_sub(s_start);
        if q_span != s_span || q_end > query_nomask.len() || s_end > subject.len() {
            continue;
        }

        for (&q, &s) in query_nomask[q_offset..q_end]
            .iter()
            .zip(subject[s_start..s_end].iter())
        {
            if q == s {
                num_ident += 1;
                num_positives += 1;
            } else if blastp_positive_score::<BLOSUM62>(matrix, q, s) {
                num_positives += 1;
            }
        }
        hsp.length = q_span;
        hsp.num_ident = num_ident;
        hsp.num_positives = num_positives;
        hsp.mismatch = q_span.saturating_sub(num_ident) as i64;
        hsp.gapopen = 0;
        hsp.identity = if q_span > 0 {
            (num_ident as f64 * 100.0) / (q_span as f64)
        } else {
            0.0
        };
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:147-158
/// `s_AdjustEvaluesForComposition`: the E-values of a match's HSPs combined with the composition
/// P-value of the pair. `query_context` is what NCBI's `query_context` holds
/// ([`evalue_adjust_contexts`]: the next query's values, or the query's own for the last query
/// of a batch); its effective search space is above 0 (a batch where it is 0 is refused), so the
/// scale is finite and above 0.
fn adjust_evalues_for_composition(
    hsp_list: &mut BlastpHspList,
    comp_p_value: f64,
    query_context: &EvalueAdjustContext,
    subject_length: i32,
    lambda_ratio: f64,
) {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:205-209
    let _ = lambda_ratio;
    let query_eff = (query_context.query_length - query_context.length_adjustment).max(1) as f64;
    let subject_eff = ((subject_length - query_context.length_adjustment) as f64).max(1.0);
    let dblen_eff = query_context.eff_searchsp / query_eff;
    let db_to_sequence_scale = subject_eff / dblen_eff;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:172-176
    // `stats::e_to_p` / `stats::p_to_e` are `BLAST_KarlinEtoP` / `BLAST_KarlinPtoE` with
    // NCBI's `BLAST_Expm1` / `BLAST_Log1p` (not libm's `exp_m1` / `ln_1p`).
    let mut best_evalue = f64::MAX;
    for hsp in &mut hsp_list.hsps {
        hsp.e_value *= db_to_sequence_scale;
        let align_p_value = e_to_p(hsp.e_value);
        let combined_p_value = blast_overall_p_value(comp_p_value, align_p_value);
        hsp.e_value = p_to_e(combined_p_value);
        hsp.e_value /= db_to_sequence_scale;
        if hsp.e_value < best_evalue {
            best_evalue = hsp.e_value;
        }
    }
    hsp_list.best_evalue = best_evalue;
}

/// `Blast_HSPListGetEvalues` of the redo (no linking): the e-values from the scaled `kbp_gap`,
/// with Spouge's finite-size correction when the search has a Gumbel block.
#[inline]
fn get_evalues(
    hsp_list: &mut BlastpHspList,
    eval_params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    query_workspace: &BlastpKappaQueryWorkspace,
    subject_length: i32,
) {
    match gumbel_params {
        Some(gumbel_params) => {
            for hsp in &mut hsp_list.hsps {
                hsp.e_value = blast_spouge_stoe(
                    hsp.raw_score,
                    eval_params,
                    gumbel_params,
                    query_workspace.query_length,
                    subject_length,
                );
            }
        }
        None => {
            let search_space = query_workspace.redo_query_info.eff_search_space;
            for hsp in &mut hsp_list.hsps {
                hsp.e_value = search_space
                    * (-eval_params.lambda * hsp.raw_score as f64 + eval_params.k.glibc_ln())
                        .glibc_exp();
            }
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:395-457
fn evaluate_and_purge(
    hsp_list: &mut BlastpHspList,
    params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    query_workspace: &BlastpKappaQueryWorkspace,
    subject_length: i32,
    expect_value: f64,
    pvalue_for_this_pair: Option<f64>,
    lambda_ratio: Option<f64>,
    score_divisor: f64,
    do_link_hsps: bool,
) -> (i32, f64) {
    let mut best_evalue = f64::MAX;
    let mut best_score = 0;

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:414-427
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:429-435
    let eval_params = scaled_karlin_params(params, score_divisor);
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:411-415
    // psiblast's `-sum_stats true`: the linked sets' e-values from the scaled `kbp_gap`, with the
    // traceback's linking parameters (`EvenGapLinkParams::traceback`); `BLAST_LinkHsps` sorts the
    // list by score and sets the best e-value.
    if do_link_hsps {
        super::hsp::link_hsp_list(
            hsp_list,
            crate::core::link_hsps::LinkContext {
                lambda: eval_params.lambda,
                log_k: eval_params.k.glibc_ln(),
                query_length: query_workspace.query_length,
                length_adjustment: query_workspace.length_adjustment,
                eff_searchsp: query_workspace.redo_query_info.eff_search_space as i64,
            },
            subject_length,
            &crate::core::link_hsps::EvenGapLinkParams::traceback(),
        );
    } else {
        get_evalues(
            hsp_list,
            &eval_params,
            gumbel_params,
            query_workspace,
            subject_length,
        );
        update_best_evalue(hsp_list);
    }
    if let Some(comp_p_value) = pvalue_for_this_pair.filter(|p| (0.0..=1.0).contains(p)) {
        adjust_evalues_for_composition(
            hsp_list,
            comp_p_value,
            &query_workspace.evalue_adjust_context,
            subject_length,
            lambda_ratio.unwrap_or(1.0),
        );
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:439-446
    reap_hsplist_by_evalue(hsp_list, expect_value);
    if let Some(first_hsp) = hsp_list.hsps.first() {
        best_evalue = hsp_list.best_evalue;
        best_score = first_hsp.raw_score;
    }

    (best_score, best_evalue)
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3658-3710
pub(crate) fn postprocess_preliminary_hits(
    preliminary_hits: &[BlastpPreliminaryHsp],
    query_workspace: &BlastpKappaQueryWorkspace,
    composition_workspace: &mut BlastCompositionWorkspace,
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    karlin_params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    redo_align_params: &BlastRedoAlignParams,
    expect_value: f64,
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1719-1720
    gap_scratch: &mut GapAlignScratch,
    subject_range_cache: &mut BlastpKappaSubjectRangeCache,
    // The Seq-aligns (outfmt 8, 12) are requested (`ReportNeeds::seqalign`): the HSPs get
    // NCBI's engine positives (`recount_engine_positives`).
    engine_positives: bool,
) -> Result<BlastpPostprocessResult> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3635-3640
    // The X-drop redo of the match (`Blast_RedoOneMatch`).
    postprocess_with(
        preliminary_hits,
        query_workspace,
        composition_workspace,
        query_nomask,
        subject,
        matrix,
        karlin_params,
        gumbel_params,
        redo_align_params,
        expect_value,
        gap_scratch,
        subject_range_cache,
        engine_positives,
        |incoming_aligns, matching_seq, lambda, composition_workspace, matrix_state| {
            blast_redo_one_match_with_workspace_queries_and_matrix(
                incoming_aligns,
                redo_align_params,
                matching_seq,
                std::slice::from_ref(&query_workspace.redo_query_info),
                lambda,
                &BLASTP_REDO_ALIGN_CALLBACKS,
                composition_workspace,
                matrix_state,
            )
        },
    )
}

/// The redo of one match with `-use_sw_tback` (item X4): `postprocess_preliminary_hits` with
/// NCBI's Smith-Waterman redo (`smith_waterman::redo_one_match_smith_waterman`) in place of the
/// X-drop redo; the steps after it are the same. `would_insert(evalue, score)` answers the test of
/// the query's heap of redone matches (`BlastCompo_HeapWouldInsert` with the match's subject
/// index): the one-thread loop asks the heap (`SwHeapAnswer::Known`); a job of the parallel redo
/// may leave it open (`Deferred`: the match ends as on FALSE, and `kept` gets the test's state),
/// and the continuation of such a match resumes from that state (`resume`; item X4-mt). The errors
/// `SwRedoFailed` and `SwRedoEndless` are NCBI's failed and endless redos.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3602-3605
// Out of line: the one-thread loop calls it where `-use_sw_tback` asks for it, and the X-drop
// redo's function keeps its own body.
#[inline(never)]
pub(crate) fn postprocess_preliminary_hits_smith_waterman(
    preliminary_hits: &[BlastpPreliminaryHsp],
    query_workspace: &BlastpKappaQueryWorkspace,
    composition_workspace: &mut BlastCompositionWorkspace,
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    karlin_params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    redo_align_params: &BlastRedoAlignParams,
    expect_value: f64,
    gap_scratch: &mut GapAlignScratch,
    subject_range_cache: &mut BlastpKappaSubjectRangeCache,
    engine_positives: bool,
    sw_scratch: &mut smith_waterman::SwMatchScratch,
    resume: Option<&smith_waterman::SwMatchResume>,
    would_insert: &mut dyn FnMut(f64, i32) -> smith_waterman::SwHeapAnswer,
    kept: &mut Option<Box<smith_waterman::SwMatchResume>>,
) -> Result<BlastpPostprocessResult> {
    postprocess_with(
        preliminary_hits,
        query_workspace,
        composition_workspace,
        query_nomask,
        subject,
        matrix,
        karlin_params,
        gumbel_params,
        redo_align_params,
        expect_value,
        gap_scratch,
        subject_range_cache,
        engine_positives,
        |incoming_aligns, matching_seq, lambda, composition_workspace, matrix_state| {
            // BLASTP's match has the HSPs of one query: its heap is the query's
            // (`significantMatches[query_index]`).
            smith_waterman::redo_one_match_smith_waterman(
                incoming_aligns,
                redo_align_params,
                matching_seq,
                std::slice::from_ref(&query_workspace.redo_query_info),
                lambda,
                &BLASTP_REDO_ALIGN_CALLBACKS,
                composition_workspace,
                matrix_state,
                expect_value,
                resume,
                &mut |_query_index, evalue, score| would_insert(evalue, score),
                kept,
                sw_scratch,
            )
        },
    )
}

/// The redo of one match and the steps after it (`postprocess_preliminary_hits`): `redo` is the
/// match's redo (`Blast_RedoOneMatch` or `Blast_RedoOneMatchSmithWaterman`) of the incoming
/// alignments, the matching sequence, the redo's scaled Lambda, the composition workspace and the
/// matrix state. One source with one instance per redo: the X-drop redo's instance is the body of
/// `postprocess_preliminary_hits` (its only caller), the Smith-Waterman redo's that of
/// `postprocess_preliminary_hits_smith_waterman`.
#[inline(always)]
fn postprocess_with<R>(
    preliminary_hits: &[BlastpPreliminaryHsp],
    query_workspace: &BlastpKappaQueryWorkspace,
    composition_workspace: &mut BlastCompositionWorkspace,
    query_nomask: &[u8],
    subject: &[u8],
    matrix: ScoringMatrix,
    karlin_params: &KarlinParams,
    gumbel_params: Option<&BlastGumbelBlk>,
    redo_align_params: &BlastRedoAlignParams,
    expect_value: f64,
    gap_scratch: &mut GapAlignScratch,
    subject_range_cache: &mut BlastpKappaSubjectRangeCache,
    engine_positives: bool,
    redo: R,
) -> Result<BlastpPostprocessResult>
where
    R: FnOnce(
        &Option<Box<BlastCompoAlignment>>,
        &BlastCompoMatchingSequence<'_>,
        f64,
        &mut BlastCompositionWorkspace,
        &mut Option<AdjustedProteinMatrix>,
    ) -> Result<
        crate::core::composition_adjustment::redo_alignment::BlastRedoOneMatchResult,
    >,
{
    if preliminary_hits.is_empty() {
        return Ok(BlastpPostprocessResult {
            hsp_list: None,
            best_score: 0,
            best_evalue: f64::MAX,
        });
    }

    let incoming_aligns =
        build_incoming_alignment_list(preliminary_hits, redo_align_params.score_divisor);
    let matching_seq = BlastCompoMatchingSequence::new(preliminary_hits[0].s_idx as i32, subject);
    if let Some(target) = blastp_trace_hsp_target() {
        if blastp_trace_matches_pair(target, preliminary_hits[0].q_idx, preliminary_hits[0].s_idx) {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3580-3581
            eprintln!(
                "[BLASTP_TRACE_HSP] kappa_input q_idx={} s_idx={} preliminary_count={} expect_value={}",
                preliminary_hits[0].q_idx,
                preliminary_hits[0].s_idx,
                preliminary_hits.len(),
                expect_value,
            );
            for (index, hit) in preliminary_hits.iter().enumerate() {
                eprintln!(
                    "[BLASTP_TRACE_HSP] preliminary[{index}] q={}..{} s={}..{} gapped_start=({}, {}) raw_score={}",
                    hit.query_start + 1,
                    hit.query_end,
                    hit.subject_start + 1,
                    hit.subject_end,
                    hit.gapped_query_start,
                    hit.gapped_subject_start,
                    hit.raw_score,
                );
            }
        }
    }
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2117-2133
    let redo_karlin_params = scaled_karlin_params(karlin_params, redo_align_params.score_divisor);
    let redo_context = BlastpRedoAlignmentContext {
        preliminary_hits,
        matrix,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:1719-1720
        gap_scratch: NonNull::from(gap_scratch),
        subject_range_cache: NonNull::from(subject_range_cache),
        _scratch_lifetime: std::marker::PhantomData,
    };
    let previous_redo_context = redo_align_params
        .gapping_params
        .context
        .replace(Some(NonNull::from(&redo_context).cast::<c_void>()));
    // The matrix that the match's adjustments leave (NCBI's thread-local `sbp->matrix->data`,
    // which `Blast_AdjustScores` overwrites), read by `recount_engine_positives`; the call is
    // that of `blast_redo_one_match_with_workspace`, which starts every match with none.
    let mut matrix_state: Option<AdjustedProteinMatrix> = None;
    let distinct_alignments = redo(
        &incoming_aligns,
        &matching_seq,
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:2120-2128
        redo_karlin_params.lambda,
        composition_workspace,
        &mut matrix_state,
    );
    redo_align_params
        .gapping_params
        .context
        .set(previous_redo_context);
    let mut distinct_alignments = distinct_alignments?;
    let query_index = usize::try_from(preliminary_hits[0].query_context)
        .expect("blastp query context is non-negative");
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3658-3667
    let mut hsp_list = distinct_alignments
        .alignments_by_query
        .get_mut(query_index)
        .and_then(Option::take)
        .map(|alignments| {
            hsplist_from_distinct_alignments(
                Some(alignments),
                &preliminary_hits[0],
                matching_seq.index as u32,
                preliminary_hits[0].q_idx,
            )
        })
        .unwrap_or_else(|| BlastpHspList {
            oid: matching_seq.index as u32,
            query_index: preliminary_hits[0].q_idx,
            hsps: Vec::new(),
            best_evalue: i32::MAX as f64,
        });

    if let Some(target) = blastp_trace_hsp_target() {
        if blastp_trace_matches_pair(target, preliminary_hits[0].q_idx, preliminary_hits[0].s_idx) {
            // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3658-3667
            eprintln!(
                "[BLASTP_TRACE_HSP] from_distinct count={} pvalue={:?} lambda_ratio={:?}",
                hsp_list.hsps.len(),
                distinct_alignments.pvalue_for_this_pair,
                distinct_alignments.lambda_ratio,
            );
            for (index, hsp) in hsp_list.hsps.iter().enumerate() {
                eprintln!(
                    "[BLASTP_TRACE_HSP] distinct_hsp[{index}] q={}..{} s={}..{} raw_score={} {}",
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.raw_score,
                    blastp_trace_gap_summary(hsp.gap_info.as_deref()),
                );
            }
        }
    }

    if hsp_list.hsps.len() > 1 {
        reap_contained_hits(&mut hsp_list);
    }
    let pvalue_for_this_pair = distinct_alignments.pvalue_for_this_pair;
    let lambda_ratio = distinct_alignments.lambda_ratio;
    let (best_score, best_evalue) = evaluate_and_purge(
        &mut hsp_list,
        karlin_params,
        gumbel_params,
        query_workspace,
        matching_seq.length,
        expect_value,
        pvalue_for_this_pair,
        lambda_ratio,
        redo_align_params.score_divisor,
        redo_align_params.do_link_hsps,
    );
    if best_evalue > expect_value || hsp_list.hsps.is_empty() {
        return Ok(BlastpPostprocessResult {
            hsp_list: None,
            best_score,
            best_evalue,
        });
    }

    normalize_scores(
        &mut hsp_list,
        redo_karlin_params.lambda,
        redo_align_params.log_k,
        redo_align_params.score_divisor,
    );
    recompute_num_identities(query_nomask, subject, &mut hsp_list.hsps, matrix);
    if engine_positives {
        if let Some(adjusted) = matrix_state.as_ref() {
            recount_engine_positives(query_nomask, subject, &mut hsp_list.hsps, adjusted);
        }
    }
    if let Some(target) = blastp_trace_hsp_target() {
        for (index, hsp) in hsp_list.hsps.iter().enumerate() {
            if blastp_trace_matches_hsp(target, hsp) {
                // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:531,3691-3701
                eprintln!(
                    "[BLASTP_TRACE_HSP] final_hsp[{index}] q={}..{} s={}..{} raw_score={} bit_score={:.12} evalue={:.12e} len={} ident={} positives={} gaps={} gap_opens={} mismatch={} identity={:.6} {}",
                    hsp.q_start,
                    hsp.q_end,
                    hsp.s_start,
                    hsp.s_end,
                    hsp.raw_score,
                    hsp.bit_score,
                    hsp.e_value,
                    hsp.length,
                    hsp.num_ident,
                    hsp.num_positives,
                    hsp.length as i64 - hsp.num_ident as i64 - hsp.mismatch,
                    hsp.gapopen,
                    hsp.mismatch,
                    hsp.identity,
                    blastp_trace_gap_summary(hsp.gap_info.as_deref()),
                );
            }
        }
    }

    // The kept list gives back the room of the HSPs removed above; the same HSPs in the same
    // order (`x_redo_fit.rs`).
    super::x_redo_fit::x_fit_redone_hsps(&mut hsp_list.hsps);
    Ok(BlastpPostprocessResult {
        hsp_list: Some(hsp_list),
        best_score,
        best_evalue,
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_kappa.c:3525-3529,3577-3585
// NCBI reference (598d8ae6): c++/src/algo/blast/composition_adjustment/redo_alignment.c:1176-1198,1239-1244
// Newton lanes, BLASTP: for a group of local matches of one query (those without
// HSPs, or already ended by `BlastCompo_EarlyTermination` with the current heaps, are skipped), the
// input of each match's first `Blast_AdjustScores` call as the redo of `postprocess_preliminary_hits`
// builds it: preliminary hits, incoming alignments, protein windows, the first window with an
// alignment, its near-identical test, `get_range` through `build_subject_range_data` with the subject
// range cache the redo will use (so the redo finds the SEG-masked range there), the compositions and
// `x_newton_input`. The inputs are solved four per vector and kept (`x_newton_lanes::prefetch`).
// No NCBI counterpart: computes ahead what the redo computes; a result is used only by a call with the
// same input bits, and the cache returns the value the redo computes itself. It does not change any
// value NCBI computes.
#[allow(clippy::too_many_arguments)]
pub(crate) fn x_newton_lanes_prefetch_matches(
    local_matches: &[BlastpHspList],
    subjects: super::x_redo_subjects::XRedoSubjects<'_>,
    heaps: &[BlastCompoHeap],
    query_workspace: &BlastpKappaQueryWorkspace,
    composition_workspace: &BlastCompositionWorkspace,
    redo_align_params: &BlastRedoAlignParams,
    subject_range_cache: &mut BlastpKappaSubjectRangeCache,
) {
    if redo_align_params.smith_waterman
        || redo_align_params.position_based
        || !redo_align_params.uses_composition_based_stats()
    {
        return;
    }
    let query_info = &query_workspace.redo_query_info;
    let mut inputs = Vec::new();
    let mut shared_q = None;
    for local_match in local_matches {
        if local_match.hsps.is_empty()
            || blast_compo_early_termination(local_match.best_evalue, heaps)
        {
            continue;
        }
        let preliminary_hits: Vec<BlastpPreliminaryHsp> = local_match
            .hsps
            .iter()
            .map(super::blast_engine::preliminary_hit_from_local_match_hsp)
            .collect();
        // Read from the database, dropped after the match (the
        // redo reads it again); on an error nothing is computed ahead, and the redo of the match
        // reaches the error itself.
        let Ok(subject) = subjects.read(local_match.oid) else {
            continue;
        };
        let subject_raw = subject.residues();
        let incoming_aligns =
            build_incoming_alignment_list(&preliminary_hits, redo_align_params.score_divisor);
        let matching_seq =
            BlastCompoMatchingSequence::new(preliminary_hits[0].s_idx as i32, subject_raw);
        let Ok(windows) =
            windows_from_protein_aligns(&incoming_aligns, query_info, matching_seq.length)
        else {
            continue;
        };
        let Some((window, head)) = windows
            .iter()
            .find_map(|window| window.align.as_deref().map(|head| (window, head)))
        else {
            continue;
        };
        let near_identical = preliminary_test_near_identical(
            query_info,
            window,
            head,
            redo_align_params.near_identical_cutoff,
        );
        let range = build_subject_range_data(
            &matching_seq,
            &window.subject_range,
            &query_info.seq,
            &window.query_range,
            query_info.words.as_deref(),
            head,
            near_identical,
            true,
            redo_align_params,
            Some(&mut *subject_range_cache),
        );
        let subject_composition = x_whole_subject_composition(&range);
        if query_info.composition.num_true_amino_acids == 0
            || subject_composition.num_true_amino_acids == 0
        {
            continue;
        }
        if let Some(input) = x_newton_input(
            &redo_align_params.matrix_info,
            &query_info.composition,
            range.query.length,
            &subject_composition,
            range.subject.length,
            redo_align_params.compo_adjust_mode,
            composition_workspace,
            &mut shared_q,
        ) {
            inputs.push(input);
        }
    }
    crate::core::composition_adjustment::x_newton_lanes::prefetch(inputs);
}

#[cfg(test)]
#[path = "kappa_tests.rs"]
mod tests;
