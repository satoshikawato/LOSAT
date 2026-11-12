//! The end of a TBLASTX query batch's search, after its subject loop: the batch's hit lists in
//! report order (`collect_hits`) and the batch's `LOSAT_TIMING` lines (`timing_lines`).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1549-1649
//! (`BlastHitList2SeqAlign_OMF`), which takes the hit lists in their order. The timing lines
//! have no NCBI counterpart.
//!
//! This module is a child of `run_impl` for the timer it keeps private (`x_search_timing`). Its
//! functions are statements of `search_query_batch` moved here in their order; what they read is
//! the writer thread and the collector, which they take over, or a struct of references and
//! `Copy` values.

use super::{x_search_timing, AtomicOrdering, AtomicU64, Instant, TblastxHsp};
use crate::algorithm::tblastx::hit_collector::HitCollector;
use crate::algorithm::tblastx::x_link_timing::XLinkTimer;
use std::time::Duration;

/// The batch's hit lists in report order: from the writer thread (`writer`, joined here once the
/// subject loop's senders are dropped), or from the collector the subject loop fed.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577 (in
/// `BlastHitList2SeqAlign_OMF`).
pub(super) fn collect_hits(
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))] writer: Option<
        std::thread::JoinHandle<Vec<TblastxHsp>>,
    >,
    local_collector: std::sync::Mutex<Option<HitCollector>>,
) -> Vec<TblastxHsp> {
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/blast_seqalign.cpp:1569-1577
    // ```c
    // for (int index = 0; index < hit_list->hsplist_count; index++) {
    //     BlastHSPList* hsp_list = hit_list->hsplist_array[index];
    //     if (!hsp_list)
    //         continue;
    //     Blast_HSPListSortByEvalue(hsp_list);
    // }
    // ```
    // Each traversal ends with the batch's hit lists in report order (tblastx/hit_collector.rs):
    // from the writer thread, or from the collector the subject loop fed.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let collected_hits: Option<Vec<TblastxHsp>> = writer.map(|writer| writer.join().unwrap());
    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    let collected_hits: Option<Vec<TblastxHsp>> = None;
    match collected_hits {
        Some(hits) => hits,
        None => local_collector
            .into_inner()
            .unwrap()
            .expect("a subject loop without a writer thread keeps the collector")
            .finish(),
    }
}

/// What `timing_lines` reads: the flag, the clocks and the counters of the batch's search (the
/// atomics are those `process_subject` adds to).
pub(super) struct TimingInputs<'a> {
    pub(super) timing_enabled: bool,
    pub(super) t_search_start: Instant,
    pub(super) scan_ns: &'a AtomicU64,
    pub(super) scan_calls: &'a AtomicU64,
    pub(super) ungapped_ns: &'a AtomicU64,
    pub(super) ungapped_calls: &'a AtomicU64,
    pub(super) reeval_ns: &'a AtomicU64,
    pub(super) reeval_calls: &'a AtomicU64,
    pub(super) linking_ns: &'a AtomicU64,
    pub(super) linking_calls: &'a AtomicU64,
    pub(super) t_read_queries: Duration,
    pub(super) t_build_lookup: Duration,
    pub(super) t_read_subjects: Duration,
    pub(super) t_total: Instant,
    pub(super) x_times: &'a x_search_timing::XSearchTimes,
    pub(super) x_t_subjects: Option<Instant>,
    pub(super) x_t_collect: Option<Instant>,
    pub(super) x_link_timer: &'a XLinkTimer,
}

// No NCBI counterpart: the lines of LOSAT_TIMING read clocks and counters that the search does
// not read.
/// The batch's `LOSAT_TIMING` lines: the `[TIMING]` lines of the phases and stages, then the stage
/// lines of x_search_timing.rs and the linking lines of x_link_timing.rs; none without
/// `timing_enabled`. With it, the search's interval is also added to the common line
/// (utils/xtiming.rs).
pub(super) fn timing_lines(inputs: &TimingInputs<'_>) -> Vec<String> {
    let &TimingInputs {
        timing_enabled,
        t_search_start,
        scan_ns,
        scan_calls,
        ungapped_ns,
        ungapped_calls,
        reeval_ns,
        reeval_calls,
        linking_ns,
        linking_calls,
        t_read_queries,
        t_build_lookup,
        t_read_subjects,
        t_total,
        x_times,
        x_t_subjects,
        x_t_collect,
        x_link_timer,
    } = inputs;
    let mut x_timing: Vec<String> = Vec::new();
    if timing_enabled {
        let t_search = t_search_start.elapsed();
        // Diagnostic (LOSAT_TIMING): the same interval for the common line (utils/xtiming.rs).
        // No NCBI counterpart.
        crate::utils::xtiming::add_duration(crate::utils::xtiming::Phase::Search, t_search);
        let scan_s = scan_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let scan_n = scan_calls.load(AtomicOrdering::Relaxed);
        let ungapped_s = ungapped_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let ungapped_n = ungapped_calls.load(AtomicOrdering::Relaxed);
        let reeval_s = reeval_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let reeval_n = reeval_calls.load(AtomicOrdering::Relaxed);
        let linking_s = linking_ns.load(AtomicOrdering::Relaxed) as f64 / 1e9;
        let linking_n = linking_calls.load(AtomicOrdering::Relaxed);

        x_timing.push(format!(
            "[TIMING] read_queries: {:.3}s",
            t_read_queries.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] build_lookup: {:.3}s",
            t_build_lookup.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] read_subjects: {:.3}s",
            t_read_subjects.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] scan_subject: {:.3}s (calls={})",
            scan_s, scan_n
        ));
        x_timing.push(format!(
            "[TIMING] ungapped_extend: {:.3}s (calls={})",
            ungapped_s, ungapped_n
        ));
        x_timing.push(format!(
            "[TIMING] reevaluate: {:.3}s (calls={})",
            reeval_s, reeval_n
        ));
        x_timing.push(format!(
            "[TIMING] sum_stats_linking: {:.3}s (calls={})",
            linking_s, linking_n
        ));
        x_timing.push(format!(
            "[TIMING] search_total: {:.3}s",
            t_search.as_secs_f64()
        ));
        x_timing.push(format!(
            "[TIMING] total: {:.3}s",
            t_total.elapsed().as_secs_f64()
        ));
        // No NCBI counterpart: the stage lines of x_search_timing.rs, after the lines above.
        if let Some(lines) = x_times.lines(t_search_start, x_t_subjects, x_t_collect, t_search) {
            x_timing.extend(lines);
        }
        // No NCBI counterpart: the linking lines of x_link_timing.rs, after the stage lines.
        if let Some(lines) = x_link_timer.lines() {
            x_timing.extend(lines);
        }
    }
    x_timing
}
