//! The subject loop of a TBLASTX query batch: the worker state of each thread, the writer thread
//! that collects the subjects' hits (native builds with `parallel`), and the loop that runs
//! `process_subject` for every subject of the batch, on the search pool's threads or on the calling
//! thread.
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607
//! (`BLAST_PreliminarySearchEngine`, whose subject loop is at :1411-1584), which every thread of
//! c++/src/algo/blast/api/prelim_stage.cpp:141-203
//! (`CBlastPrelimSearch::x_LaunchMultiThreadedSearch`) runs, at :145-188.
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_subject_next`,
//! `x_batch_pipe`). Its functions are statements of `search_query_batch` moved here in their
//! order; `process_subject` stays a closure of `search_query_batch`, which hands it to
//! `run_subject_loop`. What the loop reads is a struct of references or `Copy` values.

use super::{
    x_batch_pipe, x_subject_next, DiagStruct, OffsetPair, Subjects, TblastxHsp, TbxScore,
    WorkerState,
};
use crate::algorithm::tblastx::x_work_board::XBoard;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use crate::utils::threading::SearchPool;
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use rayon::iter::{IntoParallelIterator, ParallelIterator};
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
use {
    super::TblastxArgs, crate::algorithm::tblastx::hit_collector::HitCollector,
    crate::blastinput::fasta_reader::FastaRecord, std::sync::Arc,
};

// NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
// ```c
// while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//        != BLAST_SEQSRC_EOF) {
//    ...
//    status = s_BlastSearchEngineCore(...);
// }
// ```
/// Runs `body` for the subjects `0..subjects` in their order on this thread, with one worker state
/// made by `init`.
fn for_each_subjects<FInit, FBody, S>(subjects: usize, init: FInit, mut body: FBody)
where
    FInit: FnOnce() -> WorkerState<S>,
    FBody: FnMut(&mut WorkerState<S>, usize),
{
    let mut state = init();
    for s_idx in 0..subjects {
        body(&mut state, s_idx);
    }
}

impl<S> WorkerState<S> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
    // ```c
    //     /* iterate over all subject sequences */
    //     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //            != BLAST_SEQSRC_EOF) {
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
    // ```c
    //     aux_struct->offset_pairs =
    //       (BlastOffsetPair*) malloc(offset_array_size * sizeof(BlastOffsetPair));
    // ```
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:60-63
    // ```c
    //                 diag_table->diag_array_length = diag_array_length;
    //                 diag_table->diag_mask = diag_array_length-1;
    //                 diag_table->multiple_hits = multiple_hits;
    //                 diag_table->offset = window_size;
    // ```
    // Each thread of NCBI's search sets up its own offset pairs and diagonal table
    // (`s_BlastSetUpAuxStructures`) before its iterations of the subject loop.
    /// The state of one thread's iterations of the subject loop: the channel to the writer thread
    /// (`tx`, `None` when the loop keeps the collector), no hits yet, `offset_pairs_len` offset
    /// pairs, a diagonal table of `diag_cells` cells with its offset at `window`, and the scores
    /// of the search's matrix.
    fn new(
        tx: Option<std::sync::mpsc::Sender<(usize, Vec<TblastxHsp>)>>,
        offset_pairs_len: usize,
        diag_cells: usize,
        window: i32,
        sc: S,
    ) -> Self {
        WorkerState {
            tx,
            hits: Vec::new(),
            offset_pairs: vec![OffsetPair::default(); offset_pairs_len],
            diag_array: vec![DiagStruct::default(); diag_cells],
            diag_offset: window,
            sc,
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
// ```c
//          /* Save the results. */
//          status = BlastHSPStreamWrite(hsp_stream, &hsp_list);
// ```
/// Starts the writer thread of a parallel subject loop: it takes each subject's message from `rx`
/// until the channel closes and returns the batch's hit lists in report order.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
pub(super) fn spawn_writer(
    rx: std::sync::mpsc::Receiver<(usize, Vec<TblastxHsp>)>,
    args: &TblastxArgs,
    query_records: &[FastaRecord],
) -> std::thread::JoinHandle<Vec<TblastxHsp>> {
    // NCBI reference: ncbi-blast/c++/include/algo/blast/core/blast_hits.h:153-166
    // ```c
    // typedef struct BlastHSPList {
    //    Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
    //    Int4 query_index; /**< Index of the query which this HSPList corresponds to.
    //                       Set to 0 if not applicable */
    // } BlastHSPList;
    // ```
    // The writer thread is NCBI's collector (tblastx/hit_collector.rs): it takes each
    // subject's message and returns the batch's hit lists in report order after the join,
    // for the one place that writes every output.
    let mut collector = HitCollector::new(args, query_records);
    std::thread::spawn(move || -> Vec<TblastxHsp> {
        while let Ok((s_idx, hits)) = rx.recv() {
            collector.push(s_idx, hits);
        }
        collector.finish()
    })
}

/// What `run_subject_loop` reads.
pub(super) struct SubjectLoop<'a, S> {
    pub(super) use_parallel: bool,
    pub(super) use_serial_scan_chunks: bool,
    pub(super) subjects: Subjects<'a>,
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    pub(super) parallel_pool: &'a SearchPool<'a>,
    pub(super) offset_pairs_len: usize,
    pub(super) x_worker_diag_cells: usize,
    pub(super) window: i32,
    pub(super) sc: S,
    /// The sender of the channel to the writer thread, `None` without one.
    pub(super) tx_opt: &'a Option<std::sync::mpsc::Sender<(usize, Vec<TblastxHsp>)>>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): the search's board and the batch on it.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    pub(super) x_board_batch: &'a Option<(
        &'a XBoard,
        Arc<crate::algorithm::tblastx::x_work_board::XBatch>,
    )>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): whether the batch puts units on the board.
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    pub(super) x_board_unit_path: bool,
    /// EXPERIMENT (LOSAT_X_TBXPIPE): the batch pipe's view of this search's seed stage.
    pub(super) x_seed_end: Option<&'a x_batch_pipe::XSeedEnd<'a>>,
    /// EXPERIMENT (LOSAT_X_TBXBOARD): the search's board.
    pub(super) x_board: Option<&'a XBoard>,
}

/// Runs `process_subject` for every subject of the batch, each thread with its own worker state:
/// on the search pool's threads when the search is parallel and has several subjects (TBXNEXT's
/// loop, or rayon's `for_each_init` without it), otherwise on this thread in subject order. A
/// batch started early by the batch pipe waits for the batch before it first (`wait_gate`).
/// `process_subject` is `Sync + Send` because the threads of the parallel loops share it.
///
/// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1584 (the subject loop
/// of `BLAST_PreliminarySearchEngine`) and c++/src/algo/blast/api/prelim_stage.cpp:145-188
/// (`CBlastPrelimSearch::x_LaunchMultiThreadedSearch`, which runs it on every thread).
pub(super) fn run_subject_loop<S: TbxScore, F>(lp: &SubjectLoop<'_, S>, process_subject: F)
where
    F: Fn(&mut WorkerState<S>, usize) + Sync + Send,
{
    let &SubjectLoop {
        use_parallel,
        use_serial_scan_chunks,
        subjects,
        #[cfg(all(
            feature = "parallel",
            any(not(target_arch = "wasm32"), feature = "wasm-threads")
        ))]
        parallel_pool,
        offset_pairs_len,
        x_worker_diag_cells,
        window,
        sc,
        tx_opt,
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        x_board_batch,
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        x_board_unit_path,
        x_seed_end,
        x_board,
    } = lp;
    // NCBI reference: ncbi-blast/c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // ```c
    // if (num_threads > 1) {
    //     SetNumberOfThreads(num_threads);
    // }
    // ```
    // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1411-1475
    // ```c
    // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
    //        != BLAST_SEQSRC_EOF) {
    //    ...
    //    status = s_BlastSearchEngineCore(...);
    // }
    // ```
    // LOSAT_X_TBXNEXT: its loop takes the subjects of a parallel search of several subjects.
    if x_subject_next::on() {
        crate::utils::xreach::x_reach!(
            "LOSAT_X_TBXNEXT",
            use_parallel && !use_serial_scan_chunks && subjects.count() > 1
        );
    }

    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 {
        // NCBI reference: ncbi-blast/c++/src/algo/blast/core/blast_engine.c:1409-1427
        // ```c
        // itr = BlastSeqSrcIteratorNewEx(MAX(BlastSeqSrcGetNumSeqs(seq_src)/100,1));
        // /* iterate over all subject sequences */
        // while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
        //        != BLAST_SEQSRC_EOF) {
        //    if (BlastSeqSrcGetSequence(seq_src, &seq_arg) < 0) {
        //        continue;
        //    }
        // ```
        // The threads' subjects go to the collector the subject loop keeps (`emit_subject_hits`;
        // no writer thread on this target), which takes them in OID order.
        if x_subject_next::on() {
            // EXPERIMENT (LOSAT_X_TBXNEXT): as in native builds below, each thread of the pool
            // takes the next subject when it is free, with one worker state per thread made when
            // it takes its first subject (tblastx/x_subject_next.rs). Same `process_subject`.
            parallel_pool.install(|| {
                x_subject_next::for_each_subject(
                    subjects.count(),
                    &|| WorkerState::new(None, offset_pairs_len, x_worker_diag_cells, window, sc),
                    &process_subject,
                    None,
                );
            });
        } else {
            parallel_pool.install(|| {
                (0..subjects.count()).into_par_iter().for_each_init(
                    || WorkerState::new(None, offset_pairs_len, x_worker_diag_cells, window, sc),
                    process_subject,
                );
            });
        }
    } else {
        for_each_subjects(
            subjects.count(),
            || {
                WorkerState::new(
                    tx_opt.clone(),
                    offset_pairs_len,
                    x_worker_diag_cells,
                    window,
                    sc,
                )
            },
            &process_subject,
        );
    }

    // EXPERIMENT (LOSAT_X_TBXBOARD): a thread of the subject loop without a subject helps with the
    // batch's units (tblastx/x_work_board.rs `XBoard::help`).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    let x_board_help = x_board_batch.as_ref().map(|(board, batch)| {
        move |stopped: &std::sync::atomic::AtomicBool| {
            board.help(batch, stopped);
            // Every unit of the batch has been taken: the next batch may build its tables. A batch
            // without the unit path has no units; its next batch starts at its seed-stage end.
            if let Some(seed_end) = x_seed_end.filter(|_| x_board_unit_path) {
                seed_end.units_taken();
            }
        }
    });
    // EXPERIMENT (LOSAT_X_TBXBOARD): a batch started early waits here, after its tables, until the
    // batch before has ended its seed stage (tblastx/x_batch_pipe.rs `XGateView`).
    if let Some(seed_end) = x_seed_end {
        seed_end.wait_gate(x_board);
    }
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 && x_subject_next::on() {
        // EXPERIMENT (LOSAT_X_TBXNEXT): each thread of the pool takes the next subject when it is
        // free, with one worker state per thread, made by the same initializer as below when the
        // thread takes its first subject (tblastx/x_subject_next.rs). Same `process_subject`.
        parallel_pool.install(|| {
            x_subject_next::for_each_subject(
                subjects.count(),
                &|| {
                    WorkerState::new(
                        tx_opt.clone(),
                        offset_pairs_len,
                        x_worker_diag_cells,
                        window,
                        sc,
                    )
                },
                &process_subject,
                x_board_help
                    .as_ref()
                    .map(|help| help as &(dyn Fn(&std::sync::atomic::AtomicBool) + Sync)),
            );
        });
    } else if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 {
        parallel_pool.install(|| {
            (0..subjects.count()).into_par_iter().for_each_init(
                || {
                    WorkerState::new(
                        tx_opt.clone(),
                        offset_pairs_len,
                        x_worker_diag_cells,
                        window,
                        sc,
                    )
                },
                process_subject,
            );
        });
    } else {
        for_each_subjects(
            subjects.count(),
            || {
                WorkerState::new(
                    tx_opt.clone(),
                    offset_pairs_len,
                    x_worker_diag_cells,
                    window,
                    sc,
                )
            },
            &process_subject,
        );
    }

    #[cfg(any(
        not(feature = "parallel"),
        all(target_arch = "wasm32", not(feature = "wasm-threads"))
    ))]
    {
        for_each_subjects(
            subjects.count(),
            || {
                WorkerState::new(
                    tx_opt.clone(),
                    offset_pairs_len,
                    x_worker_diag_cells,
                    window,
                    sc,
                )
            },
            &process_subject,
        );
    }
}
