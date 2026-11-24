//! The subject loop of a TBLASTX query batch: the worker state of each thread, the writer thread
//! that collects the subjects' hits (native builds with `parallel`), the loop that runs
//! `process_subject` for every subject of the batch, on the search pool's threads or on the calling
//! thread, the fetch at the start of each subject's re-evaluation, which the negative seqid list
//! may exclude (`SubjectEmit::reevaluation_fetch_keeps`), and the end of each subject
//! (`SubjectEmit::emit`: its final HSPs to the writer thread or to the collector the loop keeps).
//!
//! NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1334-1607
//! (`BLAST_PreliminarySearchEngine`, whose subject loop is at :1411-1584), which every thread of
//! c++/src/algo/blast/api/prelim_stage.cpp:141-203
//! (`CBlastPrelimSearch::x_LaunchMultiThreadedSearch`) runs, at :145-188;
//! c++/src/algo/blast/core/blast_hits.c:2651-2665 (the fetch that a negative seqid list excludes)
//! and c++/src/algo/blast/core/blast_engine.c:1553-1554 (`BlastHSPStreamWrite`).
//!
//! This module is a child of `run_impl` for the modules it keeps private (`x_subject_next`,
//! `x_batch_pipe`). Its functions are statements of `search_query_batch` moved here in their
//! order; `process_subject` stays a closure of `search_query_batch`, which hands it to
//! `run_subject_loop`. What the loop reads is a struct of references or `Copy` values, as is
//! `SubjectEmit`, whose `emit` (the closure `emit_subject_hits` before) `process_subject` calls at
//! the end of a subject.

use super::{
    x_batch_pipe, x_subject_next, DiagStruct, OffsetPair, Subjects, TblastxHsp, TbxScore,
    WorkerState,
};
use crate::algorithm::common::database::FetchExcluded;
use crate::algorithm::tblastx::hit_collector::HitCollector;
use crate::algorithm::tblastx::subjects::XReadFailure;
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
use {super::TblastxArgs, crate::blastinput::fasta_reader::FastaRecord, std::sync::Arc};

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1475
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1040-1041
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_extend.c:60-63
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
/// Starts the writer thread of a parallel subject loop: it takes each subject's message from `rx`
/// until the channel closes and returns the batch's hit lists in report order.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
pub(super) fn spawn_writer(
    rx: std::sync::mpsc::Receiver<(usize, Vec<TblastxHsp>)>,
    args: &TblastxArgs,
    query_records: &[FastaRecord],
    sorts_results: bool,
) -> std::thread::JoinHandle<Vec<TblastxHsp>> {
    // NCBI reference (598d8ae6): c++/include/algo/blast/core/blast_hits.h:153-166
    // The writer thread is NCBI's collector (tblastx/hit_collector.rs): it takes each
    // subject's message and returns the batch's hit lists in report order after the join,
    // for the one place that writes every output.
    let mut collector = HitCollector::new(args, query_records, sorts_results);
    std::thread::spawn(move || -> Vec<TblastxHsp> {
        while let Ok((s_idx, hits)) = rx.recv() {
            collector.push(s_idx, hits);
        }
        collector.finish()
    })
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2651-2665
// The re-evaluation's fetch of a subject with HSPs excludes one that the negative seqid
// list keeps without a line (`FetchExcluded`): its HSPs are freed before anything else of the
// re-evaluation reads it (`SubjectEmit::reevaluation_fetch_keeps`), one check for the batch.
/// What the re-evaluation's fetch (`reevaluation_fetch_keeps`) and the end of each subject (`emit`)
/// read of the batch: the negative seqid list's check of that fetch (`fetch_excluded`, `None` where
/// the fetches exclude nothing), the batch's read failure, and the collector of the hit lists that
/// a subject loop without a writer thread keeps (`local_collector`).
pub(super) struct SubjectEmit<'a> {
    pub(super) fetch_excluded: Option<&'a FetchExcluded<'a>>,
    pub(super) x_read_failure: &'a XReadFailure,
    pub(super) local_collector: &'a std::sync::Mutex<Option<HitCollector>>,
}

impl SubjectEmit<'_> {
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1553-1554
    // NCBI writes a subject's list when it holds HSPs (:1480); the collector takes one message per
    // subject, with none too, to keep the subjects' order.
    /// One message per subject of the loop, with its final HSPs or none, so that the collector
    /// reaches every OID (`HitCollector::push`): `st.hits` holds the subject's HSPs until here.
    pub(super) fn emit<S>(&self, st: &mut WorkerState<S>, s_idx: usize) {
        let hits = std::mem::take(&mut st.hits);
        if let Some(tx) = &st.tx {
            tx.send((s_idx, hits)).unwrap();
        } else {
            self.local_collector
                .lock()
                .unwrap()
                .as_mut()
                .expect("a subject loop without a writer thread keeps the collector")
                .push(s_idx, hits);
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2659-2665
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/seqsrc_seqdb.cpp:306-311
    // NCBI's re-evaluation (`Blast_HSPListReevaluateUngapped`, called for every subject that has
    // HSPs after the first linking call and the preliminary E-value, blast_engine.c:1480-1497)
    // fetches the subject first, reading its Seq-ids before its letters: one that the negative
    // seqid list excludes loses its HSPs there, before its letters, its translation or anything
    // else of the re-evaluation is read; one whose header cannot be read ends NCBI's search.
    /// Whether the re-evaluation's fetch of subject `s_idx`, which still has HSPs, keeps it: not
    /// when the negative seqid list excludes it (`FetchExcluded`), nor when its header cannot be
    /// read, which fails the batch with LOSAT's refusal after the reports of the batches before
    /// (`x_fetch_keeps`; `XReadFailure`, as a subject whose read failed fails it).
    #[inline]
    pub(super) fn reevaluation_fetch_keeps(&self, s_idx: usize) -> bool {
        match (self.fetch_excluded, u32::try_from(s_idx)) {
            (Some(excluded), Ok(ordinal)) => x_fetch_keeps(excluded, self.x_read_failure, ordinal),
            _ => true,
        }
    }

    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_hits.c:2672-2677,2700
    // After the fetch (`reevaluation_fetch_keeps`, then the subject's letters), NCBI's
    // re-evaluation reads the translation of the first HSP's frame, which it never makes for a
    // subject (its `-subject_loc` range) of fewer than 12 bases
    // (`core::blast_hits::target_translation_never_made`): NULL + 1, from which
    // `Blast_HSPReevaluateWithAmbiguitiesUngapped` reads the HSP's letters (blast_hits.c:692-700),
    // and NCBI BLAST+ 2.17.0 dies with SIGSEGV. Every HSP that reaches LOSAT's re-evaluation has a
    // length of 1 or more (those of length 0 or less, which read nothing in NCBI's, are dropped
    // after the first linking call, `subject_link::link_subject_hits_int4`, or without sum
    // statistics before this refusal, `subject_link::drop_unread_without_linking`), so NCBI dies
    // at the first HSP of every such list. The query batch fails with LOSAT's refusal
    // (PD-LOSAT-NCBI-DEFECTS 1.51, rule 3) after the reports of the batches before
    // (`ErrorAfterBatchReports`), as the earliest such subject's error (`XReadFailure`).
    /// The refusal of subject `s_idx`, whose HSPs reach the re-evaluation's read of the translation
    /// of a subject of fewer than 12 bases. The caller clears the subject's list and goes on to the
    /// subject's end, which emits it without HSPs: one call out of line and no return of its own,
    /// so the subject's path holds no second copy of its ending.
    #[cold]
    #[inline(never)]
    pub(super) fn refuse_untranslated_frame(&self, s_idx: usize) {
        self.x_read_failure.fail(
            s_idx,
            anyhow::Error::new(super::ErrorAfterBatchReports::new(
                crate::core::blast_hits::x_untranslated_frame_refusal(
                    "TBLASTX",
                    "the ungapped re-evaluation",
                ),
            )),
        );
    }
}

// NCBI reference (598d8ae6): c++/src/app/blast/tblastx_app.cpp:203-209
// A header that NCBI's fetch cannot read throws in `GetSeqIDs` (a `CSerialException`), which
// NCBI reports as a `BLAST engine error` (exit 3) after it has printed the query batches before;
// LOSAT's refusal comes there too (`ErrorAfterBatchReports`; the Owner's decision of 2026-11-16
// on line RG's Q-X2c: "the results of the query batches before are written"). The run passes on
// the refusal itself (`FetchRefusal`).
/// The check of `reevaluation_fetch_keeps` under a negative seqid list: whether the fetch
/// of subject `ordinal` keeps it; a header that cannot be read is the batch's failure, after the
/// reports of the batches before.
#[inline(never)]
fn x_fetch_keeps(
    excluded: &FetchExcluded<'_>,
    x_read_failure: &XReadFailure,
    ordinal: u32,
) -> bool {
    match excluded.excludes(ordinal) {
        Ok(excluded) => !excluded,
        Err(error) => {
            x_read_failure.fail(
                ordinal as usize,
                anyhow::Error::new(super::ErrorAfterBatchReports::new(FetchExcluded::refusal(
                    (ordinal, error),
                    "TBLASTX",
                ))),
            );
            false
        }
    }
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
    // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:82-88
    // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1475
    // LOSAT_X_TBXNEXT: its loop takes the subjects of a parallel search of several subjects.
    if x_subject_next::on() {
        crate::utils::xreach::x_reach!(
            "LOSAT_X_TBXNEXT",
            use_parallel && !use_serial_scan_chunks && subjects.count() > 1
        );
    }

    #[cfg(all(feature = "parallel", target_arch = "wasm32", feature = "wasm-threads"))]
    if use_parallel && !use_serial_scan_chunks && subjects.count() > 1 {
        // NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1409-1427
        // The threads' subjects go to the collector the subject loop keeps (`SubjectEmit::emit`;
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
