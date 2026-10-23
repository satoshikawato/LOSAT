//! PROTOTYPE (v4-A s6, blastx memory at many threads; switches read as plain variables, opt-in):
//!
//! - `LOSAT_X_BXDIAGTLS=1`: the diagonal tables of the pool's workers (subject waves, the subject
//!   flow's helpers, the chunks' subject tasks below) are one table per thread, kept by the thread
//!   between searches, instead of `threads` new tables made by the searching thread for each query
//!   batch or chunk searched with the pool inside.
//! - `LOSAT_X_BXUNITCAP=K`: at more than K threads, at most K query batches (`LOSAT_X_BXFLOW`) or
//!   chunks of a split query are searched at a time, each with the pool inside (its subjects as
//!   tasks), instead of one batch or chunk per thread.
//!
//! A table's history does not change what a subject computes while no diagonal offset passes
//! 2^30 (`x_subject_flow` module comment); both apply only where `x_subject_flow` allows that.

use super::{compute_subject, SubjectComputation};
use crate::algorithm::blastx::{
    args::ResolvedOptions,
    parameters::ContextParameters,
    query_setup::PreparedQueryBatch,
    seed::{Diagonals, Lookup},
    subjects::Subjects,
};
use anyhow::Result;
use std::cell::RefCell;
use std::sync::OnceLock;

// No NCBI counterpart: reads an experiment switch.
fn env_on(name: &str) -> bool {
    std::env::var_os(name).is_some_and(|v| !v.is_empty() && v != "0")
}

/// PROTOTYPE (LOSAT_X_BXDIAGTLS): whether the workers keep one diagonal table per thread.
pub(crate) fn diag_tls(options: &ResolvedOptions, subjects: Subjects<'_>) -> bool {
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| env_on("LOSAT_X_BXDIAGTLS"))
        && super::x_subject_flow::offsets_stay_below_2_30(
            options.window_size,
            subjects.longest_length(),
        )
}

/// PROTOTYPE (LOSAT_X_BXUNITCAP=K): `Some(K)` where more than K threads search.
pub(crate) fn unit_cap(threads: usize) -> Option<usize> {
    static K: OnceLock<Option<usize>> = OnceLock::new();
    K.get_or_init(|| {
        std::env::var("LOSAT_X_BXUNITCAP")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|k| (2..=4096).contains(k))
    })
    .filter(|&k| threads > k)
}

/// PROTOTYPE (LOSAT_X_BXUNITCAP): subjects per task of a chunk searched with the pool inside.
fn subjects_per_task() -> usize {
    static K: OnceLock<usize> = OnceLock::new();
    *K.get_or_init(|| {
        std::env::var("LOSAT_X_BXUNITCAPTASK")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
            .filter(|k| (1..=65_536).contains(k))
            .unwrap_or(16)
    })
}

thread_local! {
    // No NCBI counterpart: the diagonal table this thread keeps between searches (buffer reuse).
    static X_DIAG: RefCell<Option<Diagonals>> = const { RefCell::new(None) };
}

/// No NCBI counterpart (buffer reuse): a worker's diagonal table, taken on first use from the
/// thread that uses it (its kept table where the length and window fit, else a new one) and given
/// back to the thread that drops it.
pub(crate) struct XDiagSlot {
    table: Option<Diagonals>,
    length: usize,
    window: i32,
}

impl XDiagSlot {
    pub(crate) fn new(length: usize, window: i32) -> Self {
        Self {
            table: None,
            length,
            window,
        }
    }

    pub(crate) fn get(&mut self) -> Result<&mut Diagonals> {
        if self.table.is_none() {
            let kept = X_DIAG.with(|kept| kept.borrow_mut().take());
            self.table = Some(match kept {
                Some(table) if table.x_fits(self.length, self.window) => table,
                _ => Diagonals::new(self.length, self.window)?,
            });
        }
        Ok(self.table.as_mut().expect("set above"))
    }
}

impl Drop for XDiagSlot {
    fn drop(&mut self) {
        if let Some(table) = self.table.take() {
            let _ = X_DIAG.try_with(|kept| *kept.borrow_mut() = Some(table));
        }
    }
}

/// No NCBI counterpart (buffer reuse): the diagonal table a subject-flow task runs on: the
/// owner's, a helper's own, or a helper's thread's (`XDiagSlot`).
pub(crate) enum XScratch<'a> {
    Borrowed(&'a mut Diagonals),
    Owned(Diagonals),
    Slot(XDiagSlot),
}

impl XScratch<'_> {
    pub(crate) fn get(&mut self) -> Result<&mut Diagonals> {
        match self {
            Self::Borrowed(table) => Ok(table),
            Self::Owned(table) => Ok(table),
            Self::Slot(slot) => slot.get(),
        }
    }
}

// NCBI reference (598d8ae6): c++/src/algo/blast/core/blast_engine.c:1411-1413
// ```c
//     /* iterate over all subject sequences */
//     while ( (seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr))
//            != BLAST_SEQSRC_EOF) {
// ```
/// PROTOTYPE (LOSAT_X_BXUNITCAP): `x_compute_chunk`'s subject loop as tasks of
/// `subjects_per_task()` consecutive subjects on the pool, each on its thread's diagonal table; the
/// subjects with something for the sink, in subject order, or the error of the first subject that
/// failed.
#[allow(clippy::too_many_arguments)]
pub(super) fn chunk_subjects_par(
    subjects: Subjects<'_>,
    batch: &PreparedQueryBatch,
    params: &[ContextParameters],
    options: &ResolvedOptions,
    lookup: &Lookup,
    xdrop: i32,
    db_length: i64,
    x_containment: bool,
) -> Result<Vec<(usize, SubjectComputation)>> {
    use rayon::prelude::*;
    let last = batch.contexts.last().expect("BLASTX contexts");
    let length = last.offset + last.length;
    let n = subjects.len();
    let k = subjects_per_task();
    let parts: Vec<Result<Vec<(usize, SubjectComputation)>>> = (0..n.div_ceil(k))
        .into_par_iter()
        .with_max_len(1)
        .map(|task| {
            let (first, end) = (task * k, n.min((task + 1) * k));
            let block = subjects.x_read_range(first, end)?;
            let subjects = subjects.x_with_block(block.as_ref());
            let mut slot = XDiagSlot::new(length, options.window_size);
            let table = slot.get()?;
            let mut kept = Vec::new();
            for oid in first..end {
                let mut computed = compute_subject(
                    subjects,
                    batch,
                    params,
                    options,
                    lookup,
                    table,
                    xdrop,
                    false,
                    db_length,
                    oid,
                    x_containment,
                )?;
                computed.seeds = Vec::new();
                computed.initial = Vec::new();
                if !computed.is_empty() {
                    kept.push((oid, computed));
                }
            }
            Ok(kept)
        })
        .collect();
    let mut results = Vec::new();
    for part in parts {
        results.extend(part?);
    }
    Ok(results)
}
