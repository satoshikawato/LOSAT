//! EXPERIMENT (LOSAT_X_BPQSLOTS): blastp-fast's query blocks (`LOSAT_X_BPQBLOCK`, `x_query_blocks.rs`)
//! searched several at a time on the search's pool, their reports written in block order.
//!
//! The block loop searches one block at a time on the whole pool: the block's lookup table, its scan
//! of the subjects, its gapped stage, its redo and its report, each stage a parallel loop that ends
//! in a join, and the serial parts between them (the report, the set-up). In the human x yeast
//! blastp-fast search at 16 threads (Species4 x Species6 `-subject`, 254 blocks of about 18 ms) about
//! 1.8 of the 16 threads waited all the time, and a drain at most block ends left fewer than 8 busy
//! (line R `REPORT_R17.md` §3.1: 2.44 of the 3.9 points that S(16) lacked to 16).
//!
//! With the switch on, `k` drivers run on the search's pool (one `rayon::scope`); each takes the next
//! block not yet taken and searches it with the block loop's call (`run_resolved_in_pool` with the
//! block's queries, ranges and report flags, on the same pool), its report and messages written into
//! a buffer of its own. A thread that waits in one block's join takes the tasks of another block. A
//! block's buffers are written to the outputs when the blocks before it are written, by the driver that
//! finished last, so the outputs get the blocks' bytes in block order.
//!
//! Why the outputs are those of the block loop:
//! - a block's search reads the query records, ranges and report flags of its block (the first block
//!   writes the outfmt 0 prolog, the last the epilog, outfmt 7's count and outfmt 5's iteration numbers
//!   start at the block's first query: `BlastpReports`), the search's arguments, the subjects encoded
//!   once for all blocks and, with `-db`, the database (read only), and writes only its outputs and
//!   messages. It reads nothing that another block's search writes: the lookup table, the diagonal
//!   table, the hit lists and the reports are the block's own, and the spare vectors of the flat tables
//!   (`XFlatSpare`) are moved under their locks, their values never read (`x_clut_par.rs`). The flat
//!   table does not depend on the threads that build it (`x_clut_par.rs`), nor any value of the search
//!   on the thread count (the block loop gives the same bytes at every thread count);
//! - the run-level state that outlives a block's report is declined: outfmt 13 to 16 (`Xml2Carry`),
//!   outfmt 11 (the archive handle carries the search space and the messages from document to
//!   document, `x_query_blocks.rs`), several outputs, the hit records and an observer (the observer's
//!   indices count the HSPs of the blocks before). Those take the block loop;
//! - the buffers of the blocks are written one after the other in block order, the messages of a block
//!   before its report, and the report's output flushed after each block as the block loop flushes it
//!   at the end of each block's report (a file is created at block 0's flush, and a failed write fails at
//!   the block where the block loop fails), so the bytes on each stream are the block loop's. Where the
//!   report and the messages go to one stream (`2>&1`) and a block's report passes the output's 8 KiB
//!   buffer, the block loop's writes inside the block can come between that block's messages; here all of
//!   a block's messages come before its report.
//!
//! Errors: the first block in block order whose search fails ends the search. No block after it is
//! started once the failure is known; the blocks before it are written, then the failed block's
//! messages and the bytes its report had written, then its error is returned, as the block loop
//! returns it after the same bytes. A failed write to an output ends the search at that block. The
//! release profile has `panic = "abort"`: a panic in a block's search ends the process, where the
//! blocks written so far may be fewer than in the block loop (blocks after it may have been searched,
//! blocks before it not yet written). Where panics unwind (tests, debug builds) a panicking driver stops
//! the other (`XPanicStop`), and the scope resumes the panic.
//!
//! Waiting: a driver waits (a `Condvar`) only between blocks, before it takes one, never inside a block's
//! search. Both drivers are spawned before any block's work, and a rayon worker takes the oldest job of
//! another's deque first, so each driver runs at the bottom of an idle worker's stack, not inside another
//! block's join; nothing on the blastp-fast path waits for every worker (no `broadcast`, no barrier).
//! Spawning a driver later, or inside a block's scope, would break this.
//!
//! Memory: each block in flight holds its lookup table and its hit lists until its report is written,
//! and then gives the table's vectors to the block search's spare (`LOSAT_X_BPLKPFLAT`), which a later
//! block's build takes (`XAheadBlock::none`: a later block may take it, so `LOSAT_X_BPLKDROP` keeps it
//! through the redo). The block loop holds one block and, with `LOSAT_X_BPLKPAHEAD`, the next block's
//! table built ahead: two blocks' tables, `2 x LOSAT_X_BPQBLOCK` letters. The blocks of `k` drivers hold
//! `2 x LOSAT_X_BPQBLOCK / (k + 1)` letters each (`x_slot_letters`), so the `k` tables in flight and the
//! build of one more hold the letters of the block loop's two (the build's transient vectors are as
//! large as a table: 78 appends of 4 bytes per query letter recorded beside the 4-byte hits, line BPF
//! s7). No block is built ahead. A driver takes a block only while fewer than `2 k` blocks are searched
//! or wait to be written: at most `k` finished reports wait behind a slow block (the block loop holds
//! none past its flush, but holds that block's hit lists, of which a report is made).
//!
//! Where it acts: the block search with the compressed lookup table (`CompressedAaLookupTable`:
//! `-task blastp-fast`, or a word size above 4) on a pool of more than one thread, with one output that is
//! not outfmt 11 or 13 to 16, no hit records and no observer, in builds with threads (native, or wasm32
//! with `wasm-threads`). The table `AaLookupTable` (exact blastp) keeps the block loop (its in-place
//! builds share neighbour lists that only one build uses at a time, `x_lookup_ahead.rs`).
//! `LOSAT_TIMING`'s lines of the blocks in flight come interleaved (each block prints its own, and the
//! redo's range counters that each block's start resets are shared by the two blocks).
//!
//! One form, not tuned: two blocks in flight (`X_SLOTS`), the block loop's block and the next one, as
//! `LOSAT_X_BPLKPAHEAD` keeps the next block's table beside the block (here its whole search), each of
//! the letters that keep the tables in flight to the block loop's two. `LOSAT_X_BPQSLOTS` is a switch of
//! the strict set (`utils/xswitch.rs`), on unless it is `0` (the block loop).
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup()) {
//!     ...
//!                 results = lcl_blast.Run();
//!     ...
//!                     formatter.PrintOneResultSet(**result, query_batch);
//! ```
//! NCBI searches and reports one batch at a time; a block is whole batches, and the blocks' reports
//! are written in the order of the batches.

use std::sync::atomic::{AtomicU64, Ordering};
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
use {
    crate::api::local_blast::{FormatOutput, OutputSink, ReportOutputs},
    crate::utils::threading::SearchPool,
    anyhow::Result,
    std::collections::BTreeMap,
    std::io::Write,
    std::sync::atomic::AtomicBool,
    std::sync::{Condvar, Mutex, PoisonError},
};

/// `LOSAT_X_STATS`: the block searches that took the slots, and their blocks.
static X_SEARCHES: AtomicU64 = AtomicU64::new(0);
static X_BLOCKS: AtomicU64 = AtomicU64::new(0);

// No NCBI counterpart: diagnostics only.
/// `LOSAT_X_STATS`: `[X_STATS] BPQSLOTS searches=.. blocks=..` at exit, when a count is not zero.
pub(crate) fn x_bp_qslots_print_summary() {
    if std::env::var_os("LOSAT_X_STATS").is_some() {
        let (searches, blocks) = (
            X_SEARCHES.load(Ordering::Relaxed),
            X_BLOCKS.load(Ordering::Relaxed),
        );
        if searches != 0 || blocks != 0 {
            eprintln!("[X_STATS] BPQSLOTS searches={searches} blocks={blocks}");
        }
    }
}

/// EXPERIMENT (LOSAT_X_BPQSLOTS): the blocks in flight (module comment).
#[cfg_attr(
    not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )),
    allow(dead_code)
)]
pub(crate) const X_SLOTS: usize = 2;

/// PROBE (v4-BPF s8): a numeric environment variable.
pub(crate) fn x_s8_env(name: &str) -> Option<usize> {
    std::env::var(name).ok()?.parse().ok()
}

/// PROBE (v4-BPF s8): with LOSAT_BPF_S8LOG=1, one stderr line per event: microseconds since the
/// first event, the event, two numbers and the process's RssAnon (kB).
pub(crate) fn x_s8_log(event: &str, a: usize, b: usize) {
    use std::sync::OnceLock;
    static ON: OnceLock<bool> = OnceLock::new();
    static START: OnceLock<std::time::Instant> = OnceLock::new();
    if !*ON.get_or_init(|| std::env::var("LOSAT_BPF_S8LOG").is_ok_and(|v| v == "1")) {
        return;
    }
    let t = START.get_or_init(std::time::Instant::now).elapsed().as_micros();
    let anon = std::fs::read_to_string("/proc/self/status")
        .ok()
        .and_then(|status| {
            status
                .lines()
                .find(|line| line.starts_with("RssAnon:"))
                .and_then(|line| line.split_whitespace().nth(1)?.parse::<usize>().ok())
        })
        .unwrap_or(0);
    eprintln!("S8\t{t}\t{event}\t{a}\t{b}\t{anon}");
}

// No NCBI counterpart: reads the LOSAT_X_BPQSLOTS switch once; it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPQSLOTS): the number of blocks in flight (`X_SLOTS`), or `None` (off, or a
/// build without threads). A strict-set switch, on unless it is `0`.
pub(crate) fn x_bp_qslots() -> Option<usize> {
    #[cfg(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    ))]
    {
        use std::sync::OnceLock;
        static ON: OnceLock<bool> = OnceLock::new();
        ON.get_or_init(|| crate::utils::xswitch::on("LOSAT_X_BPQSLOTS"))
            .then(|| x_s8_env("LOSAT_BPF_S8K").unwrap_or(X_SLOTS))
    }
    #[cfg(not(all(
        feature = "parallel",
        any(not(target_arch = "wasm32"), feature = "wasm-threads")
    )))]
    None
}

// No NCBI counterpart: the size of LOSAT's blocks (module comment); it does not change any value NCBI
// computes.
/// EXPERIMENT (LOSAT_X_BPQSLOTS): the letters of a block (`x_query_blocks`'s `block_letters`) when
/// `slots` blocks are in flight: `2 x block_letters / (slots + 1)`, at least 1.
pub(crate) fn x_slot_letters(block_letters: usize, slots: usize) -> usize {
    if let Some(letters) = x_s8_env("LOSAT_BPF_S8L") {
        return letters;
    }
    (block_letters.saturating_mul(2) / (slots + 1)).max(1)
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
/// The blocks searched and not yet written, and the outputs they are written to in block order.
struct XOrdered<'s> {
    /// The next block to write.
    next: usize,
    /// The blocks searched, by block: the report's bytes, the messages and the search's error.
    done: BTreeMap<usize, XDone>,
    /// The blocks taken: those below are searched, being searched or written.
    taken: usize,
    sink: &'s mut (dyn Write + Send),
    diagnostics: &'s mut (dyn Write + Send),
    /// The first failure in block order: its block and error (a failed search or write).
    failed: Option<(usize, anyhow::Error)>,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
/// A block's search, kept until the blocks before it are written.
struct XDone {
    report: Vec<u8>,
    messages: Vec<u8>,
    error: Option<anyhow::Error>,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl XOrdered<'_> {
    /// Keeps `failure` at `index` when it comes before the failure kept.
    fn fail(&mut self, index: usize, failure: anyhow::Error) {
        if self.failed.as_ref().is_none_or(|&(at, _)| index < at) {
            self.failed = Some((index, failure));
        }
    }

    /// Writes the blocks searched from `next` on, in block order, up to the first that is not
    /// searched yet, or through the first that failed.
    fn write_ready(&mut self) {
        while self.failed.as_ref().is_none_or(|&(at, _)| self.next < at) {
            let Some(done) = self.done.remove(&self.next) else {
                return;
            };
            let at = self.next;
            let written = self
                .diagnostics
                .write_all(&done.messages)
                .and_then(|()| self.sink.write_all(&done.report))
                .and_then(|()| self.sink.flush());
            self.next += 1;
            if let Some(error) = done.error {
                self.fail(at, error);
            } else if let Err(error) = written {
                self.fail(at, error.into());
            }
        }
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
/// Stops the drivers when the one that holds it panics (module comment): `stop` set under the lock that
/// the waiting driver reads it under, then the waiter woken.
#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
struct XPanicStop<'a, 's> {
    stop: &'a AtomicBool,
    ordered: &'a Mutex<XOrdered<'s>>,
    written: &'a Condvar,
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
impl Drop for XPanicStop<'_, '_> {
    fn drop(&mut self) {
        if std::thread::panicking() {
            let state = self.ordered.lock().unwrap_or_else(PoisonError::into_inner);
            self.stop.store(true, Ordering::Relaxed);
            drop(state);
            self.written.notify_all();
        }
    }
}

#[cfg(all(
    feature = "parallel",
    any(not(target_arch = "wasm32"), feature = "wasm-threads")
))]
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:257-276 (module comment): the batches
// searched and reported in order. Scheduling only: the blocks' reports are written in that order.
/// EXPERIMENT (LOSAT_X_BPQSLOTS): the blocks `0..blocks` searched by `slots` drivers on `pool`
/// (`search_block(index, outputs)` writes block `index`'s report of `outfmt` and its messages to
/// `outputs`), the reports written to `sink` and the messages to `diagnostics` in block order; the
/// first failure in block order is returned after the bytes before it (module comment).
pub(crate) fn x_search_in_slots(
    blocks: usize,
    slots: usize,
    outfmt: &str,
    pool: &SearchPool<'_>,
    sink: &mut (dyn Write + Send),
    diagnostics: &mut (dyn Write + Send),
    search_block: &(dyn Fn(usize, &mut ReportOutputs<'_>) -> Result<()> + Sync),
) -> Result<()> {
    let window = slots.saturating_mul(2).max(1);
    let stop = AtomicBool::new(false);
    let ordered = Mutex::new(XOrdered {
        next: 0,
        done: BTreeMap::new(),
        taken: 0,
        sink,
        diagnostics,
        failed: None,
    });
    let written = Condvar::new();
    // PROBE (v4-BPF s8): every driver starts before any takes a block (the fix of fix/bp-slots-start).
    let drivers = slots.min(blocks).min(pool.threads()).max(1);
    let started = std::sync::atomic::AtomicUsize::new(0);
    let wait_for_every_driver = || {
        started.fetch_add(1, Ordering::AcqRel);
        let mut idle = 0u32;
        while started.load(Ordering::Acquire) < drivers {
            idle = idle.saturating_add(1);
            if idle < 128 {
                std::hint::spin_loop();
            } else {
                std::thread::yield_now();
            }
        }
    };
    let driver = || loop {
        let _stop = XPanicStop {
            stop: &stop,
            ordered: &ordered,
            written: &written,
        };
        let index = {
            let mut state = ordered.lock().unwrap_or_else(PoisonError::into_inner);
            loop {
                if stop.load(Ordering::Relaxed) || state.taken >= blocks {
                    return;
                }
                if state.taken < state.next + window {
                    break;
                }
                state = written.wait(state).unwrap_or_else(PoisonError::into_inner);
            }
            state.taken += 1;
            state.taken - 1
        };
        x_s8_log("take", index, 0);
        X_BLOCKS.fetch_add(1, Ordering::Relaxed);
        let mut report = Vec::new();
        let mut messages = Vec::new();
        let result = {
            let mut outputs = ReportOutputs {
                formats: vec![FormatOutput {
                    outfmt,
                    sink: OutputSink::Writer(&mut report),
                }],
                diagnostics: &mut messages,
                hits: None,
                observer: None,
            };
            search_block(index, &mut outputs)
        };
        x_s8_log("done", index, report.len());
        let error = result.err();
        if error.is_some() {
            stop.store(true, Ordering::Relaxed);
        }
        let mut state = ordered.lock().unwrap_or_else(PoisonError::into_inner);
        state.done.insert(
            index,
            XDone {
                report,
                messages,
                error,
            },
        );
        state.write_ready();
        x_s8_log("written", state.next, state.done.len());
        if state.failed.is_some() {
            stop.store(true, Ordering::Relaxed);
        }
        drop(state);
        written.notify_all();
    };
    pool.install(|| {
        rayon::scope(|scope| {
            for _ in 0..drivers {
                scope.spawn(|_| {
                    wait_for_every_driver();
                    driver()
                });
            }
        })
    });
    let state = ordered.into_inner().unwrap_or_else(PoisonError::into_inner);
    X_SEARCHES.fetch_add(1, Ordering::Relaxed);
    match state.failed {
        Some((_, error)) => Err(error),
        None => Ok(()),
    }
}

#[cfg(all(test, feature = "parallel", not(target_arch = "wasm32")))]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    /// Runs `x_search_in_slots` on a pool of `threads` with blocks whose search writes
    /// `report(index)` and `message(index)`, block `slow` finishing only after `after` blocks have
    /// finished, and block `fail` (if any) failing after writing its report; returns the outputs and
    /// the result.
    fn run(
        threads: usize,
        slots: usize,
        blocks: usize,
        slow: usize,
        after: usize,
        fail: Option<usize>,
    ) -> (Vec<u8>, Vec<u8>, Result<()>) {
        let finished = AtomicUsize::new(0);
        let mut sink = Vec::new();
        let mut diagnostics = Vec::new();
        let result = crate::utils::threading::with_search_pool(threads, "blastp", |pool| {
            let search = |index: usize, outputs: &mut ReportOutputs<'_>| -> Result<()> {
                if index == slow {
                    while finished.load(Ordering::Acquire) < after {
                        std::hint::spin_loop();
                    }
                }
                let OutputSink::Writer(report) = &mut outputs.formats[0].sink else {
                    unreachable!("a writer of the slots")
                };
                writeln!(report, "row {index}")?;
                writeln!(outputs.diagnostics, "warning {index}")?;
                finished.fetch_add(1, Ordering::Release);
                if fail == Some(index) {
                    anyhow::bail!("block {index} failed");
                }
                Ok(())
            };
            Ok(x_search_in_slots(
                blocks,
                slots,
                "6",
                pool,
                &mut sink,
                &mut diagnostics,
                &search,
            ))
        })
        .expect("the pool");
        (sink, diagnostics, result)
    }

    fn lines(prefix: &str, indices: std::ops::Range<usize>) -> Vec<u8> {
        indices
            .map(|index| format!("{prefix} {index}\n"))
            .collect::<String>()
            .into_bytes()
    }

    #[test]
    fn the_reports_come_in_block_order_whatever_order_the_blocks_finish_in() {
        // Block 0 finishes after three others: the three wait behind it.
        let (sink, diagnostics, result) = run(4, 4, 10, 0, 3, None);
        assert!(result.is_ok());
        assert_eq!(sink, lines("row", 0..10));
        assert_eq!(diagnostics, lines("warning", 0..10));
        // Two drivers: block 0 finishes after the three blocks the window lets the other driver take,
        // which then waits for block 0 to be written.
        let (sink, diagnostics, result) = run(4, 2, 12, 0, 3, None);
        assert!(result.is_ok());
        assert_eq!(sink, lines("row", 0..12));
        assert_eq!(diagnostics, lines("warning", 0..12));
        // One driver per block and more drivers than blocks.
        for (slots, blocks) in [(3, 3), (8, 2)] {
            let (sink, diagnostics, result) = run(4, slots, blocks, usize::MAX, 0, None);
            assert!(result.is_ok());
            assert_eq!(sink, lines("row", 0..blocks));
            assert_eq!(diagnostics, lines("warning", 0..blocks));
        }
    }

    #[test]
    fn a_failed_block_ends_the_search_after_the_blocks_before_it() {
        // Block 2 fails while block 1 is still searched; blocks after 2 that finished are not
        // written, and the failed block's bytes come before its error.
        let (sink, diagnostics, result) = run(4, 4, 12, 1, 2, Some(2));
        assert_eq!(result.unwrap_err().to_string(), "block 2 failed");
        assert_eq!(sink, lines("row", 0..3));
        assert_eq!(diagnostics, lines("warning", 0..3));
        // The first block fails.
        let (sink, _, result) = run(4, 2, 5, usize::MAX, 0, Some(0));
        assert_eq!(result.unwrap_err().to_string(), "block 0 failed");
        assert_eq!(sink, lines("row", 0..1));
    }

    /// An output that takes `limit` bytes and fails the write that would pass them.
    struct Full {
        bytes: Vec<u8>,
        limit: usize,
    }

    impl Write for Full {
        fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
            if self.bytes.len() + buf.len() > self.limit {
                return Err(std::io::Error::other("the output is full"));
            }
            self.bytes.extend_from_slice(buf);
            Ok(buf.len())
        }

        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn a_failed_write_ends_the_search_at_its_block() {
        // Rows of 6 bytes into 20: block 3's row fails, after its messages.
        let mut sink = Full {
            bytes: Vec::new(),
            limit: 20,
        };
        let mut diagnostics = Vec::new();
        let result = crate::utils::threading::with_search_pool(4, "blastp", |pool| {
            let search = |index: usize, outputs: &mut ReportOutputs<'_>| -> Result<()> {
                let OutputSink::Writer(report) = &mut outputs.formats[0].sink else {
                    unreachable!("a writer of the slots")
                };
                writeln!(report, "row {index}")?;
                writeln!(outputs.diagnostics, "warning {index}")?;
                Ok(())
            };
            Ok(x_search_in_slots(
                8,
                2,
                "6",
                pool,
                &mut sink,
                &mut diagnostics,
                &search,
            ))
        })
        .expect("the pool");
        assert_eq!(result.unwrap_err().to_string(), "the output is full");
        assert_eq!(sink.bytes, lines("row", 0..3));
        assert_eq!(diagnostics, lines("warning", 0..4));
    }

    #[test]
    fn the_slot_letters_keep_the_letters_of_two_blocks() {
        assert_eq!(x_slot_letters(50_000, 2), 33_333);
        assert_eq!(x_slot_letters(50_000, 3), 25_000);
        assert_eq!(x_slot_letters(1, 3), 1);
        assert_eq!(x_slot_letters(usize::MAX, 2), usize::MAX / 3);
    }
}
