//! EXPERIMENT (LOSAT_X_QBATCHPAR): blastn's query batches
//! searched side by side, ahead of the batch loop of `run_in_pool`.
//!
//! NCBI's blastn reads its queries in batches (`CBatchSizeMixer`, blastn_app.cpp:261-300): the
//! size of each batch after the first comes from the successful initial extensions
//! (`good_init_extends`, "gie") of the batch before, so batch i+1 cannot be read before batch i
//! is searched. LOSAT's loop does the same and searches one batch at a time on the whole search
//! pool. A search of many queries against one long subject or a few has little parallel work
//! inside a batch (one subject chunk at a time), so most of the pool waits.
//!
//! With the switch on, the loop stays the only place that decides the batch boundaries, the
//! errors, the warnings, the appends and the order of the output. What changes is where a
//! batch's result comes from:
//! - a planner keeps a copy of the loop's mixer, advanced with *predicted* gie (the hit rate per
//!   letter of the last batch of the same size, times the planned batch's letters; 0 for a batch
//!   that will be split into query chunks, as the split batch reports), and plans the keys
//!   `(first record, end record, batch size)` of the next `2 k` batches;
//! - `k` slots search the planned batches ahead of the loop, each with `m` threads (`k x m <= n`):
//!   with `m = 1` on a plain thread with a one-thread pool (`SearchPool::x_serial`), with `m > 1`
//!   on a private pool of `m` threads whose first thread is the slot's (`with_search_pool`, built
//!   once per slot and search), so a batch sees only its own slot's pool;
//! - the loop computes each batch's key from the true mixer, as before, and takes a slot's result
//!   only under that same key (`Driver::take`); after the batch it reports the true gie
//!   (`Driver::commit`) and the planner plans again from the true mixer. A planned key that the
//!   new plan no longer holds is dropped (queued), cancelled (running) or discarded (finished).
//!
//! Why the output is the loop's: a batch's result is a function of its key (its records, its
//! batch size, which enters only an error message) and of the search's fixed inputs; no state
//! crosses batches (`-max_target_seqs`, culling and the hit lists are per query, and a query
//! never spans batches), thread-local scratch is reset per call, and the result does not depend
//! on the thread count. The first batch is never searched ahead (it is the only one that writes
//! to the outputs, `BatchOutputs`), and the loop takes a result only under the key it computed
//! itself, so it takes NCBI's batches in NCBI's order, with their errors at the same batch.
//! A slot's error is a returned value, handled by the loop's error paths when it takes the batch;
//! a slot's panic is kept and resumed only if the loop takes that batch (its message is printed
//! then, not when it happens). That needs unwinding: the release profile has `panic = "abort"`,
//! so there a panic in a batch a slot searches ends the process with its message. To bound that,
//! no slot runs while the loop writes the reports after the loop or what an error path writes:
//! `Driver::finish` cancels and joins every slot before them, and once a slot's batch returns an
//! error the batches planned after it are no longer started. The exception is
//! `LOSAT_X_RPTSTREAM` (`x_rpt_stream.rs`, on by default): there the loop writes each batch's reports
//! when it takes the batch, while the slots search the batches ahead (the slots themselves write
//! nothing), so a slot's panic in a release build ends the process after the reports of the
//! batches taken so far; its error paths still call `finish` before they write. The batches still exposed past an
//! error at batch j: before the error is known, the slots may have searched up to the plan's
//! window (2k - 1 batches) past j; once it is known, at most k - 1 still run, and they are
//! cancelled when the loop meets j. Slot panics are counted (`panicked=`). `LOSAT_TIMING` blocks
//! are text in the batch's result, printed when the loop takes it.
//!
//! Cancellation: a slot's search gets a flag through its pool (`SearchPool::x_with_ext`,
//! `Cancel`), which `search_query_batch` reads at its start, between subjects, between subject
//! chunks of 5 Mb, between query chunks, and before its traceback; set, the search ends with
//! `Cancelled`, which is never taken. The loop's own searches have no flag. The flag is set for a
//! running batch that a new plan drops, when the speculation stops, and in `finish`, so the
//! wait for a cancelled batch is at most one subject chunk's search (or one query chunk's) and
//! its batch's set-up.
//!
//! Memory: each planned batch carries an estimate of its memory while it is searched
//! (`batch_estimate_bytes`), and a slot starts a batch only while the estimates of the batches
//! running plus the result sizes of the finished batches the loop has not taken yet stay within
//! the budget plus the loop's next batch, which is always admitted (the batch the loop needs).
//!
//! Waste: two mispredicted next keys within `GUARD_COMMITS` batches, or more batches discarded
//! (dropped finished results and cancelled searches) than `WASTE_SLOTS x k` plus a
//! `1 / WASTE_SHARE` share of the batches taken, stop the speculation; the loop then searches the
//! rest itself, on the whole pool.
//!
//! Active threads: the slots use `k x m <= n` threads. The loop's thread waits while a slot
//! searches its batch. The loop searches a batch on the whole pool only when no slot is running
//! and none can start (the first batch, and every batch once the search does not speculate).
//!
//! Off (the loop as before): switch off; a build without `parallel`, or wasm32; one thread;
//! a BLAST database (`-db`); `-verbose`, `BLEMIR_DEBUG`, `LOSAT_DEBUG_BLASTN`,
//! `LOSAT_DEBUG_COORDS[_START]`, the blastn trace or `LOSAT_WASI_THREADS_DEBUG`; subjects of n or
//! more work units (5 Mb subject chunks; today's subject-parallel batch already uses the pool);
//! fewer than two batches, or fewer than two batch units, after the first; a memory budget of
//! fewer than two slots; slots that the throughput model does not predict at least 3 % faster
//! than the loop on the whole pool.
//!
//! Slots: per batch in flight `4^w x 4 + 8 x letters` bytes of lookup table (w: the batch's
//! lookup width, two contexts per letter), query-side copies of 13 bytes per query letter
//! (blastn, blastn-short) or 5 (megablast, dc-megablast), and 220 bytes per predicted initial
//! extension for the hits (calibrated to the measured 50 MiB per 1 Mb blastn batch and 65-81 MiB
//! per extra 5 Mb dc-megablast slot). The budget for the batches in flight beyond the loop's one
//! is `LOSAT_X_QBATCHPARMEM` MiB, or by default 96 MiB (blastn, blastn-short) or 28 MiB
//! (megablast, dc-megablast) per thread beyond the first: 0.8 of NCBI BLAST+ 2.17.0's measured
//! growth of max RSS per thread on the audit input. NCBI's max RSS grows with its threads only
//! where its blastn splits the search by queries (`-mt_mode 0` with `-db`: a query file larger
//! than 10,000,000 bytes, or 1,500,000 at a word size below 10, and a database shorter than
//! 6,000,000,000 letters, or 150,000,000; `ncbi_splits_by_queries`); elsewhere its max RSS stays
//! that of one thread, so the default budget is 0 there (no slots) and LOSAT's max RSS keeps to
//! that of the loop (`Profile::of_search`; the subjects of `-subject` stand for the database).
//! The slots come from a throughput model
//! (`choose_slots`): the units of work are the planned batches of at least half the
//! letter-weighted median planned batch's letters (`batch_units`); for
//! `k = 1 ..= min(1 + budget / estimate, n, units)` and `m = n / k` the predicted span is
//! rounds(k) x contention(k) / (m x E(m)), with E(m) the measured efficiency of one batch on m
//! threads, measured up to 16 and with the speedup kept flat beyond (`Profile`,
//! `predicted_span`); `k = 1` is the loop on the whole pool. Slots start only when the best
//! `k >= 2` predicts less than 0.97 of the loop's span. `LOSAT_X_QBATCHPARSLOTS=k` (tuning) sets k
//! directly.
//!
//! `LOSAT_X_DCTWOSLOTS` (a strict-set switch: on unless it is `0`): for `-task dc-megablast` only
//! (megablast with template options runs the same stages but keeps the model), where the model
//! takes three slots or more (k x m), two slots of `n / 2` threads search instead when the
//! two-slot rule holds (`SlotRule::TwoSlots`, `two_slots_win`, `choose_slots`): 12 to 33 threads
//! (the measured range: at most 16 threads per slot; 64 and 128 threads keep the model), an even
//! count of batch units (an odd count leaves the last round to one slot), and either the model's
//! slots are narrow and leave threads idle (`m <= 3` and `k x m < 2 x (n / 2)`), or all units fit
//! in its first round (`units <= k`), or its last round is at most half full within five rounds
//! (`0 < units mod k <= k / 2` and `ceil(units / k) <= 5`). Where it takes two slots or keeps the
//! loop, or the rule declines, nothing changes. The E(m)
//! points predate LOSAT_X_PAIRPARDISCARR and LOSAT_X_PAIRPARHOME (2026-10-20 and 2026-10-24),
//! which made one dc-megablast batch scale better on 8 to 32 threads; the rule is drawn from a map
//! of `LOSAT_X_QBATCHPARSLOTS=2` against the model on dev d92512d0 (ws, 2026-10-25, line D
//! session 14 window w2: dc-megablast of the first 6,900 / 10,400 / 15,600 / 20,800 / 31,168
//! D. melanogaster cDNAs and of 5,000 cDNAs against 3 Mb of chromosome 2L, `-subject`, at 12 / 16 /
//! 24 / 32 threads, 4 rounds; `the_two_slot_rule_holds_only_where_two_slots_won_on_the_w2_map`):
//! two slots win where the model's narrow slots leave threads idle (31,168 cDNAs at 32 threads,
//! 9 x 3 on 27 threads: x0.909; 15,600 at 24 threads, 9 x 2 on 18: x0.892), where its batches fit
//! in one round (10,400 at 24 / 32 threads, 6 x 4 / 6 x 5: x0.946 / x0.943) and where its last
//! round is half empty in a short search (31,168 at 16 threads, 4 x 4 for 18 units: x0.983); they
//! lose with an odd unit count (5,000 cDNAs, 3 units: x1.10-1.16 at every thread count) and where
//! the model runs two full rounds or more of 4-thread slots (20,800 cDNAs at 12 / 16 / 24
//! threads: x1.006-1.018; 31,168 at 12 / 24 threads: x1.016 / x1.031). The rule takes every
//! winning cell of that map but five: 15,600 cDNAs at 24 / 32 threads (9 units, odd), and 10,400
//! at 12 / 16 threads and 20,800 at 32 threads (two full rounds of 4- and 5-thread slots); and
//! none of the losing ones. Narrow slots that use every thread, and last rounds more than half
//! full, were not measured, so the rule keeps the model there (review 2). With the final rule
//! (ws, 2026-10-25, line D session 15 windows w3 to w5 against dev 416c427f and e6429b65, 6 to
//! 8 rounds, max RSS by `time %M`): 31,168 cDNAs at 16 threads x0.976-0.982 of dev's time, max
//! RSS 394-398 -> 222-235 MiB (NCBI `-subject`, one thread: 313 MiB); at 32 threads
//! x0.911-0.916, 824 -> 330-335 MiB; 6,900 at 16 threads x0.965-0.980, 347 -> 208-209 MiB;
//! 10,400 at 24 threads x0.967-0.973, 510 -> 233-238 MiB. Where the rule declines (5,000,
//! 15,600 and 20,800 cDNAs at 16 threads, 31,168 at 24) the model's layout runs, with dev's
//! cycles and instructions. With the rule `rule_k`, `rule_m` and `predicted` of the summary
//! line still state the model's layout, and `k`, `m` the slots that ran. The tuning variable sets
//! k directly, as before.
//!
//! `LOSAT_X_BNEIGHTSLOTS` (a strict-set switch: on unless it is `0`): for `-task blastn` only
//! (blastn-short was not measured), where the model takes more than eight slots on 12 to 64
//! threads, its best layout of at most eight slots runs instead (`SlotRule::EightSlots`,
//! `eight_slots_apply`, `choose_slots`: the least `predicted_span` of 2 to 8 slots, the fewest on a
//! tie; it always pays the 0.97 gate there). Each slot holds a batch's lookup table, query copies
//! and hits (about 38 MiB for a 1 Mb blastn batch, measured below), so the slot count sets the memory: for
//! `-subject` searches of many cDNAs NCBI's max RSS stays that of one thread, and the model's 16
//! one-thread slots at 16 threads held LOSAT above it. Where the model takes eight slots or fewer,
//! on fewer than 12 or more than 64 threads, and with the tuning variable, nothing changes. The
//! rule comes from slot counts of `LOSAT_X_QBATCHPARSLOTS` against the model on dev 51234c2a (ws,
//! 2026-10-25, line D session 15 window w6: blastn of the first 31,168 D. melanogaster cDNAs
//! (89 batch units) and of 5,000 (14 units) against 3 Mb of chromosome 2L, `-subject`, 8 / 16 / 32
//! threads, 4-6 rounds): 31,168 cDNAs at 16 threads, model 16 x 1 (max RSS 702 MiB; NCBI
//! `-subject`, one thread: 630 MiB), 8 x 2 x0.981 (400 MiB), 4 x 4 x1.022 (232), 2 x 8 x1.184
//! (149), the loop x1.274 (107); at 32 threads, model 32 x 1 (1,439 MiB), 8 x 4 x0.932 (515), 4 x 8
//! x0.994 (350); 5,000 cDNAs at 16 / 32 threads, model 14 x 1 / 14 x 2 (560 / 619 MiB), 8 slots
//! x0.891 / x0.878 (363 / 380), 4 x0.883 / x0.872 (225 / 240). Windows w1 and w3 (the same binary,
//! session 16, 4-6 rounds) added 6 / 10 / 12 / 16 slots, 12 / 24 / 64 threads, 4,000 / 6,900 /
//! 15,600 cDNAs and 10 / 30 Mb subjects (chromosome 2L, then 2R): 31,168 cDNAs at 12 threads 6 x 2
//! x0.991 of the model's 12 x 1 (307 / 542 MiB), at 24 threads 8 x 3 x0.946 (427 / 1,041), at 64
//! threads 8 x 8 x0.874 of 32 x 2 (575 / 1,563); slots that leave threads idle lose (16 threads: 6
//! x 2 x1.237, 10 x 1 x1.448). Of the 16 inputs and thread counts where the rule applies
//! (4,000-31,168 cDNAs, 3-30 Mb, 12-64 threads), 10 were measured with the rule's own layout:
//! x0.755-0.991 of the model's time, within x1.054 of the best measured count (4,000 cDNAs (12
//! units) at 16 threads 4 x 4 x0.755, 214 / 452 MiB; 6,900 cDNAs x 10 Mb at 32 threads 8 x 4
//! x0.841, 661 / 686); for 14 units the rule runs 7 x m (a tie with 8 x m, the fewest slots),
//! measured there through 8 x m (x0.768-0.901); 4,000 cDNAs at 32 threads (6 x 5) was not measured.
//! Refits of the E(m) points (with contention, Amdahl's serial share, or a loop-thread cost per
//! batch) to these 145 cells missed the best count by up to x1.07 (the best form) to x1.36, so the
//! points stay and the rule caps the slots. Where the model takes eight slots or fewer it can still
//! miss (5,000 cDNAs x 30 Mb at 16 threads: 7 x 2 at x1.320 of 2 x 8's time); the rule leaves those
//! layouts as they are.
//! With the rule `rule_k`, `rule_m` and `predicted` of the summary line still state the model's
//! layout, and `k`, `m` the slots that ran.
//!
//! Diagnostics: with the switch and `LOSAT_X_STATS` set, one line on standard error at the end of
//! the search, `[X_STATS] QBATCHPAR batches= precomputed= loop= mispredicted= discarded=
//! cancelled= panicked= waste_bound= k= m= idle= units= rule_k= rule_m= predicted= budget_mib=
//! estimate_mib= peak_mib= off=` (precomputed: batches taken from a slot; loop: batches the loop
//! searched itself; discarded: finished or cancelled batches dropped, `cancelled` of them
//! cancelled; waste_bound: the discards that stop the speculation, at the end; idle: threads no
//! slot uses; units: the batch units after the first batch, 0 when no batch follows it; rule_k,
//! rule_m: the model's best layout of two slots or more, or the layout `LOSAT_X_QBATCHPARSLOTS`
//! sets, 0 when the model did not run; predicted: that layout's predicted span over the loop's,
//! printed to three decimals, the value the 0.97 gate compares (`slots_pay`), `-` when the model
//! did not run; peak_mib: the most memory admitted at once; off: why the search did not, or stopped to,
//! speculate, `none` if it never stopped).
//! Without `LOSAT_X_STATS` the switch writes nothing of its own. `LOSAT_X_QBATCHPARMISS=N` (test
//! only) corrupts every N-th prediction so that the discard, replan and guard paths run.
//!
//! NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:277-300
//! ```c++
//!         for (; !input.End(); formatter.ResetScopeHistory(), QueryBatchCleanup() ) {
//!             CRef<CBlastQueryVector> query_batch(input.GetNextSeqBatch(*scope));
//!     ...
//!                 results = lcl_blast.Run();
//!                 if (!batch_size)
//!                     input.SetBatchSize(mixer.GetBatchSize(lcl_blast.GetNumExtensions()));
//! ```

#![cfg_attr(
    not(all(feature = "parallel", not(target_arch = "wasm32"))),
    allow(dead_code, unused_imports, unused_variables)
)]

use std::collections::{HashMap, HashSet, VecDeque};
use std::ops::Range;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::{Arc, Condvar, Mutex, MutexGuard, OnceLock, PoisonError};

use crate::algorithm::blastn::query_split::{calculate_num_chunks, SplitSizes};
use crate::blastinput::query_batch::BatchSizeMixer;
use crate::blastinput::seq_range::{next_ranged_batch, QueryInput};
use crate::utils::threading::SearchPool;

/// Bytes in a MiB.
pub(super) const MIB: u64 = 1 << 20;
/// Query-side bytes per query letter of a blastn or blastn-short batch in flight (the batch's
/// copy of its records, reverse complements, contexts with and without sentinels, masks).
const QUERY_BYTES_PER_LETTER: u64 = 13;
/// The same for megablast and dc-megablast: each extra 5 Mb dc-megablast slot measured 65-81 MiB
/// after #217's query-copy cuts (perf-mb split3, 2026-10-17: Ad 4 threads 2 slots 303 MiB against
/// 229 in the loop, Ad 16 threads 4 slots 502 against 258, (b) 4 threads 2 slots 176 against
/// 110), where 13 bytes per letter estimate 117 MiB and 5 estimate 79.
const MEGABLAST_QUERY_BYTES_PER_LETTER: u64 = 5;
/// Hit bytes per predicted successful initial extension of a batch in flight.
const HIT_BYTES_PER_EXTENSION: u64 = 220;
/// Mispredicted next keys within `GUARD_COMMITS` batches that stop the speculation.
const GUARD_COMMITS: u64 = 8;
/// Discarded batches that stop the speculation: `WASTE_SLOTS x k` (a plan's window) plus one
/// `WASTE_SHARE`-th of the batches taken from the slots.
const WASTE_SLOTS: u64 = 2;
const WASTE_SHARE: u64 = 4;
/// The stack of a slot's thread: that of the main thread, where the loop's searches (and their
/// top level) ran before.
const SLOT_STACK_BYTES: usize = 8 << 20;

/// No NCBI counterpart: the switches of this path, read once per search; they choose a LOSAT
/// code path and change no value NCBI computes.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Switches {
    /// `LOSAT_X_QBATCHPAR` (a strict-set switch: on unless it is `0`).
    pub on: bool,
    /// `LOSAT_X_STATS` (any value): the summary line.
    pub stats: bool,
    /// `LOSAT_X_QBATCHPARMEM`: the budget of the batches in flight beyond the loop's, in MiB.
    pub mem_mib: Option<u64>,
    /// `LOSAT_X_QBATCHPARSLOTS` (tuning): the number of slots, at most the thread count.
    pub slots: Option<usize>,
    /// `LOSAT_X_QBATCHPARMISS` (test only): corrupt every N-th prediction.
    pub miss_every: Option<u64>,
    /// `LOSAT_X_DCTWOSLOTS` (a strict-set switch: on unless it is `0`): the profile's two-slot
    /// rule (`SlotRule::TwoSlots`, `two_slots_win`).
    pub dc_two: bool,
    /// `LOSAT_X_BNEIGHTSLOTS` (a strict-set switch: on unless it is `0`): the profile's rule of
    /// at most eight slots (`SlotRule::EightSlots`, `eight_slots_apply`).
    pub bn_eight: bool,
}

impl Switches {
    /// The switches as the environment sets them. A value variable that is not a number (or is
    /// 0 for the slots and the corruption) counts as not set.
    pub(super) fn from_env() -> Self {
        let number = |name: &str| -> Option<u64> {
            std::env::var(name)
                .ok()
                .and_then(|value| value.trim().parse::<u64>().ok())
        };
        Switches {
            on: crate::utils::xswitch::on("LOSAT_X_QBATCHPAR"),
            stats: std::env::var_os("LOSAT_X_STATS").is_some(),
            mem_mib: number("LOSAT_X_QBATCHPARMEM"),
            slots: number("LOSAT_X_QBATCHPARSLOTS")
                .filter(|&slots| slots > 0)
                .map(|slots| usize::try_from(slots).unwrap_or(usize::MAX)),
            miss_every: number("LOSAT_X_QBATCHPARMISS").filter(|&every| every > 0),
            dc_two: crate::utils::xswitch::on("LOSAT_X_DCTWOSLOTS"),
            bn_eight: crate::utils::xswitch::on("LOSAT_X_BNEIGHTSLOTS"),
        }
    }
}

/// Why a search does not speculate, or stopped to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Off {
    /// The switch is off (`LOSAT_X_QBATCHPAR=0`).
    NotSet,
    /// A build without the `parallel` feature, or wasm32.
    Build,
    /// One thread.
    Threads,
    /// A BLAST database (`-db`).
    Database,
    /// Verbose, debug or trace output.
    Debug,
    /// Subjects of as many work units as threads.
    Subjects,
    /// The search ended before its first batch was searched.
    Unstarted,
    /// Fewer than two batches, or fewer than two batch units, after the first.
    Single,
    /// A budget of fewer than two slots (by default, a query file that NCBI's blastn does not
    /// split by queries: `Profile::of_search`).
    Budget,
    /// The throughput model predicts no gain: its best layout of two slots or more needs 0.97 of
    /// the loop's span or more.
    Model,
    /// A slot's thread or pool could not be started.
    Spawn,
    /// Two mispredicted next keys within `GUARD_COMMITS` batches.
    Guard,
    /// More discarded batches than the waste bound.
    Waste,
    /// The loop's next key was not planned (not expected to happen).
    Unplanned,
    /// EXPERIMENT (LOSAT_X_BNQDROP): the loop takes its batches' records out of their store.
    Taken,
}

impl Off {
    fn name(self) -> &'static str {
        match self {
            Off::NotSet => "not_set",
            Off::Build => "build",
            Off::Threads => "threads",
            Off::Database => "db",
            Off::Debug => "debug",
            Off::Subjects => "subjects",
            Off::Unstarted => "unstarted",
            Off::Single => "single",
            Off::Budget => "budget",
            Off::Model => "model",
            Off::Spawn => "spawn",
            Off::Guard => "guard",
            Off::Waste => "waste",
            Off::Unplanned => "unplanned",
            Off::Taken => "taken",
        }
    }
}

/// The key of a query batch: its input records `input_start..input_end` and the batch size the
/// loop read them with (it enters only an error message). The loop computes it from the true
/// mixer; a result is used only under the key it was searched for.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) struct BatchKey {
    pub input_start: usize,
    pub input_end: usize,
    pub batch_size: u32,
}

/// A batch for a slot: its key, its searched queries, and its estimated memory while it is
/// searched (`batch_estimate_bytes`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct Job {
    pub key: BatchKey,
    pub searched: Range<usize>,
    pub bytes: u64,
}

/// The query input as the batch loop reads it.
#[derive(Clone, Copy)]
pub(super) struct PlanInput<'a> {
    pub query_input: &'a QueryInput,
    /// The length of each searched query.
    pub lengths: &'a [usize],
    /// The reading ends with a reader error after the last record (`QueryReading::Error`).
    pub reading_error: bool,
    /// `BATCH_SIZE` is set: every batch has that size.
    pub fixed_batch_size: bool,
    pub split: SplitSizes,
}

/// A planned batch.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Planned {
    pub key: BatchKey,
    pub searched: Range<usize>,
    /// Its searched letters, and its longest query.
    pub letters: u64,
    pub max_len: usize,
    /// The predicted successful initial extensions.
    pub gie: u64,
}

/// No NCBI counterpart: the speculation of NCBI's `CBatchSizeMixer` (scheduling only).
#[derive(Debug, Default)]
pub(super) struct Planner {
    /// The hit rate per letter of the last committed batch of each size.
    rates: HashMap<u32, f64>,
    /// That of the last committed batch.
    last_rate: Option<f64>,
    /// Predictions made (for `miss_every`).
    predictions: u64,
    miss_every: Option<u64>,
}

impl Planner {
    pub(super) fn new(miss_every: Option<u64>) -> Self {
        Planner {
            miss_every,
            ..Planner::default()
        }
    }

    /// A committed batch of `batch_size` with `letters` searched letters had `gie`.
    pub(super) fn record(&mut self, batch_size: u32, gie: u64, letters: u64) {
        let rate = gie as f64 / letters.max(1) as f64;
        self.rates.insert(batch_size, rate);
        self.last_rate = Some(rate);
    }

    /// The predicted gie of a batch of `batch_size` and `letters`; 0 for a batch that NCBI
    /// splits into query chunks (its own counter is not used; the chunks have theirs).
    fn predict(&self, batch_size: u32, letters: u64, split: bool) -> u64 {
        if split {
            return 0;
        }
        let rate = self
            .rates
            .get(&batch_size)
            .copied()
            .or(self.last_rate)
            .unwrap_or(0.0);
        (rate * letters as f64) as u64
    }

    /// Up to `window` batches from input record `start`, the first of `batch_size` (the true
    /// next size), the next ones of the sizes that `mixer` (the loop's mixer after the last
    /// batch; `None` with `BATCH_SIZE`) gives for the predicted gie. Stops after the batch that
    /// reaches the last record, and before a batch that ends the run with an error in the loop
    /// (only skipped records, only queries without letters, or the reader error after the last
    /// record). `speculate` counts the predictions (and corrupts every `miss_every`-th).
    pub(super) fn plan(
        &mut self,
        input: &PlanInput<'_>,
        start: usize,
        batch_size: u32,
        mixer: Option<&BatchSizeMixer>,
        window: usize,
        speculate: bool,
    ) -> Vec<Planned> {
        let records = input.query_input.input_lengths.len();
        let mut mixer = mixer.cloned();
        let mut planned = Vec::new();
        let (mut input_start, mut size) = (start, batch_size);
        while planned.len() < window && input_start < records {
            let (input_end, reached_size) = next_ranged_batch(
                &input.query_input.input_lengths,
                &input.query_input.skipped,
                input_start,
                size,
            );
            if input_end == records && !reached_size && input.reading_error {
                break;
            }
            let searched = input.query_input.searched_in(input_start..input_end);
            let lengths = &input.lengths[searched.clone()];
            if lengths.iter().all(|&length| length == 0) {
                break;
            }
            let letters: u64 = lengths.iter().map(|&length| length as u64).sum();
            let max_len = lengths.iter().copied().max().unwrap_or(0);
            let split = calculate_num_chunks(input.split, letters as usize).0 > 1;
            let gie = self.predict(size, letters, split);
            planned.push(Planned {
                key: BatchKey {
                    input_start,
                    input_end,
                    batch_size: size,
                },
                searched,
                letters,
                max_len,
                gie,
            });
            if input_end == records || planned.len() == window {
                break;
            }
            // The loop's `mixer.batch_size(Some(gie as u32 as i32))`.
            let mut next = match mixer.as_mut() {
                Some(mixer) if !input.fixed_batch_size => {
                    mixer.batch_size(Some(gie as u32 as i32)) as u32
                }
                _ => size,
            };
            if speculate {
                self.predictions += 1;
                if self
                    .miss_every
                    .is_some_and(|every| self.predictions % every == 0)
                {
                    next = if next > 1 { next / 2 } else { 2 };
                }
            }
            input_start = input_end;
            size = next;
        }
        planned
    }
}

/// The default budget per thread beyond the first, in MiB (0.8 x NCBI BLAST+ 2.17.0's max-RSS
/// growth per thread on the audit input: blastn (1888 - 127) / 15, dc-megablast (679 - 150) /
/// 15; megablast uses dc-megablast's).
pub(super) fn allowance_mib(task: &str) -> u64 {
    match task {
        "megablast" | "dc-megablast" => 28,
        _ => 96,
    }
}

/// No NCBI counterpart (scheduling only): the rule a task applies to the slot model's layout
/// while the rule's switch is on (module header, `choose_slots`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) enum SlotRule {
    /// The model's layout: no rule, or its switch is `0`.
    #[default]
    Model,
    /// `LOSAT_X_DCTWOSLOTS` (`-task dc-megablast`): two slots of half the threads where the
    /// model takes more and `two_slots_win` holds.
    TwoSlots,
    /// `LOSAT_X_BNEIGHTSLOTS` (`-task blastn`): where the model takes more than `BN_SLOTS` slots
    /// and `eight_slots_apply` holds, its best layout of at most `BN_SLOTS` slots.
    EightSlots,
}

/// No NCBI counterpart (scheduling only): what the slot choice knows of a task.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Profile {
    /// The default budget per thread beyond the first, in MiB (`allowance_mib`).
    pub allowance_mib: u64,
    /// Query-side bytes per query letter of a batch in flight.
    pub query_bytes_per_letter: u64,
    /// The in-batch efficiency E(m) (`EFFICIENCY_MEGABLAST`, `EFFICIENCY_BLASTN`).
    pub efficiency: &'static [(usize, f64)],
    /// The rule applied to the model's layout while its switch is on: `-task dc-megablast` the
    /// two-slot rule of `LOSAT_X_DCTWOSLOTS`, `-task blastn` the eight-slot rule of
    /// `LOSAT_X_BNEIGHTSLOTS` (module header).
    pub rule: SlotRule,
}

impl Profile {
    /// The profile of `-task` `task`: megablast and dc-megablast (5 Mb batches) take the
    /// efficiency measured on dc-megablast; blastn and blastn-short that of blastn. Only
    /// dc-megablast takes the two-slot rule, and only blastn the eight-slot rule (module header).
    pub(super) fn of_task(task: &str) -> Self {
        let megablast = matches!(task, "megablast" | "dc-megablast");
        Profile {
            allowance_mib: allowance_mib(task),
            query_bytes_per_letter: if megablast {
                MEGABLAST_QUERY_BYTES_PER_LETTER
            } else {
                QUERY_BYTES_PER_LETTER
            },
            efficiency: if megablast {
                EFFICIENCY_MEGABLAST
            } else {
                EFFICIENCY_BLASTN
            },
            rule: match task {
                "dc-megablast" => SlotRule::TwoSlots,
                "blastn" => SlotRule::EightSlots,
                _ => SlotRule::Model,
            },
        }
    }

    /// No NCBI counterpart (scheduling only): the profile of a search of `-task` `task` with
    /// word size `word_size`, a query file of `query_file_bytes` (`None`: standard input or no
    /// file) and subjects of `subject_letters` letters: that of `of_task`, with no default budget
    /// (`allowance_mib` 0, so no slots unless `LOSAT_X_QBATCHPARMEM` sets a budget) where NCBI's
    /// blastn would not split the search by queries (`ncbi_splits_by_queries`). There NCBI's max
    /// RSS does not grow with its threads (perf-mb defbundle, 2026-10-17: plain blastn of 1,775
    /// cDNAs, a 5.4 MB query file, against 3 Mb: NCBI `-db` 123-136 MiB at 1-16 threads, while 5
    /// slots at 16 threads took LOSAT from 122 to 332 MiB), so slots there would hold LOSAT above
    /// it. The slots change no output; this chooses only whether they run.
    pub(super) fn of_search(
        task: &str,
        word_size: usize,
        query_file_bytes: Option<u64>,
        subject_letters: u64,
    ) -> Self {
        let mut profile = Profile::of_task(task);
        if !ncbi_splits_by_queries(word_size, query_file_bytes, subject_letters) {
            profile.allowance_mib = 0;
        }
        profile
    }
}

/// NCBI's cut-offs of its split by queries for blastn: the database letters below which, and
/// the query file bytes above which, `-mt_mode 0` with more than one thread searches the queries
/// split between the threads, by word size.
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:96-114
/// ```c++
/// static void s_GetMT1Cutoffs(const int word_size, Int8 & max_db_size, Int8 & min_q_size)
/// {
///     const Int8 kMT1MaxDbSz_Blastn = 150000000; // bases
///     const Int8 kMT1MinQSz_Blastn = 1500000; //bytes
///     // Choose a bigger min quer size so we can support larger DB
///     const Int8 kMT1MaxDbSz_Megablast = 6000000000;
///     const Int8 kMT1MinQSz_Megablast = 10000000;
///     if (word_size < 10) {
///         // CompressedAa
///         max_db_size = kMT1MaxDbSz_Blastn;
///         min_q_size = kMT1MinQSz_Blastn;
///     }
///     else {
///         max_db_size = kMT1MaxDbSz_Megablast;
///         min_q_size = kMT1MinQSz_Megablast;
///     }
/// }
/// ```
pub(super) fn ncbi_mt1_cutoffs(word_size: usize) -> (u64, u64) {
    if word_size < 10 {
        (150_000_000, 1_500_000)
    } else {
        (6_000_000_000, 10_000_000)
    }
}

/// Whether NCBI's blastn, with `-mt_mode 0` (its default) and more than one thread, splits the
/// search of a query file of `query_file_bytes` against a database of `subject_letters` letters
/// by queries (each thread then searches its own query batches, and its max RSS grows with the
/// threads). A query that is not a file of known size (`None`: standard input) is not split.
/// LOSAT uses it only to choose `LOSAT_X_QBATCHPAR`'s default budget (`Profile::of_search`).
///
/// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:134-153
/// ```c++
///     int num_threads = m_CmdLineArgs->GetNumThreads();
///     int mt_mode = m_CmdLineArgs->GetMTMode();
///     if (!m_CmdLineArgs->ExecuteRemotely() && (num_threads > 1) &&
///        (mt_mode != CMTArgs::eSplitByDB)){
/// ...
///             s_GetMT1Cutoffs(word_size, max_db_size, min_q_size);
///             CRef<CSeqDB> seqdb = sdb->GetSeqDb();
///             Uint8 total_length = seqdb->GetTotalLength();
///             if (mt_mode == CMTArgs::eSplitAuto){
///                 if (total_length < max_db_size) {
///                     if (args.Exist(kArgQuery) && args[kArgQuery].HasValue()) {
///                         CFile file( args[kArgQuery].AsString());
///                         if (file.GetLength() > min_q_size) {
///                             m_UsageReport.AddParam(CBlastUsageReport::eMTMode, CMTArgs::eSplitByQueries);
///                             return x_RunMTBySplitQuery();
/// ```
pub(super) fn ncbi_splits_by_queries(
    word_size: usize,
    query_file_bytes: Option<u64>,
    subject_letters: u64,
) -> bool {
    let (max_db_size, min_q_size) = ncbi_mt1_cutoffs(word_size);
    subject_letters < max_db_size && query_file_bytes.is_some_and(|bytes| bytes > min_q_size)
}

/// The estimated memory of one batch in flight: lookup table of width `lut_width` over the two
/// contexts of `letters`, the query-side copies of `query_bytes_per_letter` per letter
/// (`Profile`), and the hits of `gie` initial extensions.
pub(super) fn batch_estimate_bytes(
    lut_width: usize,
    letters: u64,
    gie: u64,
    query_bytes_per_letter: u64,
) -> u64 {
    let table = 4u64.saturating_pow(lut_width as u32).saturating_mul(4);
    table
        .saturating_add(letters.saturating_mul((2 * 4u64).saturating_add(query_bytes_per_letter)))
        .saturating_add(gie.saturating_mul(HIT_BYTES_PER_EXTENSION))
}

// No NCBI counterpart (scheduling only): the parameters of the slot model, all in one place.
//
// E(m) = the speedup of one batch on m threads over one thread, divided by m, as (m, E) points,
// log2-linear between them. Measured up to 16 threads per batch, on one machine with the SMT
// siblings idle: the workstation (perf-mb split3 window 1, 2026-10-17 10:45-10:58: 2x Xeon Gold
// 6338, node 1 CPUs 32-47, D. melanogaster cDNAs against 3 Mb of chromosome 2L with `-subject`,
// every search with LOSAT_X_PAIRPARDISC, DIAGLV, PAIRPARFOLD, QDUSTPAR and DISCLUTPAR), from the
// batch times of `logs/split3_report.md` section 1c, dc-megablast 5 Mb batches (inputs (a), (b)
// and Ad: 1,775, 3,492 and 31,168 cDNAs) and blastn 1 Mb batches (inputs (c) and Ab: 1,775 and
// 31,168 cDNAs). Per point (reconstructed from that table; batch time at 1 thread: dc 1.812 s,
// Ad alone, blastn 0.804 s, (c) alone): dc E(2) from Ad 4 threads in 2 slots (beside one other
// slot), E(4) from Ad 4 threads in the loop (alone), E(8) from (b) 16 threads in 2 slots (beside
// one other slot), E(16) from the loop at 16 threads (alone); blastn E(2) from (c) 4 threads in 2
// slots (beside one), E(4) and E(16) from (c) in the loop at 4 and 16 threads (alone), E(8)
// between (c) 8 threads in the loop and 16 threads in 2 slots. So the points measured beside a
// slot already hold some contention, which the contention term then counts again for those m.
// Beyond the last point the speedup is kept flat (E(m) = E(16) x 16 / m): nothing was measured
// there, and blastn's speedup grew only 1.16-fold from 8 to 16 threads, so the model must not
// assume that the loop keeps scaling on 32-128 threads (slots are exact and stay within the
// budget, so erring towards slots costs memory, not output). On machines whose `-num_threads`
// counts SMT siblings, E(m) for large m is lower still.
//
// Contention: a batch runs `CONTENTION_PER_SLOT` slower per further slot beside it, up to
// `CONTENTION_SLOTS` slots (the same window: with 4 dc-megablast slots the seed stage, 58-66 % of
// a batch, took 25-34 % longer than alone, 18 % more from 2 to 8 slots). It was measured on
// dc-megablast; blastn's is larger (a 1 Mb batch in 5 slots of 3 threads took 0.51-0.53 s at
// (c) 16 threads, against the model's 0.42), which changes no measured cell's choice.
//
// The model picked the measured best layout of every cell of that report's section 3 (window
// 2: (a) 4 / 16 threads off, (b) 4: 2 x 2, (b) 16: off, (c) 4: off, (c) 8: 2 x 4, (c) 16: 5 x 3,
// Ad 4: 2 x 2, Ad 16: 4 x 4, Ab 4: 4 x 1, Ab 16: 16 x 1; `cells_of_the_split3_report`).
// Sensitivity (perf-mb `scripts/qbprule_model.py sens`): every E(m >= 2) times one factor in
// 0.8 ..= 1.2 changes three choices: (c) 8 threads 2 x 4 -> 5 x 1 below 0.93, Ab 4 threads
// 4 x 1 -> 2 x 2 above 1.06, Ab 16 threads 16 x 1 -> 8 x 2 above 1.12. One point times 0.8 or
// 1.2 alone changes (b) and Ad at 4 threads (2 x 2 -> off with E(2) x 0.8 or E(4) x 1.2: their
// predicted 0.956 is near the 0.97 gate), and for some points (b) 16, (c) 4 and 8, Ad 16 and
// Ab 4 and 16; never (a) or (c) 16. A contention of 0.075 per slot turns (b) and Ad at 4 threads
// off; one of 0 gives (b) 16 threads 2 x 8 and (c) 8 threads 5 x 1.

/// E(m) points of megablast and dc-megablast (measured on dc-megablast; see above).
const EFFICIENCY_MEGABLAST: &[(usize, f64)] =
    &[(1, 1.0), (2, 0.89), (4, 0.81), (8, 0.51), (16, 0.49)];
/// E(m) points of blastn and blastn-short (measured on blastn; see above).
const EFFICIENCY_BLASTN: &[(usize, f64)] = &[(1, 1.0), (2, 0.85), (4, 0.71), (8, 0.50), (16, 0.29)];
/// A batch's slowdown per further slot beside it, and the slots beyond which it grows no more.
const CONTENTION_PER_SLOT: f64 = 0.05;
const CONTENTION_SLOTS: usize = 5;
/// Slots start only when the model predicts their span below this share of the loop's.
const MODEL_GAIN: f64 = 0.97;
/// The two-slot rule of dc-megablast (`LOSAT_X_DCTWOSLOTS`, module header, `two_slots_win`): the
/// slots it runs where the model takes more, the fewest threads and the most threads per slot
/// it was measured with, the slot width of the model's layout at or below which it holds (where
/// that layout leaves threads idle), and the rounds within which a last round at most half full
/// makes it hold.
const DC_SLOTS: usize = 2;
const DC_MIN_THREADS: usize = 12;
const DC_SLOT_THREADS: usize = 16;
const DC_NARROW_SLOT: usize = 3;
const DC_SHORT_SEARCH_ROUNDS: usize = 5;
/// The eight-slot rule of blastn (`LOSAT_X_BNEIGHTSLOTS`, module header, `eight_slots_apply`):
/// the most slots it lets run, and the fewest and the most threads it was measured with.
const BN_SLOTS: usize = 8;
const BN_MIN_THREADS: usize = 12;
const BN_MAX_THREADS: usize = 64;

/// A predicted ratio in thousandths, as the summary line prints it (`predicted=`).
fn milli(ratio: f64) -> u64 {
    (ratio * 1000.0).round() as u64
}

/// Whether a layout of this predicted span over the loop's is worth its slots: its ratio as the
/// summary line prints it, in thousandths, strictly below `MODEL_GAIN` (so the printed value
/// decides: `predicted=0.970` keeps the loop, `0.969` takes the slots).
fn slots_pay(ratio: f64) -> bool {
    milli(ratio) < milli(MODEL_GAIN)
}

/// The reference batch of planned batches of `letters` each: the letter-weighted median, the
/// batch at which the letters of the batches sorted by size reach half of all letters (0 for
/// none). One oversized batch (a long record, which never spans batches) does not move it, as
/// the largest batch would.
fn reference_letters(letters: &[u64]) -> u64 {
    let mut sorted = letters.to_vec();
    sorted.sort_unstable();
    let total: u64 = sorted
        .iter()
        .fold(0u64, |sum, &batch| sum.saturating_add(batch));
    let mut sum = 0u64;
    for batch in sorted {
        sum = sum.saturating_add(batch);
        if sum.saturating_mul(2) >= total {
            return batch;
        }
    }
    0
}

/// The batch units of planned batches of `letters` each: the batches of at least half the
/// reference batch's letters (`reference_letters`), each one unit, oversized ones too. A tiny
/// batch of the mixer's alternation (100 letters beside 5 Mb ones) or a short last batch is no
/// unit of work.
pub(super) fn batch_units(letters: &[u64]) -> usize {
    let reference = reference_letters(letters);
    letters
        .iter()
        .filter(|&&batch| batch.saturating_mul(2) >= reference)
        .count()
}

/// E(m) of `points` (`EFFICIENCY_*`): log2-linear between the points; beyond the last point the
/// speedup `m x E(m)` stays at the last point's (see the parameters above).
fn efficiency(points: &[(usize, f64)], m: usize) -> f64 {
    let m = m.max(1);
    let x = (m as f64).log2();
    for pair in points.windows(2) {
        let ((m0, e0), (m1, e1)) = (pair[0], pair[1]);
        let (x0, x1) = ((m0 as f64).log2(), (m1 as f64).log2());
        if x <= x1 {
            return e0 + (e1 - e0) * (x - x0) / (x1 - x0);
        }
    }
    points
        .last()
        .map_or(1.0, |&(last, e)| e * last as f64 / m as f64)
}

/// The predicted span of `units` equal batches on `k` slots of `m = threads / k` threads, in
/// one batch's one-thread time: rounds x contention / (m x E(m)), with `ceil(units / k)` rounds
/// while `units <= 2k`, else `units / k` plus half a round when `k` does not divide `units`, and a
/// contention of `1 + CONTENTION_PER_SLOT x (min(k, CONTENTION_SLOTS) - 1)`. The half round is a
/// tail allowance (the prototype's), not a count of rounds: it makes the span non-monotone in
/// the units (8 units on 3 slots predict 3.17 rounds, 9 predict 3), so the choice can depend on
/// the parity of the batch count; it changes no measured cell's choice.
pub(super) fn predicted_span(
    efficiency_points: &[(usize, f64)],
    threads: usize,
    k: usize,
    units: usize,
) -> f64 {
    let k = k.max(1);
    let units = units.max(1);
    let m = (threads / k).max(1);
    let rounds = if units <= 2 * k {
        units.div_ceil(k) as f64
    } else {
        units as f64 / k as f64 + if units % k == 0 { 0.0 } else { 0.5 }
    };
    let contention = 1.0 + CONTENTION_PER_SLOT * (k.min(CONTENTION_SLOTS) - 1) as f64;
    rounds * contention / (m as f64 * efficiency(efficiency_points, m))
}

/// The two-slot rule of `LOSAT_X_DCTWOSLOTS` (module header): whether `units` batch units on
/// `threads` threads, where the model takes `k` slots of `m` threads (`k > DC_SLOTS`), run as
/// `DC_SLOTS` slots of `threads / DC_SLOTS` threads instead: `DC_MIN_THREADS` threads or more and
/// at most `DC_SLOT_THREADS` threads per slot, an even unit count, and either slots of at most
/// `DC_NARROW_SLOT` threads in the model's layout that leave threads idle which the two slots
/// use (`k x m < DC_SLOTS x (threads / DC_SLOTS)`), or every unit in its first round
/// (`units <= k`), or a last round at most half full within `DC_SHORT_SEARCH_ROUNDS` rounds
/// (`0 < units mod k <= k / 2` and `ceil(units / k)` at most that; the rounds are counted whole
/// here, not with `predicted_span`'s half-round tail allowance).
fn two_slots_win(threads: usize, k: usize, m: usize, units: usize) -> bool {
    let k = k.max(1);
    let half = threads / DC_SLOTS;
    let last = units % k;
    threads >= DC_MIN_THREADS
        && half <= DC_SLOT_THREADS
        && units % DC_SLOTS == 0
        && ((m <= DC_NARROW_SLOT && k.saturating_mul(m) < DC_SLOTS * half)
            || units <= k
            || (last != 0 && 2 * last <= k && units.div_ceil(k) <= DC_SHORT_SEARCH_ROUNDS))
}

/// The eight-slot rule of `LOSAT_X_BNEIGHTSLOTS` (module header): whether, on `threads` threads
/// where the model takes `k` slots, its best layout of at most `BN_SLOTS` slots runs instead:
/// more than `BN_SLOTS` slots, on `BN_MIN_THREADS` to `BN_MAX_THREADS` threads.
fn eight_slots_apply(threads: usize, k: usize) -> bool {
    k > BN_SLOTS && (BN_MIN_THREADS..=BN_MAX_THREADS).contains(&threads)
}

/// A layout of the slot model: `k` slots of `m` threads, and its predicted span over the loop's
/// (`k = 1` on the whole pool).
#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Layout {
    pub k: usize,
    pub m: usize,
    pub ratio: f64,
}

/// No NCBI counterpart (scheduling only): the slots and the threads per slot for `threads`
/// threads, a budget of `budget_bytes` for the batches in flight beyond the loop's one,
/// `estimate_bytes` per batch, `units` batch units of work left (`batch_units`), the task's
/// E(m) points, the slot count of `LOSAT_X_QBATCHPARSLOTS` and the task's rule while its switch is
/// on (`SlotRule`); with the layout the model judged: its best of two slots or more, or the one
/// the tuning variable sets (`None` when the model did not run: too small a budget, or fewer
/// than two units).
///
/// `k` runs over `1 ..= min(1 + budget / estimate, threads, units)` and `m = threads / k`; the
/// least `predicted_span` of `k >= 2` (the fewest slots on a tie) is taken when its ratio to
/// `k = 1`'s pays (`slots_pay`: below 0.970 in thousandths), else the search does not speculate
/// (`Off::Model`). With `SlotRule::TwoSlots` (LOSAT_X_DCTWOSLOTS) a layout of more than
/// `DC_SLOTS` slots runs as `DC_SLOTS` slots of `threads / DC_SLOTS` threads where the two-slot
/// rule holds (`two_slots_win`); with `SlotRule::EightSlots` (LOSAT_X_BNEIGHTSLOTS) a layout of
/// more than `BN_SLOTS` slots runs as the least `predicted_span` of `2 ..= BN_SLOTS` slots (the
/// fewest on a tie) where the eight-slot rule holds (`eight_slots_apply`); either way the layout
/// judged stays the model's. The tuning variable's slots are not subject to the rules.
///
/// Exactness: this chooses only `k` and `m`. Which batches a slot searches, and that the loop
/// takes a slot's result only under the key it computed itself, do not depend on them (module
/// header), so the output is the loop's for every layout, as with the earlier rule.
pub(super) fn choose_slots(
    threads: usize,
    budget_bytes: u64,
    estimate_bytes: u64,
    slots: Option<usize>,
    rule: SlotRule,
    units: usize,
    efficiency_points: &[(usize, f64)],
) -> (Result<(usize, usize), Off>, Option<Layout>) {
    let layout = |k: usize| Layout {
        k,
        m: threads / k,
        ratio: predicted_span(efficiency_points, threads, k, units)
            / predicted_span(efficiency_points, threads, 1, units),
    };
    if let Some(slots) = slots {
        let k = slots.clamp(1, threads.max(1));
        return (Ok((k, threads / k)), Some(layout(k)));
    }
    let by_budget = 1 + budget_bytes / estimate_bytes.max(1);
    if by_budget < 2 {
        return (Err(Off::Budget), None);
    }
    let k_max = usize::try_from(by_budget)
        .unwrap_or(usize::MAX)
        .min(threads)
        .min(units);
    if k_max < 2 {
        return (Err(Off::Single), None);
    }
    let best_of = |k_max: usize| {
        let mut best = layout(2);
        for k in 3..=k_max {
            let candidate = layout(k);
            if candidate.ratio < best.ratio {
                best = candidate;
            }
        }
        best
    };
    let best = best_of(k_max);
    if slots_pay(best.ratio) {
        let k = match rule {
            SlotRule::TwoSlots
                if best.k > DC_SLOTS && two_slots_win(threads, best.k, best.m, units) =>
            {
                DC_SLOTS
            }
            SlotRule::EightSlots if eight_slots_apply(threads, best.k) => best_of(BN_SLOTS).k,
            _ => best.k,
        };
        (Ok((k, threads / k)), Some(best))
    } else {
        (Err(Off::Model), Some(best))
    }
}

/// The discarded batches beyond which the speculation stops, with `k` slots after `taken`
/// batches taken from them: a plan's window of waste, plus a quarter of the work taken.
pub(super) fn waste_bound(k: usize, taken: u64) -> u64 {
    WASTE_SLOTS * k as u64 + taken / WASTE_SHARE
}

/// What the executor needs to know of a slot's search result.
pub(super) trait SlotResult {
    /// The search ended with an error (not a cancellation): the loop stops at this batch if it
    /// takes it, so the batches planned after it are not started.
    fn failed(&self) -> bool;
    /// The search was cancelled (`Cancelled`); such a result is never taken.
    fn cancelled(&self) -> bool;
    /// The memory the result holds until the loop takes it.
    fn held_bytes(&self) -> u64;
}

/// No NCBI counterpart: the error of a slot's search that its cancellation ended (`Cancel`). It
/// is never taken by the loop.
#[derive(Debug)]
pub(super) struct Cancelled;

impl std::fmt::Display for Cancelled {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("a query batch searched ahead of the loop was cancelled")
    }
}

impl std::error::Error for Cancelled {}

/// No NCBI counterpart: the cancellation flag of a slot, attached to the pool its searches get
/// (`SearchPool::x_with_ext`); `search_query_batch` reads it (`cancelled`).
#[derive(Debug, Default)]
pub(super) struct Cancel(AtomicBool);

impl Cancel {
    pub(super) fn is_set(&self) -> bool {
        self.0.load(Ordering::Relaxed)
    }

    fn set(&self) {
        self.0.store(true, Ordering::Relaxed);
    }

    fn clear(&self) {
        self.0.store(false, Ordering::Relaxed);
    }
}

/// Whether the search on `pool` is a slot's whose batch was cancelled (never for the loop's own
/// searches, whose pool has no flag).
pub(super) fn cancelled(pool: &SearchPool<'_>) -> bool {
    pool.x_ext::<Cancel>().is_some_and(Cancel::is_set)
}

/// What a slot made of a batch: its value, or its panic with the messages the panic hook
/// kept instead of printing them.
pub(super) enum Outcome<T> {
    Value(T),
    Panic(Box<dyn std::any::Any + Send>, Vec<String>),
}

impl<T: SlotResult> Outcome<T> {
    /// The value; a panic is printed (its kept messages) and resumed on this thread.
    fn resume(self) -> T {
        match self {
            Outcome::Value(value) => value,
            Outcome::Panic(payload, messages) => {
                for message in &messages {
                    eprint!("{message}");
                }
                std::panic::resume_unwind(payload)
            }
        }
    }

    /// An error or a panic: the loop stops at this batch if it takes it.
    fn failed(&self) -> bool {
        match self {
            Outcome::Value(value) => value.failed(),
            Outcome::Panic(..) => true,
        }
    }
}

/// The batches of the slots.
struct State<T> {
    /// Planned batches that no slot has started, in batch order.
    queue: VecDeque<Job>,
    /// The batch each slot is searching, with its estimate.
    running: Vec<Option<(BatchKey, u64)>>,
    /// Finished batches of the current plan, with the memory their results hold.
    done: HashMap<BatchKey, (Outcome<T>, u64)>,
    /// The keys of the current plan.
    wanted: HashSet<BatchKey>,
    /// The loop's next key and its estimate: always admitted.
    next: Option<(BatchKey, u64)>,
    /// The memory the batches in flight may use beyond the loop's next batch.
    budget: u64,
    shutdown: bool,
    /// Finished or cancelled batches dropped unused, and how many of them were cancelled.
    discarded: u64,
    cancelled: u64,
    /// Slot batches that panicked (kept for the loop or dropped).
    panicked: u64,
    /// The most slots searching at once, and the most memory admitted at once.
    peak_busy: usize,
    peak_bytes: u64,
    /// Slot threads that have not ended.
    alive: usize,
}

impl<T> State<T> {
    fn new() -> Self {
        State {
            queue: VecDeque::new(),
            running: Vec::new(),
            done: HashMap::new(),
            wanted: HashSet::new(),
            next: None,
            budget: 0,
            shutdown: false,
            discarded: 0,
            cancelled: 0,
            panicked: 0,
            peak_busy: 0,
            peak_bytes: 0,
            alive: 0,
        }
    }

    /// The memory of the batches running and of the results not taken yet.
    fn in_use(&self) -> u64 {
        let running: u64 = self.running.iter().flatten().map(|(_, bytes)| bytes).sum();
        let held: u64 = self.done.values().map(|(_, bytes)| bytes).sum();
        running + held
    }

    /// Whether a slot may start `job` now: the loop's next batch always, any other while the
    /// memory in use with it stays within the budget plus the next batch's estimate.
    fn admits(&self, job: &Job) -> bool {
        match self.next {
            Some((key, _)) if key == job.key => true,
            next => self.in_use() + job.bytes <= self.budget + next.map_or(0, |(_, bytes)| bytes),
        }
    }

    fn note_peak(&mut self) {
        self.peak_bytes = self.peak_bytes.max(self.in_use());
        let busy = self.running.iter().flatten().count();
        self.peak_busy = self.peak_busy.max(busy);
    }
}

/// The counters of the summary lines.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Stats {
    pub batches: u64,
    pub precomputed: u64,
    pub loop_searched: u64,
    pub mispredicted: u64,
    pub k: usize,
    pub m: usize,
    pub idle: usize,
    /// The batch units after the first batch (`batch_units`; 0 when no batch follows it).
    pub units: usize,
    /// The layout the slot model judged (`choose_slots`), 0 when it did not run.
    pub rule_k: usize,
    pub rule_m: usize,
    /// Its predicted span over the loop's, in thousandths (`None` when the model did not run).
    pub predicted_milli: Option<u64>,
    pub budget_mib: u64,
    pub estimate_mib: u64,
    pub waste_bound: u64,
    pub off: Option<Off>,
}

/// What the run's counters were at the end, with the executor's own.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(super) struct Summary {
    pub stats: Stats,
    pub discarded: u64,
    pub cancelled: u64,
    pub panicked: u64,
    pub peak_busy: usize,
    pub peak_bytes: u64,
}

/// A slot's kept panic messages.
type PanicSink = Mutex<Vec<String>>;

/// What the loop and the slots share: declared before the thread scope of the loop, and the
/// summary lines at its end.
pub(super) struct Shared<T> {
    switches: Switches,
    state: Mutex<State<T>>,
    wake: Condvar,
    sinks: OnceLock<Vec<Arc<PanicSink>>>,
    cancels: OnceLock<Vec<Cancel>>,
    stats: Mutex<Stats>,
}

impl<T> Shared<T> {
    pub(super) fn new(switches: Switches) -> Self {
        Shared {
            switches,
            state: Mutex::new(State::new()),
            wake: Condvar::new(),
            sinks: OnceLock::new(),
            cancels: OnceLock::new(),
            stats: Mutex::new(Stats {
                off: Some(if switches.on {
                    Off::Unstarted
                } else {
                    Off::NotSet
                }),
                ..Stats::default()
            }),
        }
    }

    pub(super) fn switches(&self) -> Switches {
        self.switches
    }

    fn lock(&self) -> MutexGuard<'_, State<T>> {
        self.state.lock().unwrap_or_else(PoisonError::into_inner)
    }

    /// PROBE (LOSAT_PROBE_BNAHEAD): the loop's thread is thread 0 of the pool the slots search
    /// on; it runs the pool's jobs while it waits for them.
    fn x_wait_loop<'g>(&'g self, guard: MutexGuard<'g, State<T>>) -> MutexGuard<'g, State<T>> {
        if !super::x_probe_ahead::on() {
            return self.wait(guard);
        }
        drop(guard);
        #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
        if rayon::yield_now() == Some(rayon::Yield::Executed) {
            return self.lock();
        }
        let guard = self.lock();
        self.wake
            .wait_timeout(guard, std::time::Duration::from_micros(200))
            .unwrap_or_else(PoisonError::into_inner)
            .0
    }

    fn wait<'g>(&self, guard: MutexGuard<'g, State<T>>) -> MutexGuard<'g, State<T>> {
        self.wake
            .wait(guard)
            .unwrap_or_else(PoisonError::into_inner)
    }

    /// Set the cancellation flag of every slot whose running batch `drop` selects.
    fn cancel_running(&self, state: &State<T>, drop: impl Fn(&BatchKey) -> bool) {
        if let Some(cancels) = self.cancels.get() {
            for (slot, running) in state.running.iter().enumerate() {
                if running.as_ref().is_some_and(|(key, _)| drop(key)) {
                    cancels[slot].set();
                }
            }
        }
    }

    /// The counters so far.
    pub(super) fn summary(&self) -> Summary {
        let stats = *self.stats.lock().unwrap_or_else(PoisonError::into_inner);
        let state = self.lock();
        Summary {
            stats,
            discarded: state.discarded,
            cancelled: state.cancelled,
            panicked: state.panicked,
            peak_busy: state.peak_busy,
            peak_bytes: state.peak_bytes,
        }
    }

    /// `[X_STATS] QBATCHPAR ...`
    pub(super) fn stats_line(&self) -> String {
        let summary = self.summary();
        let stats = summary.stats;
        format!(
            "[X_STATS] QBATCHPAR batches={} precomputed={} loop={} mispredicted={} discarded={} cancelled={} panicked={} waste_bound={} k={} m={} idle={} units={} rule_k={} rule_m={} predicted={} budget_mib={} estimate_mib={} peak_mib={} off={}",
            stats.batches,
            stats.precomputed,
            stats.loop_searched,
            stats.mispredicted,
            summary.discarded,
            summary.cancelled,
            summary.panicked,
            stats.waste_bound,
            stats.k,
            stats.m,
            stats.idle,
            stats.units,
            stats.rule_k,
            stats.rule_m,
            stats
                .predicted_milli
                .map_or_else(|| "-".to_string(), |milli| format!(
                    "{}.{:03}",
                    milli / 1000,
                    milli % 1000
                )),
            stats.budget_mib,
            stats.estimate_mib,
            summary.peak_bytes.div_ceil(MIB),
            stats.off.map_or("none", Off::name)
        )
    }
}

impl<T> Drop for Shared<T> {
    // No NCBI counterpart: the summary line at the end of a search: the counters with
    // LOSAT_X_STATS (as the other counters of the tree).
    fn drop(&mut self) {
        if !self.switches.on || std::thread::panicking() {
            return;
        }
        if self.switches.stats {
            eprintln!("{}", self.stats_line());
        }
    }
}

/// Whether the search speculates.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Phase {
    /// Before the first batch is committed.
    Pending,
    Active,
    Off(Off),
}

/// The loop's side: the planner, the plan, and the slots' threads (started in the thread scope
/// of the loop when the first batch is committed).
pub(super) struct Driver<'scope, 'env, T, F> {
    shared: &'env Shared<T>,
    search: &'env F,
    scope: &'scope std::thread::Scope<'scope, 'env>,
    input: PlanInput<'env>,
    threads: usize,
    profile: Profile,
    /// The lookup width of a batch of these letters and longest query.
    lut_width: Box<dyn Fn(u64, usize) -> usize + 'env>,
    planner: Planner,
    phase: Phase,
    /// The keys of the current plan; the first is the next key the loop takes.
    plan: Vec<BatchKey>,
    window: usize,
    commits: u64,
    last_miss: Option<u64>,
    stats: Stats,
    slots: Vec<std::thread::ScopedJoinHandle<'scope, ()>>,
    finished: bool,
}

impl<'scope, 'env, T, F> Driver<'scope, 'env, T, F>
where
    T: SlotResult + Send + 'env,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync + 'env,
{
    /// The driver of one search; `off` is why the search cannot speculate at all (`None` when
    /// it may, once its first batch is committed).
    #[allow(clippy::too_many_arguments)]
    pub(super) fn new(
        shared: &'env Shared<T>,
        search: &'env F,
        scope: &'scope std::thread::Scope<'scope, 'env>,
        input: PlanInput<'env>,
        threads: usize,
        profile: Profile,
        lut_width: Box<dyn Fn(u64, usize) -> usize + 'env>,
        off: Option<Off>,
    ) -> Self {
        let switches = shared.switches;
        let off = if switches.on { off } else { Some(Off::NotSet) };
        if switches.on && off.is_some() {
            crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", false);
        }
        let phase = match off {
            Some(off) => Phase::Off(off),
            None => Phase::Pending,
        };
        Driver {
            shared,
            search,
            scope,
            input,
            threads,
            profile,
            lut_width,
            planner: Planner::new(switches.miss_every),
            phase,
            plan: Vec::new(),
            window: 0,
            commits: 0,
            last_miss: None,
            stats: Stats {
                off: off.or(Some(Off::Unstarted)),
                ..Stats::default()
            },
            slots: Vec::new(),
            finished: false,
        }
    }

    /// Whether a slot may still search a batch of this search: `false` once the driver is off
    /// (nothing starts any more, and `take` has waited for the slots to be idle), and before the
    /// first commit when no slot can start: without `LOSAT_X_QBATCHPARSLOTS`, `choose_slots`
    /// answers a budget of 0 with `Off::Budget` and one thread with `Off::Budget` or `Off::Single`
    /// at that commit. LOSAT_X_MBLUTKEEP keeps the loop's lookup table between its batches only
    /// then (`lookup::x_mb_lut_keep`).
    pub(super) fn may_speculate(&self) -> bool {
        match self.phase {
            Phase::Off(_) => false,
            Phase::Active => true,
            Phase::Pending => {
                self.shared.switches.slots.is_some() || (self.budget_mib() > 0 && self.threads > 1)
            }
        }
    }

    /// EXPERIMENT (LOSAT_X_BNQDROP): the loop takes its batches' records out of their store, which
    /// no slot can read: a search whose slots could start does not speculate (`Off::Taken`; the rule
    /// of `x_query_drop` acts only where `search_cli` found that none could).
    pub(super) fn x_refuse_taken(&mut self) {
        if self.may_speculate() {
            crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", false);
            self.stop(Off::Taken);
        }
    }

    /// The budget of the batches in flight beyond the loop's (MiB): `LOSAT_X_QBATCHPARMEM`, or the
    /// profile's allowance per thread beyond the first.
    fn budget_mib(&self) -> u64 {
        self.shared
            .switches
            .mem_mib
            .unwrap_or(self.profile.allowance_mib * (self.threads as u64).saturating_sub(1))
    }

    // No NCBI counterpart: tells the report writer whether the search's pool is idle; it does not
    // change any value NCBI computes.
    /// EXPERIMENT (LOSAT_X_PAIRWISEPAR): whether no slot's thread exists (none started, or all
    /// joined by `finish`), so that the search's pool is idle between the loop's searches.
    pub(super) fn x_slots_idle(&self) -> bool {
        self.slots.is_empty()
    }

    /// The result of batch `key` if a slot searched it (waiting for it if a slot is searching
    /// it or will); `None` when the loop searches it itself, on the whole pool: then no slot is
    /// running and none can start. The first batch is never searched ahead.
    pub(super) fn take(&mut self, key: BatchKey) -> Option<T> {
        if key.input_start != 0 && self.phase == Phase::Active {
            let mut state = self.shared.lock();
            let off = loop {
                if let Some((outcome, _)) = state.done.remove(&key) {
                    drop(state);
                    // Its memory is free for the next batches.
                    self.shared.wake.notify_all();
                    self.stats.precomputed += 1;
                    self.flush();
                    return Some(outcome.resume());
                }
                let planned = state
                    .running
                    .iter()
                    .flatten()
                    .any(|(running, _)| *running == key)
                    || state.queue.iter().any(|job| job.key == key);
                if !planned {
                    break Off::Unplanned;
                }
                // Every slot's thread ended outside a batch (after a panic it kept): nothing
                // will search the queued key.
                if state.alive == 0 {
                    break Off::Spawn;
                }
                state = self.shared.x_wait_loop(state);
            };
            drop(state);
            self.stop(off);
        }
        self.wait_idle();
        self.stats.loop_searched += 1;
        self.flush();
        None
    }

    /// Batch `key` is done in the loop with `gie` successful initial extensions; `mixer` is the
    /// loop's mixer after it (`None` with `BATCH_SIZE`) and `next_size` the next batch's size.
    pub(super) fn commit(
        &mut self,
        key: BatchKey,
        gie: u64,
        mixer: Option<&BatchSizeMixer>,
        next_size: u32,
    ) {
        self.stats.batches += 1;
        self.commits += 1;
        let searched = self
            .input
            .query_input
            .searched_in(key.input_start..key.input_end);
        let letters = self.input.lengths[searched]
            .iter()
            .map(|&length| length as u64)
            .sum();
        self.planner.record(key.batch_size, gie, letters);
        if self.phase == Phase::Pending {
            self.activate(key, mixer, next_size);
        }
        if self.phase == Phase::Active {
            let next = self.planner.plan(
                &self.input,
                key.input_end,
                next_size,
                mixer,
                self.window,
                true,
            );
            let predicted = if self.plan.first() == Some(&key) {
                self.plan.get(1).copied()
            } else {
                None
            };
            let mut off = None;
            if let (Some(predicted), Some(actual)) = (predicted, next.first()) {
                if predicted != actual.key {
                    self.stats.mispredicted += 1;
                    let commit = self.commits;
                    if self
                        .last_miss
                        .is_some_and(|last| commit - last < GUARD_COMMITS)
                    {
                        off = Some(Off::Guard);
                    }
                    self.last_miss = Some(commit);
                }
            }
            match off {
                Some(off) => self.stop(off),
                None => self.replan(next),
            }
            // The waste bound: discarded (dropped or cancelled) batches.
            self.stats.waste_bound = waste_bound(self.stats.k, self.stats.precomputed);
            if self.phase == Phase::Active && self.shared.lock().discarded > self.stats.waste_bound
            {
                self.stop(Off::Waste);
            }
        }
        self.flush();
    }

    /// The end of the speculation before the loop writes output (the reports after the loop, or
    /// what an error writes; `LOSAT_X_RPTSTREAM`'s per-batch reports are written while the slots
    /// run): nothing more starts, running batches are cancelled, and every slot's thread is
    /// joined. Also run when the driver is dropped.
    pub(super) fn finish(&mut self) {
        if self.finished {
            return;
        }
        self.finished = true;
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.shutdown = true;
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
        let started = !self.slots.is_empty();
        for slot in self.slots.drain(..) {
            // A slot's thread catches the panics of its batches; a panic outside them (not
            // expected) has ended the thread, and the batches were dropped above.
            let _ = slot.join();
            UNJOINED.with(|unjoined| unjoined.set(unjoined.get().saturating_sub(1)));
        }
        if started && !std::thread::panicking() {
            check_no_global_pool();
        }
        self.flush();
    }

    fn flush(&self) {
        *self
            .shared
            .stats
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = self.stats;
    }

    fn set_off(&mut self, off: Off) {
        self.phase = Phase::Off(off);
        self.stats.off = Some(off);
    }

    /// After the first batch: plan the whole input once, choose the slots, start them.
    fn activate(&mut self, key: BatchKey, mixer: Option<&BatchSizeMixer>, next_size: u32) {
        let full = self.planner.plan(
            &self.input,
            key.input_end,
            next_size,
            mixer,
            usize::MAX,
            false,
        );
        let letters: Vec<u64> = full.iter().map(|planned| planned.letters).collect();
        let units = batch_units(&letters);
        self.stats.units = units;
        if full.len() < 2 {
            crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", false);
            self.set_off(Off::Single);
            return;
        }
        let estimate = full
            .iter()
            .map(|planned| self.estimate(planned))
            .max()
            .unwrap_or(0);
        let switches = self.shared.switches;
        let budget_mib = self.budget_mib();
        self.stats.budget_mib = budget_mib;
        self.stats.estimate_mib = estimate.div_ceil(MIB);
        let rule = match self.profile.rule {
            SlotRule::TwoSlots if switches.dc_two => SlotRule::TwoSlots,
            SlotRule::EightSlots if switches.bn_eight => SlotRule::EightSlots,
            _ => SlotRule::Model,
        };
        let (choice, layout) = choose_slots(
            self.threads,
            budget_mib.saturating_mul(MIB),
            estimate,
            switches.slots,
            rule,
            units,
            self.profile.efficiency,
        );
        // LOSAT_X_DCTWOSLOTS and LOSAT_X_BNEIGHTSLOTS decide where the model itself takes slots
        // for a task with the rule (marked once the slots started): the rule's own side where
        // fewer slots run than the model judged, the reference side where the model's layout runs
        // (as few slots of its own, or the rule declined).
        let capped = match (rule, switches.slots, choice, layout) {
            (SlotRule::Model, ..) => None,
            (rule, None, Ok((k, _)), Some(judged)) => Some((rule, k < judged.k)),
            _ => None,
        };
        if let Some(layout) = layout {
            self.stats.rule_k = layout.k;
            self.stats.rule_m = layout.m;
            self.stats.predicted_milli = Some(milli(layout.ratio));
        }
        let (k, m) = match choice {
            Ok(slots) => slots,
            Err(off) => {
                crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", false);
                self.set_off(off);
                return;
            }
        };
        // PROBE (LOSAT_PROBE_BNAHEAD): one-thread set-up slots.
        let m = if super::x_probe_ahead::on() { 1 } else { m };
        self.stats.k = k;
        self.stats.m = m;
        self.stats.idle = self.threads.saturating_sub(k * m);
        self.shared.lock().budget = budget_mib.saturating_mul(MIB);
        match self.start_slots(k, m) {
            Ok(()) => {
                crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", true);
                match capped {
                    Some((SlotRule::TwoSlots, fewer)) => {
                        crate::utils::xreach::x_reach!("LOSAT_X_DCTWOSLOTS", fewer);
                    }
                    Some((SlotRule::EightSlots, fewer)) => {
                        crate::utils::xreach::x_reach!("LOSAT_X_BNEIGHTSLOTS", fewer);
                    }
                    _ => {}
                }
                self.phase = Phase::Active;
                self.stats.off = None;
                self.window = 2 * k;
            }
            Err(off) => {
                crate::utils::xreach::x_reach!("LOSAT_X_QBATCHPAR", false);
                self.stop(off);
            }
        }
    }

    /// The memory estimate of a planned batch.
    fn estimate(&self, planned: &Planned) -> u64 {
        batch_estimate_bytes(
            (self.lut_width)(planned.letters, planned.max_len),
            planned.letters,
            planned.gie,
            self.profile.query_bytes_per_letter,
        )
    }

    #[cfg(not(all(feature = "parallel", not(target_arch = "wasm32"))))]
    fn start_slots(&mut self, _k: usize, _m: usize) -> Result<(), Off> {
        Err(Off::Build)
    }

    /// `k` slot threads (named `losat-qbp-<slot>`), each with a pool of `m` threads when `m > 1`
    /// (`slot_main`).
    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    fn start_slots(&mut self, k: usize, m: usize) -> Result<(), Off> {
        install_panic_hook();
        let sinks: Vec<Arc<PanicSink>> = (0..k).map(|_| Arc::default()).collect();
        let cancels: Vec<Cancel> = (0..k).map(|_| Cancel::default()).collect();
        if self.shared.sinks.set(sinks).is_err() || self.shared.cancels.set(cancels).is_err() {
            return Err(Off::Spawn);
        }
        self.shared.lock().running = vec![None; k];
        for slot in 0..k {
            let (shared, search) = (self.shared, self.search);
            self.shared.lock().alive += 1;
            let spawned = std::thread::Builder::new()
                .name(format!("losat-qbp-{slot}"))
                .stack_size(SLOT_STACK_BYTES)
                .spawn_scoped(self.scope, move || slot_main(shared, search, slot, m));
            match spawned {
                Ok(handle) => {
                    self.slots.push(handle);
                    UNJOINED.with(|unjoined| unjoined.set(unjoined.get() + 1));
                }
                Err(_) => {
                    self.shared.lock().alive -= 1;
                    return Err(Off::Spawn);
                }
            }
        }
        Ok(())
    }

    /// Stop speculating: nothing more is started, running batches are cancelled, finished
    /// results are dropped.
    fn stop(&mut self, off: Off) {
        self.set_off(off);
        self.plan.clear();
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
    }

    /// The plan `next` replaces the current one: its keys that no slot has started or finished
    /// are queued in order (up to a finished batch that failed: the loop stops there if it
    /// reaches it); finished batches outside it are dropped and running ones cancelled.
    fn replan(&mut self, next: Vec<Planned>) {
        let jobs: Vec<Job> = next
            .iter()
            .map(|planned| Job {
                key: planned.key,
                searched: planned.searched.clone(),
                bytes: self.estimate(planned),
            })
            .collect();
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            let mut jobs = jobs;
            if let Some(failed) = jobs.iter().position(|job| {
                state
                    .done
                    .get(&job.key)
                    .is_some_and(|(outcome, _)| outcome.failed())
            }) {
                jobs.truncate(failed + 1);
            }
            self.plan = jobs.iter().map(|job| job.key).collect();
            let wanted: HashSet<BatchKey> = self.plan.iter().copied().collect();
            let stale: Vec<BatchKey> = state
                .done
                .keys()
                .filter(|key| !wanted.contains(key))
                .copied()
                .collect();
            let unused: Vec<(Outcome<T>, u64)> = stale
                .iter()
                .filter_map(|key| state.done.remove(key))
                .collect();
            state.discarded += unused.len() as u64;
            self.shared
                .cancel_running(&state, |key| !wanted.contains(key));
            state.next = jobs.first().map(|job| (job.key, job.bytes));
            let queue: VecDeque<Job> = jobs
                .into_iter()
                .filter(|job| {
                    !state
                        .running
                        .iter()
                        .flatten()
                        .any(|(running, _)| *running == job.key)
                        && !state.done.contains_key(&job.key)
                })
                .collect();
            state.queue = queue;
            state.wanted = wanted;
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
    }

    /// Wait until no slot is searching (none starts: nothing is queued when this is called).
    fn wait_idle(&self) {
        let mut state = self.shared.lock();
        while state.running.iter().any(Option::is_some) {
            state = self.shared.x_wait_loop(state);
        }
    }
}

impl<T, F> Drop for Driver<'_, '_, T, F> {
    // The end of the loop without `finish` (an early return or a panic): the same steps.
    fn drop(&mut self) {
        if self.finished {
            return;
        }
        self.finished = true;
        let unused: Vec<(Outcome<T>, u64)> = {
            let mut state = self.shared.lock();
            state.shutdown = true;
            state.queue.clear();
            state.wanted.clear();
            state.next = None;
            let unused: Vec<(Outcome<T>, u64)> = state.done.drain().map(|(_, done)| done).collect();
            state.discarded += unused.len() as u64;
            self.shared.cancel_running(&state, |_| true);
            unused
        };
        self.shared.wake.notify_all();
        drop(unused);
        for slot in self.slots.drain(..) {
            let _ = slot.join();
            UNJOINED.with(|unjoined| unjoined.set(unjoined.get().saturating_sub(1)));
        }
        *self
            .shared
            .stats
            .lock()
            .unwrap_or_else(PoisonError::into_inner) = self.stats;
    }
}

/// A slot's thread: search queued batches until the loop ends. With `threads > 1` the slot's
/// thread is the first thread of a pool of `threads` (`with_search_pool`, as the search's own
/// pool has the loop's thread), so a batch's search runs on the slot's thread with that pool
/// around it for its parallel stages (one heap arena keeps a batch's large buffers, as the loop
/// keeps them in its own), and the pool's threads are joined when the slot ends. A pool that
/// cannot be built leaves the slot searching on its own thread. Either pool carries the slot's
/// cancellation flag.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_main<T, F>(shared: &Shared<T>, search: &F, slot: usize, threads: usize)
where
    T: SlotResult + Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    let sink = Arc::clone(
        &shared
            .sinks
            .get()
            .expect("slot sinks are set before the slots")[slot],
    );
    let cancel = &shared
        .cancels
        .get()
        .expect("slot flags are set before the slots")[slot];
    set_panic_sink(Some(Arc::clone(&sink)));
    // When the slot's thread ends (also by a panic outside a batch, not expected), the loop must
    // not wait for it.
    struct Ended<'a, T>(&'a Shared<T>, usize);
    impl<T> Drop for Ended<'_, T> {
        fn drop(&mut self) {
            let mut state = self.0.lock();
            state.alive -= 1;
            if let Some(running) = state.running.get_mut(self.1) {
                *running = None;
            }
            drop(state);
            self.0.wake.notify_all();
        }
    }
    let _ended = Ended(shared, slot);
    if threads > 1 {
        let served = crate::utils::threading::with_search_pool(threads, "blastn", |pool| {
            rayon::broadcast(|_| set_panic_sink(Some(Arc::clone(&sink))));
            slot_loop(
                shared,
                search,
                slot,
                &pool.x_with_ext(cancel),
                &sink,
                cancel,
            );
            Ok(())
        });
        if served.is_ok() {
            return;
        }
    }
    let serial = SearchPool::x_serial();
    slot_loop(
        shared,
        search,
        slot,
        &serial.x_with_ext(cancel),
        &sink,
        cancel,
    );
}

/// The batches of one slot, each searched with `pool`, until the loop ends or a batch panics
/// (its thread-local state is then not trusted further; the panic is kept for the loop).
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn slot_loop<T, F>(
    shared: &Shared<T>,
    search: &F,
    slot: usize,
    pool: &SearchPool<'_>,
    sink: &PanicSink,
    cancel: &Cancel,
) where
    T: SlotResult + Send,
    F: Fn(&Job, &SearchPool<'_>) -> T + Sync,
{
    loop {
        let job = {
            let mut state = shared.lock();
            loop {
                if state.shutdown {
                    return;
                }
                if state.queue.front().is_some_and(|job| state.admits(job)) {
                    let job = state.queue.pop_front().expect("a queued batch");
                    state.running[slot] = Some((job.key, job.bytes));
                    cancel.clear();
                    state.note_peak();
                    break job;
                }
                state = shared.wait(state);
            }
        };
        sink.lock().unwrap_or_else(PoisonError::into_inner).clear();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| search(&job, pool)));
        let panicked = result.is_err();
        let outcome = match result {
            Ok(value) => Outcome::Value(value),
            Err(payload) => Outcome::Panic(
                payload,
                std::mem::take(&mut *sink.lock().unwrap_or_else(PoisonError::into_inner)),
            ),
        };
        // The result's size, before the lock (it walks every HSP).
        let held = match &outcome {
            Outcome::Value(value) => value.held_bytes(),
            Outcome::Panic(..) => 0,
        };
        let unused = {
            let mut state = shared.lock();
            // The flag as it stands when the search has returned (it is set and cleared only
            // under this lock).
            let flagged = cancel.is_set();
            settle(&mut state, slot, job, outcome, held, flagged)
        };
        shared.wake.notify_all();
        drop(unused);
        if panicked {
            return;
        }
    }
}

/// What becomes of a slot's finished batch, under the state lock; returns what to drop outside
/// it. `held` is the result's size and `flagged` whether the slot's cancellation flag was set
/// when the search returned. A flagged search counts as cancelled whatever it returned (a search
/// may stop part-way at any of its checks, and a later edit could return a partial value
/// there), so it is never taken: it is counted, and searched again if the plan wants its key.
/// A panic is counted whether it is kept or dropped.
fn settle<T: SlotResult>(
    state: &mut State<T>,
    slot: usize,
    job: Job,
    outcome: Outcome<T>,
    held: u64,
    flagged: bool,
) -> Option<Outcome<T>> {
    state.running[slot] = None;
    if matches!(&outcome, Outcome::Panic(..)) {
        state.panicked += 1;
    }
    let wanted =
        !state.shutdown && state.wanted.contains(&job.key) && !state.done.contains_key(&job.key);
    let cancelled = flagged || matches!(&outcome, Outcome::Value(value) if value.cancelled());
    if cancelled {
        state.discarded += 1;
        state.cancelled += 1;
        if wanted {
            let at = state
                .queue
                .iter()
                .position(|queued| queued.key.input_start > job.key.input_start)
                .unwrap_or(state.queue.len());
            state.queue.insert(at, job);
        }
        return Some(outcome);
    }
    if !wanted {
        state.discarded += 1;
        return Some(outcome);
    }
    if outcome.failed() {
        // The loop stops at this batch if it reaches it: nothing after it starts.
        state
            .queue
            .retain(|queued| queued.key.input_start < job.key.input_start);
    }
    state.done.insert(job.key, (outcome, held));
    state.note_peak();
    None
}

thread_local! {
    // No NCBI counterpart: the slot threads of the search on this thread (the loop's) that have
    // not been joined; the writers of `run_in_pool` check that there are none.
    static UNJOINED: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// The writers of `run_in_pool` (the reports after the loop, and what a failed batch writes) run
/// only after `Driver::finish` has joined every slot of the search, so that no speculative batch
/// runs while they write (`LOSAT_X_RPTSTREAM`'s per-batch report writes, made while the slots
/// search ahead, do not call this check). With debug assertions this checks it, on the loop's
/// thread.
pub(super) fn debug_assert_slots_joined() {
    debug_assert_eq!(
        UNJOINED.with(std::cell::Cell::get),
        0,
        "a writer of run_in_pool runs while LOSAT_X_QBATCHPAR slots may still run (Driver::finish was not called before it)"
    );
}

thread_local! {
    // No NCBI counterpart: where the panic hook keeps the message of a panic on a slot's
    // threads (diagnostics only).
    static PANIC_SINK: std::cell::RefCell<Option<Arc<PanicSink>>> =
        const { std::cell::RefCell::new(None) };
}

fn set_panic_sink(sink: Option<Arc<PanicSink>>) {
    PANIC_SINK.with(|current| *current.borrow_mut() = sink);
}

/// The panic hook of the process, once: on a slot's threads it keeps the message (printed if
/// the loop takes the batch); elsewhere it calls the hook it replaced. Not with `panic =
/// "abort"` (the release profile): a panic cannot be kept there, and the process ends with the
/// default hook's message.
#[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
fn install_panic_hook() {
    if cfg!(panic = "abort") {
        return;
    }
    static ONCE: std::sync::Once = std::sync::Once::new();
    ONCE.call_once(|| {
        let previous = std::panic::take_hook();
        std::panic::set_hook(Box::new(move |info| {
            let sink = PANIC_SINK
                .try_with(|current| current.borrow().clone())
                .ok()
                .flatten();
            let Some(sink) = sink else {
                previous(info);
                return;
            };
            let payload = info
                .payload()
                .downcast_ref::<&str>()
                .map(|text| (*text).to_string())
                .or_else(|| info.payload().downcast_ref::<String>().cloned())
                .unwrap_or_else(|| "Box<dyn Any>".to_string());
            let location = info
                .location()
                .map(|location| {
                    format!(
                        "{}:{}:{}",
                        location.file(),
                        location.line(),
                        location.column()
                    )
                })
                .unwrap_or_default();
            let thread = std::thread::current();
            let message = format!(
                "thread '{}' panicked at {location}:\n{payload}\n",
                thread.name().unwrap_or("<unnamed>")
            );
            sink.lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(message);
        }));
    });
}

/// Design section 7.4: no slot reaches rayon's global pool. A slot of one thread runs on a
/// plain thread, where a `par_iter` or `rayon::join` without a check of
/// `rayon::current_thread_index()` would start the global pool with every CPU; the slots of
/// `m > 1` threads have their own pool. With debug assertions (and not in this crate's unit
/// tests, whose other tests may use the global pool), after the slots are joined: the global
/// pool must not exist yet, unless this check made it (a one-thread pool, so a later search's
/// slots are checked as well as the first one's are only if nothing else started it).
/// Caveats, all in debug builds only: `build_global` can also fail for another reason (a thread
/// that cannot be spawned), and a host that embeds LOSAT and has built the global pool itself
/// trips the check; the check leaves a one-thread global pool in the process. Fine for the
/// command line, whose blastn path never builds the global pool.
fn check_no_global_pool() {
    #[cfg(all(
        debug_assertions,
        not(test),
        feature = "parallel",
        not(target_arch = "wasm32")
    ))]
    {
        static OURS: AtomicBool = AtomicBool::new(false);
        let built = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .thread_name(|_| "losat-qbp-global-check".to_string())
            .build_global()
            .is_ok();
        if built {
            OURS.store(true, Ordering::Relaxed);
        }
        debug_assert!(
            built || OURS.load(Ordering::Relaxed),
            "a LOSAT_X_QBATCHPAR slot started rayon's global pool (a parallel iterator or join without a check of rayon::current_thread_index())"
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blastinput::seq_range::Placements;

    fn input_of(lengths: &[usize], skipped: &[bool]) -> QueryInput {
        let ordinals = (0..lengths.len())
            .filter(|&index| !skipped[index])
            .collect();
        QueryInput {
            placements: Placements::default(),
            ordinals,
            input_lengths: lengths.to_vec(),
            skipped: skipped.to_vec(),
            range_length: None,
        }
    }

    fn searched_lengths(input: &QueryInput) -> Vec<usize> {
        input
            .ordinals
            .iter()
            .map(|&ordinal| input.input_lengths[ordinal])
            .collect()
    }

    fn plan_input<'a>(
        query_input: &'a QueryInput,
        lengths: &'a [usize],
        fixed: bool,
        megablast: bool,
    ) -> PlanInput<'a> {
        PlanInput {
            query_input,
            lengths,
            reading_error: false,
            fixed_batch_size: fixed,
            split: SplitSizes::default_for(megablast),
        }
    }

    /// One batch of the loop: its key, letters and gie, the mixer after it, and the next size.
    type TrueBatch = (BatchKey, u64, u64, BatchSizeMixer, u32);

    /// The loop's batches for `gie` (a function of a batch's searched letters), as `run_in_pool`
    /// reads them.
    fn true_batches(
        input: &PlanInput<'_>,
        first: u32,
        max: i32,
        gie: impl Fn(u64) -> u64,
    ) -> Vec<TrueBatch> {
        let mut mixer = BatchSizeMixer::new(max);
        assert_eq!(mixer.batch_size(None) as u32, first);
        let mut size = first;
        let mut start = 0;
        let mut out = Vec::new();
        let records = input.query_input.input_lengths.len();
        while start < records {
            let (end, _) = next_ranged_batch(
                &input.query_input.input_lengths,
                &input.query_input.skipped,
                start,
                size,
            );
            let searched = input.query_input.searched_in(start..end);
            let letters: u64 = input.lengths[searched].iter().map(|&l| l as u64).sum();
            let extensions = gie(letters);
            let key = BatchKey {
                input_start: start,
                input_end: end,
                batch_size: size,
            };
            size = mixer.batch_size(Some(extensions as u32 as i32)) as u32;
            out.push((key, letters, extensions, mixer.clone(), size));
            start = end;
        }
        out
    }

    /// Replays the loop with the planner: how many next keys the plan mispredicted.
    fn mispredictions(input: &PlanInput<'_>, batches: &[TrueBatch]) -> usize {
        let mut planner = Planner::new(None);
        let mut plan: Vec<BatchKey> = Vec::new();
        let mut misses = 0;
        for (index, (key, letters, gie, mixer, next)) in batches.iter().enumerate() {
            planner.record(key.batch_size, *gie, *letters);
            let next_plan = planner.plan(input, key.input_end, *next, Some(mixer), 8, true);
            if index > 0 && plan.first() == Some(key) {
                if let (Some(predicted), Some(actual)) = (plan.get(1), next_plan.first()) {
                    if predicted != &actual.key {
                        misses += 1;
                    }
                }
            }
            plan = next_plan.iter().map(|planned| planned.key).collect();
            // The plan's first key is the loop's next key.
            if let Some((next_key, ..)) = batches.get(index + 1) {
                assert_eq!(plan.first(), Some(next_key), "batch {index}");
            }
        }
        misses
    }

    // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.cpp:65-81 (`GetBatchSize`).
    // dc-megablast on the audit input: a batch with few extensions per letter (the first one of
    // 5,000, or a full one of 4,999,000 with fewer than ~2,328) makes `(Int4)` of the next size
    // overflow to INT_MIN, so the next batch is the smallest (100, one record), and a small
    // batch is followed by a full one. The scaled rate of the first batch predicts it.
    #[test]
    fn the_planner_follows_the_mixer_alternation_of_tiny_and_full_batches() {
        let lengths: Vec<usize> = (0..12_000)
            .map(|index| 1_500 + (index * 37) % 2_000)
            .collect();
        let query_input = input_of(&lengths, &vec![false; lengths.len()]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, false, true);
        let max = SplitSizes::default_for(true).mixer_max_batch_size();
        let batches = true_batches(&input, 5_000, max, |letters| letters / 4_000);
        let sizes: Vec<u32> = batches.iter().map(|(key, ..)| key.batch_size).collect();
        assert_eq!(&sizes[..5], &[5_000, 100, 4_999_000, 100, 4_999_000]);
        assert!(sizes.len() > 5, "{sizes:?}");
        assert_eq!(mispredictions(&input, &batches), 0);
    }

    // blastn (word 11) on the audit input: every batch after the first has the maximum size
    // (999,000) and its gie stays in the clamped regime.
    #[test]
    fn the_planner_predicts_blastn_batches_without_misses() {
        let lengths: Vec<usize> = (0..3_000)
            .map(|index| 2_000 + (index * 53) % 3_000)
            .collect();
        let query_input = input_of(&lengths, &vec![false; lengths.len()]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, false, false);
        let max = SplitSizes::default_for(false).mixer_max_batch_size();
        let batches = true_batches(&input, 5_000, max, |letters| letters * 72 / 1_000);
        assert!(batches.len() > 6);
        assert!(batches[1..batches.len() - 1]
            .iter()
            .all(|(key, ..)| key.batch_size == 999_000));
        assert_eq!(mispredictions(&input, &batches), 0);
    }

    // `BATCH_SIZE`: every batch has that size, whatever its gie; no prediction can be wrong.
    #[test]
    fn a_fixed_batch_size_plans_exact_keys() {
        let lengths = [300, 300, 300, 300, 300, 300, 300];
        let query_input = input_of(&lengths, &[false; 7]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        planner.record(500, 10, 600);
        let plan = planner.plan(&input, 2, 500, None, 10, true);
        let keys: Vec<(usize, usize)> = plan
            .iter()
            .map(|planned| (planned.key.input_start, planned.key.input_end))
            .collect();
        assert_eq!(keys, [(2, 4), (4, 6), (6, 7)]);
        assert!(plan.iter().all(|planned| planned.key.batch_size == 500));
        // The window bounds the plan.
        assert_eq!(planner.plan(&input, 2, 500, None, 2, true).len(), 2);
    }

    // The plan stops before a batch that ends the run in the loop. A batch reads records until
    // their letters reach its size (skipped records and records without letters add none), so
    // such a batch is the last one.
    #[test]
    fn the_plan_stops_before_a_batch_that_ends_the_run() {
        // Records 3 and 4 are skipped (`-query_loc` past their end): the batch from record 3
        // has no searched query (`Empty CBlastQueryVector`).
        let lengths = [100, 100, 100, 50, 50];
        let skipped = [false, false, false, true, true];
        let query_input = input_of(&lengths, &skipped);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan.len(), 2);
        assert_eq!(plan[1].key.input_end, 3);
        // A skipped record inside the input joins the next batch.
        let skipped = [false, true, false, false, false];
        let query_input = input_of(&lengths, &skipped);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan[0].key.input_end, 3);
        assert_eq!(plan[0].searched, 1..2);
        // Queries without letters only (`Sequence contains no data`).
        let lengths = [100, 100, 0, 0];
        let query_input = input_of(&lengths, &[false; 4]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let plan = planner.plan(&input, 1, 1, None, 10, true);
        assert_eq!(plan.len(), 1);
        // A batch with a query with letters is searched.
        let plan = planner.plan(&input, 1, 150, None, 10, true);
        assert_eq!(plan[0].key.input_end, 4);
        // The reader error after the last record comes in the batch that does not reach its
        // size; a batch that reaches it is searched, and the error comes after it.
        let lengths = [100, 100, 100];
        let query_input = input_of(&lengths, &[false; 3]);
        let searched = searched_lengths(&query_input);
        let mut input = plan_input(&query_input, &searched, true, false);
        input.reading_error = true;
        let plan = planner.plan(&input, 0, 150, None, 10, true);
        assert_eq!(plan.len(), 1, "{plan:?}");
        let plan = planner.plan(&input, 0, 100, None, 10, true);
        assert_eq!(plan.len(), 3, "{plan:?}");
    }

    // A batch that NCBI splits into query chunks reports no extensions of its own.
    #[test]
    fn a_split_batch_is_predicted_without_extensions() {
        let lengths = [2_500_000, 100, 100];
        let query_input = input_of(&lengths, &[false; 3]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(None);
        planner.record(1, 1_000, 1_000);
        let plan = planner.plan(&input, 0, 1, None, 3, true);
        assert_eq!(plan[0].gie, 0);
        assert_eq!(plan[1].gie, 100);
    }

    // `LOSAT_X_QBATCHPARMISS=N`: every N-th prediction gives another size.
    #[test]
    fn the_test_switch_corrupts_every_nth_prediction() {
        let lengths = [100; 12];
        let query_input = input_of(&lengths, &[false; 12]);
        let searched = searched_lengths(&query_input);
        let input = plan_input(&query_input, &searched, true, false);
        let mut planner = Planner::new(Some(2));
        let plan = planner.plan(&input, 0, 200, None, 4, true);
        let sizes: Vec<u32> = plan.iter().map(|planned| planned.key.batch_size).collect();
        assert_eq!(sizes, [200, 200, 100, 100]);
        // Without counting (the whole-input plan of `activate`) nothing is corrupted.
        let plan = planner.plan(&input, 0, 200, None, 4, false);
        assert!(plan.iter().all(|planned| planned.key.batch_size == 200));
    }

    // The estimate is calibrated to the measured memory of one batch in flight: ~50 MiB for a
    // 1 Mb blastn batch (lookup width 11, ~72,000 extensions) and 65-81 MiB per extra 5 Mb
    // dc-megablast slot (width 11, ~1,000 extensions; perf-mb split3 section 3).
    #[test]
    fn the_estimate_matches_the_measured_batches() {
        let blastn = Profile::of_task("blastn").query_bytes_per_letter;
        let dc = Profile::of_task("dc-megablast").query_bytes_per_letter;
        let blastn = batch_estimate_bytes(11, 999_193, 71_763, blastn) as f64 / MIB as f64;
        let dc = batch_estimate_bytes(11, 5_001_442, 1_000, dc) as f64 / MIB as f64;
        assert!((47.0..=54.0).contains(&blastn), "{blastn}");
        assert!((65.0..=82.0).contains(&dc), "{dc}");
        // dc-megablast differs from megablast only by the two-slot rule (LOSAT_X_DCTWOSLOTS), and
        // blastn from blastn-short only by the eight-slot rule (LOSAT_X_BNEIGHTSLOTS).
        assert_eq!(
            Profile::of_task("megablast"),
            Profile {
                rule: SlotRule::Model,
                ..Profile::of_task("dc-megablast")
            }
        );
        assert_eq!(Profile::of_task("dc-megablast").rule, SlotRule::TwoSlots);
        assert_eq!(Profile::of_task("blastn").rule, SlotRule::EightSlots);
        for task in ["megablast", "blastn-short"] {
            assert_eq!(Profile::of_task(task).rule, SlotRule::Model, "{task}");
        }
        assert_eq!(
            Profile::of_task("blastn-short"),
            Profile {
                rule: SlotRule::Model,
                ..Profile::of_task("blastn")
            }
        );
        assert_eq!(allowance_mib("blastn"), 96);
        assert_eq!(allowance_mib("blastn-short"), 96);
        assert_eq!(allowance_mib("dc-megablast"), 28);
        assert_eq!(allowance_mib("megablast"), 28);
    }

    /// NCBI's split by queries (blastn_app.cpp:96-153): a query file larger than 10,000,000
    /// bytes at a word size of 10 or more (1,500,000 below 10), against fewer database letters
    /// than 6e9 (1.5e8); standard input is not split. Elsewhere the default budget is 0.
    #[test]
    fn the_default_budget_follows_ncbis_split_by_queries() {
        assert_eq!(ncbi_mt1_cutoffs(11), (6_000_000_000, 10_000_000));
        assert_eq!(ncbi_mt1_cutoffs(10), (6_000_000_000, 10_000_000));
        assert_eq!(ncbi_mt1_cutoffs(9), (150_000_000, 1_500_000));
        assert_eq!(ncbi_mt1_cutoffs(7), (150_000_000, 1_500_000));
        // The cells of the 2026-10-17 measurements: 1,775 cDNAs (5,398,904 bytes), 3,492
        // (10,779,147) and 31,168 (96,191,719) against 3 Mb.
        assert!(!ncbi_splits_by_queries(11, Some(5_398_904), 3_000_000));
        assert!(ncbi_splits_by_queries(11, Some(10_779_147), 3_000_000));
        assert!(ncbi_splits_by_queries(11, Some(96_191_719), 3_000_000));
        assert!(ncbi_splits_by_queries(28, Some(96_191_719), 3_000_000));
        // The bounds are strict, as NCBI's comparisons.
        assert!(!ncbi_splits_by_queries(11, Some(10_000_000), 0));
        assert!(ncbi_splits_by_queries(11, Some(10_000_001), 0));
        assert!(!ncbi_splits_by_queries(11, Some(10_000_001), 6_000_000_000));
        assert!(ncbi_splits_by_queries(11, Some(10_000_001), 5_999_999_999));
        assert!(!ncbi_splits_by_queries(7, Some(1_500_000), 0));
        assert!(ncbi_splits_by_queries(7, Some(1_500_001), 149_999_999));
        assert!(!ncbi_splits_by_queries(7, Some(1_500_001), 150_000_000));
        assert!(!ncbi_splits_by_queries(11, None, 0));
        for task in ["blastn", "blastn-short", "megablast", "dc-megablast"] {
            // blastn-short's word size 7 has the smaller cut-off.
            let (word, small_file) = if task == "blastn-short" {
                (7, 1_000_000)
            } else {
                (11, 5_398_904)
            };
            let split = Profile::of_search(task, word, Some(96_191_719), 3_000_000);
            assert_eq!(split, Profile::of_task(task), "{task}");
            let small = Profile::of_search(task, word, Some(small_file), 3_000_000);
            assert_eq!(small.allowance_mib, 0, "{task}");
            assert_eq!(
                Profile {
                    allowance_mib: allowance_mib(task),
                    ..small
                },
                Profile::of_task(task),
                "{task}"
            );
            let stdin = Profile::of_search(task, word, None, 3_000_000);
            assert_eq!(stdin.allowance_mib, 0, "{task}");
        }
        // No default budget: the slots do not start (a budget of fewer than two slots).
        let (choice, layout) =
            choose_slots(16, 0, 50 * MIB, None, SlotRule::Model, 5, EFFICIENCY_BLASTN);
        assert_eq!((choice, layout), (Err(Off::Budget), None));
    }

    #[test]
    fn the_efficiency_is_log2_linear_between_the_points() {
        for points in [EFFICIENCY_MEGABLAST, EFFICIENCY_BLASTN] {
            for &(m, e) in points {
                assert!((efficiency(points, m) - e).abs() < 1e-12, "{m}");
            }
            // 3 threads lie between the points of 2 and 4, log2-linear.
            let (e2, e4) = (points[1].1, points[2].1);
            let three = e2 + (e4 - e2) * (3f64.log2() - 1.0);
            assert!((efficiency(points, 3) - three).abs() < 1e-12);
            // Beyond the last point the speedup stays that of 16 threads, and 0 threads count as
            // 1.
            for m in [17, 24, 32, 64, 128] {
                let speedup = m as f64 * efficiency(points, m);
                assert!((speedup - 16.0 * points[4].1).abs() < 1e-9, "{m}");
            }
            assert_eq!(efficiency(points, 0), 1.0);
        }
    }

    #[test]
    fn units_are_the_batches_of_at_least_half_the_median_batch() {
        assert_eq!(batch_units(&[]), 0);
        assert_eq!(reference_letters(&[]), 0);
        assert_eq!(batch_units(&[100]), 1);
        assert_eq!(batch_units(&[597, 5_001_442]), 1);
        assert_eq!(batch_units(&[597, 5_001_442, 152, 4_999_591]), 2);
        // Exactly half counts; one letter less does not.
        assert_eq!(batch_units(&[1_000, 500, 499]), 2);
        // All tiny batches are all units.
        assert_eq!(batch_units(&[100; 40]), 40);
        // The reference is the letter-weighted median: the batch at which the sorted batches
        // reach half of all letters.
        assert_eq!(reference_letters(&[1, 2, 3, 4]), 3);
        assert_eq!(reference_letters(&[10, 1, 1, 1]), 10);
        // One long record (an oversized batch) moves neither the reference nor the other
        // batches' units, and counts as one unit itself.
        let mut outlier = SPLIT3_AB.to_vec();
        outlier.push(3_000_000);
        assert_eq!(batch_units(SPLIT3_AB), 89);
        assert_eq!(batch_units(&outlier), 90);
        let mut contigs = vec![5_000_000u64; 10];
        contigs.push(11_000_000);
        assert_eq!(batch_units(&contigs), 11);
        // The split3 inputs keep their units (as against the largest batch).
        for (planned, units) in [(SPLIT3_A, 1), (SPLIT3_B, 2), (SPLIT3_C, 5), (SPLIT3_AD, 18)] {
            assert_eq!(batch_units(planned), units, "{planned:?}");
        }
    }

    /// The planned batches after the first of the perf-mb split3 inputs (window 2, from the
    /// profiler's `batch_entry` events; D. melanogaster cDNAs against 3 Mb of chromosome 2L).
    /// (a): 1,775 cDNAs, dc-megablast.
    const SPLIT3_A: &[u64] = &[597, 5_001_442];
    /// (b): 3,492 cDNAs, dc-megablast.
    const SPLIT3_B: &[u64] = &[597, 5_001_442, 152, 4_999_591];
    /// (c): 1,775 cDNAs, blastn.
    const SPLIT3_C: &[u64] = &[999_193, 1_004_184, 1_007_311, 1_000_099, 991_252];
    /// Ad: 31,168 cDNAs, dc-megablast (the mixer alternates tiny and 5 Mb batches).
    const SPLIT3_AD: &[u64] = &[
        597, 5_001_442, 152, 4_999_591, 520, 5_003_054, 7_699, 5_001_159, 2_712, 5_000_068, 2_886,
        5_001_250, 4_017, 5_002_579, 3_440, 5_002_671, 7_569, 4_999_636, 6_903, 5_006_181, 7_973,
        4_999_035, 1_413, 5_002_970, 4_063, 5_000_666, 1_484, 5_000_251, 3_399, 5_002_564, 3_100,
        5_000_456, 8_889, 5_003_064, 7_149, 4_179_903,
    ];
    /// Ab: 31,168 cDNAs, blastn (in the order the slots started them).
    const SPLIT3_AB: &[u64] = &[
        1_007_311, 1_004_184, 999_193, 1_000_353, 1_000_099, 1_001_407, 1_003_745, 1_003_157,
        1_000_299, 999_894, 1_001_925, 1_001_163, 999_535, 1_006_649, 1_000_443, 999_611,
        1_000_658, 1_001_795, 1_007_392, 999_866, 999_212, 1_001_313, 1_001_504, 999_303,
        1_003_329, 999_017, 999_338, 999_798, 999_331, 999_279, 1_008_671, 1_003_640, 1_000_849,
        1_002_476, 1_001_822, 1_005_414, 1_002_084, 1_000_502, 1_000_078, 1_011_910, 1_000_776,
        1_002_900, 1_001_242, 1_003_735, 999_091, 1_000_014, 1_001_122, 1_001_639, 1_006_094,
        1_000_240, 999_462, 1_000_646, 1_000_228, 1_002_641, 999_660, 1_005_343, 999_254, 999_243,
        1_002_195, 999_687, 1_003_349, 1_001_580, 1_001_965, 999_853, 999_729, 999_390, 999_394,
        1_000_589, 1_001_936, 999_675, 1_004_210, 1_006_364, 1_003_179, 1_044_868, 999_748,
        1_007_842, 1_001_043, 1_002_396, 999_901, 1_001_106, 999_503, 999_756, 1_000_455,
        1_000_853, 1_000_622, 1_000_009, 1_002_680, 1_000_573, 999_257, 95_919,
    ];

    /// The choice for `planned` batches of `task` at `threads` threads with the default budget
    /// (the estimate of the largest batch: lookup width 11, ~72,000 extensions per 1 Mb blastn
    /// batch, ~1,000 per 5 Mb dc-megablast batch), without the two-slot rule
    /// (LOSAT_X_DCTWOSLOTS=0).
    fn split3_choice(
        task: &str,
        planned: &[u64],
        threads: usize,
    ) -> (Result<(usize, usize), Off>, Option<Layout>, usize) {
        split3_choice_ruled(task, planned, threads, SlotRule::Model)
    }

    /// `split3_choice` with the rule `rule`.
    fn split3_choice_ruled(
        task: &str,
        planned: &[u64],
        threads: usize,
        rule: SlotRule,
    ) -> (Result<(usize, usize), Off>, Option<Layout>, usize) {
        let profile = Profile::of_task(task);
        let largest = planned.iter().copied().max().unwrap();
        let gie = if task == "blastn" {
            largest * 72 / 1_000
        } else {
            1_000
        };
        let estimate = batch_estimate_bytes(11, largest, gie, profile.query_bytes_per_letter);
        let budget = profile.allowance_mib * (threads as u64 - 1) * MIB;
        let units = batch_units(planned);
        let (choice, layout) = choose_slots(
            threads,
            budget,
            estimate,
            None,
            rule,
            units,
            profile.efficiency,
        );
        (choice, layout, units)
    }

    // Every cell of perf-mb `logs/split3_report.md` section 3: the rule gives the layout that
    // was measured best there (window 1 with `LOSAT_X_QBATCHPARSLOTS`, window 2 with the
    // prototype), where the earlier rule (units = planned letters over the largest batch, k by
    // the budget) gave (a) 16 threads 2 x 8, (c) 4 threads 3 x 1, (c) 8 threads 5 x 1, (b) and
    // Ad 4 threads none (budget), Ab 16 threads 15 x 1.
    #[test]
    fn cells_of_the_split3_report() {
        // (input, task, threads, units, expected choice, predicted ratio of the judged layout)
        let cells: [(
            &str,
            &str,
            &[u64],
            usize,
            usize,
            Result<(usize, usize), Off>,
            f64,
        ); 11] = [
            (
                "a",
                "dc-megablast",
                SPLIT3_A,
                4,
                1,
                Err(Off::Single),
                f64::NAN,
            ),
            (
                "a",
                "dc-megablast",
                SPLIT3_A,
                16,
                1,
                Err(Off::Single),
                f64::NAN,
            ),
            ("b", "dc-megablast", SPLIT3_B, 4, 2, Ok((2, 2)), 0.9556),
            (
                "b",
                "dc-megablast",
                SPLIT3_B,
                16,
                2,
                Err(Off::Model),
                1.0088,
            ),
            ("c", "blastn", SPLIT3_C, 4, 5, Err(Off::Model), 1.0525),
            ("c", "blastn", SPLIT3_C, 8, 5, Ok((2, 4)), 0.8873),
            ("c", "blastn", SPLIT3_C, 16, 5, Ok((5, 3)), 0.4833),
            ("Ad", "dc-megablast", SPLIT3_AD, 4, 18, Ok((2, 2)), 0.9556),
            ("Ad", "dc-megablast", SPLIT3_AD, 16, 18, Ok((4, 4)), 0.7730),
            ("Ab", "blastn", SPLIT3_AB, 4, 89, Ok((4, 1)), 0.8348),
            ("Ab", "blastn", SPLIT3_AB, 16, 89, Ok((16, 1)), 0.3793),
        ];
        for (input, task, planned, threads, units, expected, ratio) in cells {
            let (choice, layout, got_units) = split3_choice(task, planned, threads);
            let cell = format!("({input}) {threads} threads");
            assert_eq!(got_units, units, "{cell}");
            assert_eq!(choice, expected, "{cell}: {layout:?}");
            match layout {
                None => assert!(ratio.is_nan(), "{cell}"),
                Some(layout) => {
                    assert!((layout.ratio - ratio).abs() < 1e-4, "{cell}: {layout:?}");
                    assert_eq!(slots_pay(layout.ratio), choice.is_ok(), "{cell}");
                }
            }
        }
        // (b) 16 threads: the best layout of two slots or more is 2 x 8, 1 % slower than the loop.
        let (_, layout, _) = split3_choice("dc-megablast", SPLIT3_B, 16);
        assert_eq!(layout.map(|layout| (layout.k, layout.m)), Some((2, 8)));
        // At one thread the default budget (per thread beyond the first) is 0: no second batch.
        for (task, planned) in [("blastn", SPLIT3_AB), ("dc-megablast", SPLIT3_AD)] {
            assert_eq!(
                split3_choice(task, planned, 1).0,
                Err(Off::Budget),
                "{task}"
            );
        }
    }

    // More threads than were measured (E(m) up to 16 threads per batch, the speedup kept flat
    // beyond): the cells keep slots. With E flat instead, (b) at 24-64 threads, (c) at 64 / 128
    // and Ad at 128 would keep the loop (review S1), where the earlier rule ran slots ((b) 2 x 12
    // / 2 x 16 / 2 x 32 / 2 x 64, (c) 5 x 4 / 5 x 6 / 5 x 12 / 5 x 25, Ad 18 x 7 at 128).
    #[test]
    fn cells_beyond_sixteen_threads() {
        // (input, task, planned, threads, expected choice, predicted ratio)
        let cells: [(&str, &str, &[u64], usize, (usize, usize), f64); 12] = [
            ("b", "dc-megablast", SPLIT3_B, 24, (2, 12), 0.6883),
            ("b", "dc-megablast", SPLIT3_B, 32, (2, 16), 0.5250),
            ("b", "dc-megablast", SPLIT3_B, 64, (2, 32), 0.5250),
            ("b", "dc-megablast", SPLIT3_B, 128, (2, 64), 0.5250),
            ("c", "blastn", SPLIT3_C, 24, (5, 4), 0.3921),
            ("c", "blastn", SPLIT3_C, 32, (5, 6), 0.3161),
            ("c", "blastn", SPLIT3_C, 64, (5, 12), 0.2461),
            ("c", "blastn", SPLIT3_C, 128, (5, 25), 0.2400),
            ("Ad", "dc-megablast", SPLIT3_AD, 24, (6, 4), 0.4840),
            ("Ad", "dc-megablast", SPLIT3_AD, 32, (9, 3), 0.4132),
            ("Ad", "dc-megablast", SPLIT3_AD, 64, (18, 3), 0.2066),
            ("Ad", "dc-megablast", SPLIT3_AD, 128, (18, 7), 0.1315),
        ];
        for (input, task, planned, threads, expected, ratio) in cells {
            let (choice, layout, _) = split3_choice(task, planned, threads);
            let cell = format!("({input}) {threads} threads");
            assert_eq!(choice, Ok(expected), "{cell}: {layout:?}");
            let layout = layout.expect("a judged layout");
            assert!((layout.ratio - ratio).abs() < 1e-4, "{cell}: {layout:?}");
        }
    }

    // LOSAT_X_DCTWOSLOTS on the split3 cells: dc-megablast's rule turns the model's layouts of
    // three slots or more into two slots of half the threads where it holds (Ad, 18 units: 4 x 4
    // at 16 threads, a half-empty fifth round, and 9 x 3 at 32 and 33 threads, narrow slots on 27
    // threads) and
    // leaves every other choice as it was ((a) one unit, (b) two units: the loop at 16 threads,
    // 2 x 2 and 2 x 16 where the model takes two; Ad at 12 and 24 threads: 6 and 3 full rounds of
    // 4-thread slots; at 64 and 128 threads: 32 threads per slot; Ab's 89 units, odd, at every
    // thread count), with the model's judged layout in both cases.
    #[test]
    fn the_dc_rule_takes_two_slots_on_the_split3_cells_where_it_holds() {
        // (input, planned, threads, choice with the rule)
        let cells: [(&str, &[u64], usize, Result<(usize, usize), Off>); 17] = [
            ("a", SPLIT3_A, 16, Err(Off::Single)),
            ("b", SPLIT3_B, 4, Ok((2, 2))),
            ("b", SPLIT3_B, 16, Err(Off::Model)),
            ("b", SPLIT3_B, 32, Ok((2, 16))),
            ("Ad", SPLIT3_AD, 1, Err(Off::Budget)),
            ("Ad", SPLIT3_AD, 4, Ok((2, 2))),
            ("Ad", SPLIT3_AD, 8, Ok((2, 4))),
            ("Ad", SPLIT3_AD, 12, Ok((3, 4))),
            ("Ad", SPLIT3_AD, 16, Ok((2, 8))),
            ("Ad", SPLIT3_AD, 24, Ok((6, 4))),
            ("Ad", SPLIT3_AD, 32, Ok((2, 16))),
            ("Ad", SPLIT3_AD, 33, Ok((2, 16))),
            ("Ad", SPLIT3_AD, 64, Ok((18, 3))),
            ("Ad", SPLIT3_AD, 128, Ok((18, 7))),
            ("Ab", SPLIT3_AB, 16, Ok((15, 1))),
            ("Ab", SPLIT3_AB, 32, Ok((16, 2))),
            ("Ab", SPLIT3_AB, 33, Ok((31, 1))),
        ];
        let rule = Profile::of_task("dc-megablast").rule;
        for (input, planned, threads, expected) in cells {
            let cell = format!("({input}) {threads} threads");
            let (choice, layout, units) =
                split3_choice_ruled("dc-megablast", planned, threads, rule);
            let (model, judged, model_units) = split3_choice("dc-megablast", planned, threads);
            assert_eq!(choice, expected, "{cell}: {layout:?}");
            assert_eq!((layout, units), (judged, model_units), "{cell}");
            let expected_units = match input {
                "a" => 1,
                "b" => 2,
                "Ad" => 18,
                _ => 89,
            };
            assert_eq!(units, expected_units, "{cell}");
            match model {
                Ok((k, m)) if k > DC_SLOTS => {
                    let taken = matches!((input, threads), ("Ad", 16 | 32 | 33));
                    assert_eq!(two_slots_win(threads, k, m, units), taken, "{cell}");
                    if !taken {
                        assert_eq!(choice, model, "{cell}");
                    }
                }
                _ => assert_eq!(choice, model, "{cell}"),
            }
        }
        // Megablast has no rule: Ad's batches at 16 threads keep the model's 4 x 4.
        let megablast = Profile::of_task("megablast").rule;
        assert_eq!(megablast, SlotRule::Model);
        assert_eq!(
            split3_choice_ruled("megablast", SPLIT3_AD, 16, megablast).0,
            Ok((4, 4))
        );
    }

    // The two-slot rule on the 24 cells of the map it was drawn from (line D session 14 window
    // w2, dev d92512d0, ws: `LOSAT_X_QBATCHPARSLOTS=2` over the model, the median of 4 rounds;
    // the model's layout and the batch units from its `LOSAT_X_STATS` line; module header): it
    // holds only on cells where two slots won (by 1.7 % or more) and on none where they lost; the
    // wins it forgoes are an odd unit count (15,600 cDNAs at 24 / 32 threads) and two full rounds
    // of 4- and 5-thread slots (10,400 at 12 / 16 threads, 20,800 at 32).
    #[test]
    fn the_two_slot_rule_holds_only_where_two_slots_won_on_the_w2_map() {
        // (input, threads, model k, model m, units, two slots' time over the model's, thousandths)
        let cells: [(&str, usize, usize, usize, usize, u32); 24] = [
            ("D5", 12, 3, 4, 3, 1118),
            ("D5", 16, 3, 5, 3, 1156),
            ("D5", 24, 3, 8, 3, 1106),
            ("D5", 32, 3, 10, 3, 1099),
            ("U4", 12, 4, 3, 4, 964),
            ("U4", 16, 4, 4, 4, 978),
            ("U4", 24, 4, 6, 4, 973),
            ("U4", 32, 2, 16, 4, 1005),
            ("U6", 12, 3, 4, 6, 986),
            ("U6", 16, 3, 5, 6, 957),
            ("U6", 24, 6, 4, 6, 946),
            ("U6", 32, 6, 5, 6, 943),
            ("U9", 12, 3, 4, 9, 1071),
            ("U9", 16, 3, 5, 9, 1013),
            ("U9", 24, 9, 2, 9, 892),
            ("U9", 32, 9, 3, 9, 920),
            ("U12", 12, 3, 4, 12, 1018),
            ("U12", 16, 4, 4, 12, 1008),
            ("U12", 24, 6, 4, 12, 1006),
            ("U12", 32, 6, 5, 12, 979),
            ("AD", 12, 3, 4, 18, 1016),
            ("AD", 16, 4, 4, 18, 983),
            ("AD", 24, 6, 4, 18, 1031),
            ("AD", 32, 9, 3, 18, 909),
        ];
        let taken = [
            ("U4", 12),
            ("U4", 16),
            ("U4", 24),
            ("U6", 24),
            ("U6", 32),
            ("AD", 16),
            ("AD", 32),
        ];
        let forgone = [("U6", 12), ("U6", 16), ("U9", 24), ("U9", 32), ("U12", 32)];
        for (input, threads, k, m, units, ratio) in cells {
            let cell = format!("{input} n{threads}");
            assert_eq!(threads / k, m, "{cell}");
            if k <= DC_SLOTS {
                // The model takes two slots itself (U4 at 32 threads): the rule is not consulted.
                assert_eq!((input, threads), ("U4", 32));
                continue;
            }
            let holds = two_slots_win(threads, k, m, units);
            assert_eq!(holds, taken.contains(&(input, threads)), "{cell}");
            if holds {
                assert!(ratio <= 983, "{cell}: {ratio}");
            } else if forgone.contains(&(input, threads)) {
                assert!(ratio < 1000, "{cell}: {ratio}");
            } else {
                assert!(ratio > 1000, "{cell}: {ratio}");
            }
        }
    }

    // The two-slot rule at its boundaries (`two_slots_win`): each clause on either side.
    #[test]
    fn the_two_slot_rule_at_its_boundaries() {
        // (threads, model k, model m, units, holds)
        let cases = [
            // Threads: 12 to 33 (32 and 33 give 16 per slot, 34 give 17); 11 and fewer, 64 and
            // 128 keep the model.
            (32, 9, 3, 18, true),
            (33, 9, 3, 18, true),
            (34, 9, 3, 18, false),
            (64, 9, 7, 18, false),
            (128, 9, 14, 18, false),
            (12, 4, 3, 4, true),
            (11, 3, 3, 2, false),
            (9, 3, 3, 2, false),
            // An even unit count.
            (16, 4, 4, 18, true),
            (16, 4, 4, 17, false),
            (16, 4, 4, 19, false),
            (32, 9, 3, 9, false),
            // Narrow slots that leave threads idle: 9 x 3 on 32 threads (27 used), 5 x 3 on 16
            // (15), 14 x 1 on 16; not 4 x 3 on 12, 8 x 3 on 24, 11 x 3 on 33 (every thread, or
            // one more than two slots use), nor 4-thread slots.
            (16, 5, 3, 20, true),
            (16, 14, 1, 28, true),
            (12, 4, 3, 12, false),
            (24, 8, 3, 16, false),
            (33, 11, 3, 22, false),
            (16, 4, 4, 16, false),
            // One round: units at most the slots.
            (24, 6, 4, 6, true),
            (24, 6, 4, 4, true),
            (24, 6, 4, 2, true),
            (24, 8, 3, 8, true),
            (24, 6, 4, 12, false),
            // A last round at most half full within five rounds: 18 units on 4 slots are 5
            // rounds (counted whole) with 2 of 4 in the last, 22 are 6, 20 fill 5; 14 on 4 slots
            // leave 2 of 4, 20 on 6 leave 2 of 6 and 28 on 8 leave 4 of 8; 28 on 6 leave 4 of 6,
            // 22 on 6 leave 4 of 6 and 30 on 8 leave 6 of 8, more than half.
            (16, 4, 4, 18, true),
            (16, 4, 4, 22, false),
            (16, 4, 4, 20, false),
            (16, 4, 4, 14, true),
            (32, 6, 5, 28, false),
            (32, 6, 5, 30, false),
            (32, 6, 5, 32, false),
            (32, 6, 5, 22, false),
            (32, 6, 5, 20, true),
            (32, 8, 4, 30, false),
            (32, 8, 4, 28, true),
        ];
        for (threads, k, m, units, holds) in cases {
            assert_eq!(threads / k, m, "{threads} threads, {k} slots");
            assert_eq!(
                two_slots_win(threads, k, m, units),
                holds,
                "{threads} threads, {k} x {m}, {units} units"
            );
        }
    }

    // The rule against the model without it (LOSAT_X_DCTWOSLOTS=0) on random and edge inputs:
    // the judged layout is the same; the choice is the same wherever the model keeps the loop,
    // does not run, takes at most two slots, the tuning variable sets the slots, or the rule
    // declines; elsewhere it is two slots of `threads / 2` threads. Every clause of the rule,
    // written out here as a chain, decides some cases: an odd unit count, more than 16 threads
    // per slot, fewer than 12 threads, narrow slots that leave threads idle, one round, a last
    // round at most half full within five rounds, and a declined even count.
    #[test]
    fn the_rule_changes_only_layouts_of_more_than_two_slots_where_it_holds() {
        let mut state = 0x9e37_79b9_7f4a_7c15u64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        let budgets = [0, 1, 50, 84, 196, 420, 868, 4096];
        let estimates = [1, 4 * MIB, 79 * MIB, 140 * MIB, u64::MAX];
        // kept, odd, wide, few threads, narrow and idle, one round, half-empty last round,
        // declined even
        let mut hits = [0usize; 8];
        for case in 0..20_000 {
            // Half the cases lie in the rule's range (4 to 33 threads, at most 40 units).
            let threads = match case % 4 {
                0 => [1, 2, 3, 16, 32, 33, 34, 128][next(8) as usize],
                1 => 1 + next(130) as usize,
                _ => 4 + next(30) as usize,
            };
            let budget = budgets[next(budgets.len() as u64) as usize] * MIB;
            let estimate = estimates[next(estimates.len() as u64) as usize];
            let units = match case % 4 {
                0 | 1 => next(120) as usize,
                _ => next(41) as usize,
            };
            let slots = (next(8) == 0).then(|| next(threads as u64 + 3) as usize);
            let points = [EFFICIENCY_MEGABLAST, EFFICIENCY_BLASTN][next(2) as usize];
            let what = format!("{threads} {budget} {estimate} {units} {slots:?}");
            let (old, old_layout) = choose_slots(
                threads,
                budget,
                estimate,
                slots,
                SlotRule::Model,
                units,
                points,
            );
            let (new, new_layout) = choose_slots(
                threads,
                budget,
                estimate,
                slots,
                SlotRule::TwoSlots,
                units,
                points,
            );
            assert_eq!(new_layout, old_layout, "{what}");
            match (old, old_layout) {
                (Ok((k, m)), Some(judged)) if k > DC_SLOTS && slots.is_none() => {
                    assert_eq!((k, m), (judged.k, judged.m), "{what}");
                    let two = Ok((DC_SLOTS, threads / DC_SLOTS));
                    let (clause, expected) = if units % 2 == 1 {
                        (1, old)
                    } else if threads / 2 > 16 {
                        (2, old)
                    } else if threads < 12 {
                        (3, old)
                    } else if m <= 3 && k * m < 2 * (threads / 2) {
                        (4, two)
                    } else if units <= k {
                        (5, two)
                    } else if units % k != 0 && 2 * (units % k) <= k && units.div_ceil(k) <= 5 {
                        (6, two)
                    } else {
                        (7, old)
                    };
                    hits[clause] += 1;
                    assert_eq!(new, expected, "{what}: clause {clause}");
                }
                _ => {
                    hits[0] += 1;
                    assert_eq!(new, old, "{what}");
                }
            }
        }
        assert!(hits.iter().all(|&hit| hit >= 30), "{hits:?}");
    }

    // The eight-slot rule at its boundaries (`eight_slots_apply`): more than eight slots in the
    // model's layout, on 12 to 64 threads.
    #[test]
    fn the_eight_slot_rule_at_its_boundaries() {
        // (threads, model k, applies)
        let cases = [
            (16, 16, true),
            (16, 9, true),
            (16, 8, false),
            (16, 2, false),
            (12, 12, true),
            (12, 9, true),
            (11, 11, false),
            (9, 9, false),
            (64, 32, true),
            (64, 9, true),
            (65, 33, false),
            (128, 89, false),
        ];
        for (threads, k, applies) in cases {
            assert_eq!(
                eight_slots_apply(threads, k),
                applies,
                "{threads} threads, {k} slots"
            );
        }
    }

    // LOSAT_X_BNEIGHTSLOTS on the measured cells (module header): -task blastn of 31,168 cDNAs
    // (89 units, Ab's batches) and of 5,000 (14 units of about 1.1 Mb) against 3 Mb. Where the
    // model takes more than eight slots on 12 to 64 threads the rule runs its best layout of at
    // most eight (8 x 2 / 8 x 4 for 89 units at 16 / 32 threads, 6 x 2 at 12, 8 x 3 at 24, 8 x 8
    // at 64; 7 x 2 / 7 x 4 / 7 x 9 for 14 units at 16 / 32 / 64 threads, a tie with 8 slots: the
    // fewest);
    // where it takes at most eight (8 and 12 threads for 14 units, 8 for 89) or on 128 threads,
    // the model's layout runs. The judged layout is the model's either way.
    #[test]
    fn the_eight_slot_rule_on_the_measured_cells() {
        let c5 = vec![1_100_000u64; 14];
        // (input, planned, threads, the model's layout, the layout with the rule)
        let cells: [(&str, &[u64], usize, (usize, usize), (usize, usize)); 13] = [
            ("AB", SPLIT3_AB, 8, (8, 1), (8, 1)),
            ("AB", SPLIT3_AB, 12, (12, 1), (6, 2)),
            ("AB", SPLIT3_AB, 16, (16, 1), (8, 2)),
            ("AB", SPLIT3_AB, 24, (24, 1), (8, 3)),
            ("AB", SPLIT3_AB, 32, (32, 1), (8, 4)),
            ("AB", SPLIT3_AB, 64, (32, 2), (8, 8)),
            ("AB", SPLIT3_AB, 128, (89, 1), (89, 1)),
            ("C5", &c5, 8, (7, 1), (7, 1)),
            ("C5", &c5, 12, (4, 3), (4, 3)),
            ("C5", &c5, 16, (14, 1), (7, 2)),
            ("C5", &c5, 32, (14, 2), (7, 4)),
            ("C5", &c5, 64, (14, 4), (7, 9)),
            ("C5", &c5, 128, (14, 9), (14, 9)),
        ];
        let rule = Profile::of_task("blastn").rule;
        for (input, planned, threads, model, ruled) in cells {
            let cell = format!("{input} {threads} threads");
            let (choice, layout, units) = split3_choice_ruled("blastn", planned, threads, rule);
            let (unruled, judged, model_units) = split3_choice("blastn", planned, threads);
            assert_eq!(unruled, Ok(model), "{cell}: {judged:?}");
            assert_eq!(choice, Ok(ruled), "{cell}: {layout:?}");
            assert_eq!((layout, units), (judged, model_units), "{cell}");
            assert_eq!(units, planned.len().min(89), "{cell}");
        }
        // blastn-short has no rule: Ab's batches at 16 threads keep the model's 16 x 1.
        let short = Profile::of_task("blastn-short").rule;
        assert_eq!(short, SlotRule::Model);
        assert_eq!(
            split3_choice_ruled("blastn-short", SPLIT3_AB, 16, short).0,
            Ok((16, 1))
        );
    }

    // The eight-slot rule against the model without it (LOSAT_X_BNEIGHTSLOTS=0) on random and
    // edge inputs: the judged layout is the same; the choice is the same wherever the model
    // keeps the loop, does not run, takes at most eight slots, the tuning variable sets the
    // slots, or the threads are fewer than 12 or more than 64; elsewhere it is the least
    // predicted span of two to eight slots (the fewest on a tie), which also pays. Each of the
    // four outcomes (kept, fewer than 12 threads, more than 64, the rule's layout) decides at
    // least 30 cases.
    #[test]
    fn the_eight_slot_rule_changes_only_layouts_of_more_than_eight_slots_where_it_holds() {
        let mut state = 0x2545_f491_4f6c_dd1du64;
        let mut next = |bound: u64| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state % bound
        };
        let budgets = [0, 1, 50, 84, 196, 420, 868, 4096];
        let estimates = [1, 4 * MIB, 50 * MIB, 140 * MIB, u64::MAX];
        // kept, fewer than 12 threads, more than 64, the rule's layout
        let mut hits = [0usize; 4];
        for case in 0..20_000 {
            // Half the cases lie around the rule's range (8 to 71 threads, at most 100 units).
            let threads = match case % 4 {
                0 => [1, 2, 11, 12, 16, 64, 65, 128][next(8) as usize],
                1 => 1 + next(130) as usize,
                _ => 8 + next(64) as usize,
            };
            let budget = budgets[next(budgets.len() as u64) as usize] * MIB;
            let estimate = estimates[next(estimates.len() as u64) as usize];
            let units = match case % 4 {
                0 | 1 => next(200) as usize,
                _ => next(101) as usize,
            };
            let slots = (next(8) == 0).then(|| next(threads as u64 + 3) as usize);
            let points = [EFFICIENCY_MEGABLAST, EFFICIENCY_BLASTN][next(2) as usize];
            let what = format!("{threads} {budget} {estimate} {units} {slots:?}");
            let (old, old_layout) = choose_slots(
                threads,
                budget,
                estimate,
                slots,
                SlotRule::Model,
                units,
                points,
            );
            let (new, new_layout) = choose_slots(
                threads,
                budget,
                estimate,
                slots,
                SlotRule::EightSlots,
                units,
                points,
            );
            assert_eq!(new_layout, old_layout, "{what}");
            match (old, old_layout) {
                (Ok((k, m)), Some(judged)) if k > BN_SLOTS && slots.is_none() => {
                    assert_eq!((k, m), (judged.k, judged.m), "{what}");
                    let (clause, expected) = if threads < 12 {
                        (1, old)
                    } else if threads > 64 {
                        (2, old)
                    } else {
                        let ratio = |k: usize| {
                            predicted_span(points, threads, k, units)
                                / predicted_span(points, threads, 1, units)
                        };
                        let mut fewer = 2;
                        for k in 3..=8 {
                            if ratio(k) < ratio(fewer) {
                                fewer = k;
                            }
                        }
                        assert!(slots_pay(ratio(fewer)), "{what}: {fewer}");
                        (3, Ok((fewer, threads / fewer)))
                    };
                    hits[clause] += 1;
                    assert_eq!(new, expected, "{what}: clause {clause}");
                }
                _ => {
                    hits[0] += 1;
                    assert_eq!(new, old, "{what}");
                }
            }
        }
        assert!(hits.iter().all(|&hit| hit >= 30), "{hits:?}");
    }

    // One long record among normal batches (review S2): Ab's 89 blastn batches of ~1 Mb with one
    // 3 Mb record keep slots (the largest batch's estimate, 121 MiB, caps k by the budget), and
    // ten 5 Mb dc-megablast batches with one 11 Mb contig too; against the largest batch both had
    // one unit and kept the loop.
    #[test]
    fn one_long_record_keeps_the_slots() {
        let mut outlier = SPLIT3_AB.to_vec();
        outlier.push(3_000_000);
        for (threads, expected, ratio) in [
            (4, (2, 2), 0.8771),
            (8, (4, 2), 0.6915),
            (16, (8, 2), 0.4276),
        ] {
            let (choice, layout, units) = split3_choice("blastn", &outlier, threads);
            assert_eq!((choice, units), (Ok(expected), 90), "{threads}: {layout:?}");
            let layout = layout.expect("a judged layout");
            assert!((layout.ratio - ratio).abs() < 1e-4, "{threads}: {layout:?}");
        }
        let mut contigs = vec![5_000_000u64; 10];
        contigs.push(11_000_000);
        let (choice, layout, units) = split3_choice("dc-megablast", &contigs, 16);
        assert_eq!((choice, units), (Ok((3, 5)), 11), "{layout:?}");
        assert!((layout.expect("a layout").ratio - 0.9158).abs() < 1e-4);
    }

    // The gate is strict, and a tie between layouts goes to the fewest slots (less memory in
    // flight). With an efficiency table that makes every m >= 2 slow, 18 units at 16 threads
    // predict the same span (2 rounds x 1.2) for 9 to 16 one-thread slots: 9 are taken.
    #[test]
    fn the_gate_is_strict_and_ties_take_the_fewest_slots() {
        assert!(!slots_pay(MODEL_GAIN));
        assert!(!slots_pay(0.9696));
        assert!(slots_pay(0.9694));
        assert_eq!(
            (milli(0.9694), milli(0.9696), milli(1.0088)),
            (969, 970, 1009)
        );
        let slow: &[(usize, f64)] = &[(1, 1.0), (2, 0.1), (16, 0.1)];
        for k in 10..=16 {
            assert_eq!(
                predicted_span(slow, 16, k, 18),
                predicted_span(slow, 16, 9, 18),
                "{k}"
            );
        }
        let tiny = batch_estimate_bytes(8, 100, 10, QUERY_BYTES_PER_LETTER);
        let (choice, layout) = choose_slots(16, 4096 * MIB, tiny, None, SlotRule::Model, 18, slow);
        assert_eq!(choice, Ok((9, 1)), "{layout:?}");
    }

    #[test]
    fn slots_follow_the_budget_the_units_and_the_tuning_variable() {
        let mb = EFFICIENCY_MEGABLAST;
        let dc = batch_estimate_bytes(11, 5_000_000, 1_000, MEGABLAST_QUERY_BYTES_PER_LETTER);
        // Too small a budget (or none) for a second batch.
        assert_eq!(
            choose_slots(16, 50 * MIB, dc, None, SlotRule::Model, 18, mb),
            (Err(Off::Budget), None)
        );
        assert_eq!(
            choose_slots(16, 0, dc, None, SlotRule::Model, 18, mb),
            (Err(Off::Budget), None)
        );
        // One unit of work, or none, whatever the budget.
        assert_eq!(
            choose_slots(16, 4096 * MIB, dc, None, SlotRule::Model, 1, mb),
            (Err(Off::Single), None)
        );
        assert_eq!(
            choose_slots(16, 4096 * MIB, dc, None, SlotRule::Model, 0, mb),
            (Err(Off::Single), None)
        );
        // One thread.
        assert_eq!(
            choose_slots(1, 4096 * MIB, dc, None, SlotRule::Model, 18, mb),
            (Err(Off::Single), None)
        );
        // A large budget does not give more slots than the model wants: 18 dc units at 16
        // threads stay at 4 x 4 (the earlier rule gave 9 x 1).
        assert_eq!(
            choose_slots(16, 4096 * MIB, dc, None, SlotRule::Model, 18, mb).0,
            Ok((4, 4))
        );
        // The budget caps k: with 84 MiB the best is 2 x 8, which the model judges 1 % slower
        // than the loop on 16 threads.
        let (choice, layout) = choose_slots(16, 84 * MIB, dc, None, SlotRule::Model, 18, mb);
        assert_eq!(choice, Err(Off::Model));
        let layout = layout.expect("a judged layout");
        assert_eq!((layout.k, layout.m), (2, 8));
        assert!(layout.ratio > 1.0 && layout.ratio < 1.01, "{layout:?}");
        // All tiny batches are units of work too: 40 one-query batches at 4 threads, 4 x 1.
        let tiny = batch_estimate_bytes(8, 100, 10, QUERY_BYTES_PER_LETTER);
        for points in [EFFICIENCY_BLASTN, mb] {
            assert_eq!(
                choose_slots(4, 4096 * MIB, tiny, None, SlotRule::Model, 40, points).0,
                Ok((4, 1))
            );
        }
        // The tuning variable wins, within the thread count, whatever the budget and units; the
        // layout it sets is judged too.
        let (choice, layout) = choose_slots(16, 0, dc, Some(8), SlotRule::Model, 18, mb);
        assert_eq!(choice, Ok((8, 2)));
        let layout = layout.expect("a judged layout");
        assert_eq!((layout.k, layout.m), (8, 2));
        let loop_span = predicted_span(mb, 16, 1, 18);
        assert!((layout.ratio - predicted_span(mb, 16, 8, 18) / loop_span).abs() < 1e-12);
        assert_eq!(
            choose_slots(4, 0, dc, Some(8), SlotRule::Model, 18, mb).0,
            Ok((4, 1))
        );
        assert_eq!(
            choose_slots(4, 0, dc, Some(2), SlotRule::Model, 1, mb).0,
            Ok((2, 2))
        );
        assert_eq!(
            choose_slots(1, 0, dc, Some(2), SlotRule::Model, 0, mb).0,
            Ok((1, 1))
        );
    }

    // The span model: rounds, contention and the in-batch efficiency.
    #[test]
    fn the_predicted_span_counts_rounds_contention_and_efficiency() {
        let bn = EFFICIENCY_BLASTN;
        // The loop: units / (n x E(n)).
        assert!((predicted_span(bn, 4, 1, 5) - 5.0 / (4.0 * 0.71)).abs() < 1e-12);
        // Up to 2k units, whole rounds: 5 units on 3 slots of one thread are 2 rounds.
        assert!((predicted_span(bn, 3, 3, 5) - 2.0 * 1.10).abs() < 1e-12);
        // Beyond, units / k plus half a round when k does not divide them.
        assert!((predicted_span(bn, 4, 4, 89) - (89.0 / 4.0 + 0.5) * 1.15).abs() < 1e-12);
        assert!((predicted_span(bn, 4, 4, 88) - 22.0 * 1.15).abs() < 1e-12);
        // The contention stops growing at five slots.
        assert!((predicted_span(bn, 16, 16, 16) - 1.20).abs() < 1e-12);
        // Idle threads count: 16 threads in 5 slots of 3.
        let e3 = efficiency(bn, 3);
        assert!((predicted_span(bn, 16, 5, 5) - 1.20 / (3.0 * e3)).abs() < 1e-12);
        // Degenerate arguments stay finite.
        assert!(predicted_span(bn, 0, 0, 0).is_finite());
    }

    // The waste bound: a plan's window (two batches per slot) plus a quarter of the batches
    // taken.
    #[test]
    fn the_waste_bound_grows_with_the_slots_and_the_work_taken() {
        assert_eq!(waste_bound(2, 0), 4);
        assert_eq!(waste_bound(4, 7), 9);
        assert_eq!(waste_bound(15, 100), 55);
    }

    /// A result for `settle`: `Err("cancelled")` is a cancelled search, another `Err` a failed
    /// one.
    type Settled = Result<u32, &'static str>;

    impl SlotResult for Settled {
        fn failed(&self) -> bool {
            matches!(self, Err(error) if *error != "cancelled")
        }

        fn cancelled(&self) -> bool {
            matches!(self, Err(error) if *error == "cancelled")
        }

        fn held_bytes(&self) -> u64 {
            7
        }
    }

    fn key_at(input_start: usize) -> BatchKey {
        BatchKey {
            input_start,
            input_end: input_start + 1,
            batch_size: 100,
        }
    }

    fn job_at(input_start: usize) -> Job {
        Job {
            key: key_at(input_start),
            searched: input_start..input_start + 1,
            bytes: 10,
        }
    }

    /// A state with one slot running batch 1, batches 1-4 wanted and 2-4 queued.
    fn settling_state() -> State<Settled> {
        let mut state = State::new();
        state.running = vec![Some((key_at(1), 10))];
        state.wanted = (1..5).map(key_at).collect();
        state.queue = (2..5).map(job_at).collect();
        state
    }

    // Review R2: the slot decides "cancelled" from its flag after the search, not from what the
    // search returned: a search that stopped part-way and returned a value that looks complete
    // is never taken, and its key, still wanted, is searched again.
    #[test]
    fn a_flagged_search_is_never_taken_whatever_it_returned() {
        let mut state = settling_state();
        let unused = settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, true);
        assert!(matches!(unused, Some(Outcome::Value(Ok(1)))));
        assert!(state.done.is_empty());
        assert_eq!((state.cancelled, state.discarded), (1, 1));
        let queued: Vec<usize> = state.queue.iter().map(|job| job.key.input_start).collect();
        assert_eq!(queued, [1, 2, 3, 4]);
        assert_eq!(state.running, [None]);
        // Not flagged: the same value is kept for the loop.
        let mut state = settling_state();
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, false).is_none());
        assert_eq!(state.done.get(&key_at(1)).map(|(_, held)| *held), Some(7));
        assert_eq!(state.cancelled, 0);
        // A key no plan wants is dropped, not queued again.
        let mut state = settling_state();
        state.wanted.remove(&key_at(1));
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Ok(1)), 7, true).is_some());
        assert_eq!(state.queue.len(), 3);
    }

    // Review R3: a slot's panic is counted, kept for the loop when its key is wanted (it then
    // ends the batches after it, as a failure does), dropped otherwise.
    #[test]
    fn slot_panics_are_counted_kept_or_not() {
        let mut state = settling_state();
        let panic = Outcome::Panic(Box::new("boom"), vec!["message\n".to_string()]);
        assert!(settle(&mut state, 0, job_at(1), panic, 0, false).is_none());
        assert_eq!(state.panicked, 1);
        assert!(state.done.contains_key(&key_at(1)));
        assert!(
            state.queue.is_empty(),
            "nothing after a failed batch starts"
        );
        let mut state = settling_state();
        state.wanted.clear();
        let panic = Outcome::Panic(Box::new("boom"), Vec::new());
        assert!(settle(&mut state, 0, job_at(1), panic, 0, false).is_some());
        assert_eq!((state.panicked, state.discarded), (1, 1));
        // A failed (not cancelled) result ends the queue after it too.
        let mut state = settling_state();
        assert!(settle(&mut state, 0, job_at(1), Outcome::Value(Err("e")), 7, false).is_none());
        assert!(state.queue.is_empty());
    }

    #[cfg(all(feature = "parallel", not(target_arch = "wasm32")))]
    mod executor {
        use super::*;
        use std::sync::atomic::{AtomicUsize, Ordering};
        use std::time::{Duration, Instant};

        /// A fake batch result: its key and the thread count of its pool; `Err("cancelled")`
        /// for a search that saw its cancellation.
        type Fake = Result<(BatchKey, usize), String>;

        impl SlotResult for Fake {
            fn failed(&self) -> bool {
                matches!(self, Err(error) if error != "cancelled")
            }

            fn cancelled(&self) -> bool {
                matches!(self, Err(error) if error == "cancelled")
            }

            fn held_bytes(&self) -> u64 {
                1
            }
        }

        /// What the fake searches do besides returning their key.
        #[derive(Clone, Copy, Default)]
        struct Fakes {
            /// The batch of this input start returns an error (the loop stops there).
            fail: Option<usize>,
            /// The batch of this input start panics.
            panic: Option<usize>,
            /// The failing batch waits until a slot has started the panicking one (at most 10 s;
            /// then it panics, so a run without slots fails instead of spinning).
            fail_waits_for_panic: bool,
            /// A batch of another size than `BATCH_SIZE` (a corrupted prediction) waits for its
            /// cancellation (at most 10 s), or returns at once.
            stale: Stale,
            /// The loop waits this long after each batch it takes, so the slots run ahead.
            loop_delay: Option<Duration>,
            /// After the first batch is committed (the slots start then), the loop waits until a
            /// slot has started a batch of another size (a corrupted prediction; at most 10 s, then
            /// it panics), so that the next take finds that batch running and cancels it. Without
            /// the wait a slow machine could start no such batch before the next take drops it
            /// from the queue, and nothing would be cancelled.
            loop_waits_for_stale: bool,
            /// The lookup width the estimates assume (8 if 0).
            lut: usize,
            /// The task profile (blastn's with a budget of 1 GiB per thread if `None`).
            profile: Option<Profile>,
        }

        /// The thread count a search reports when it saw its cancellation and returned anyway
        /// (`Stale::PartialOnCancel`): such a result must never be taken.
        const PARTIAL: usize = 9_999;

        #[derive(Clone, Copy, Default, PartialEq, Eq)]
        enum Stale {
            #[default]
            Normal,
            WaitsForCancel,
            /// Waits for its cancellation, then returns a value as if it were complete.
            PartialOnCancel,
            Quick,
        }

        #[derive(Debug)]
        struct Run {
            summary: Summary,
            results: Vec<Fake>,
            /// The keys the searches started, in order.
            started: Vec<BatchKey>,
            /// Searches running after `finish` (0: every slot was joined).
            active_after_finish: usize,
            line: String,
        }

        /// Runs the loop over `lengths` (fixed `batch_size`, or the mixer with a gie of one per
        /// ten letters) with `threads` threads, as `run_in_pool` does.
        fn run(
            lengths: &[usize],
            batch_size: Option<u32>,
            threads: usize,
            switches: Switches,
            fakes: Fakes,
        ) -> Run {
            let query_input = input_of(lengths, &vec![false; lengths.len()]);
            let searched = searched_lengths(&query_input);
            let input = PlanInput {
                query_input: &query_input,
                lengths: &searched,
                reading_error: false,
                fixed_batch_size: batch_size.is_some(),
                split: SplitSizes::default_for(false),
            };
            let active = AtomicUsize::new(0);
            let started = Mutex::new(Vec::new());
            let panic_started = AtomicBool::new(false);
            let stale_started = AtomicBool::new(false);
            let search = |job: &Job, pool: &SearchPool<'_>| -> Fake {
                started.lock().unwrap().push(job.key);
                active.fetch_add(pool.threads(), Ordering::SeqCst);
                let stale = batch_size.is_some_and(|size| job.key.batch_size != size);
                if stale {
                    stale_started.store(true, Ordering::SeqCst);
                }
                match fakes.stale {
                    Stale::WaitsForCancel | Stale::PartialOnCancel if stale => {
                        let since = Instant::now();
                        while !cancelled(pool) && since.elapsed() < Duration::from_secs(10) {
                            std::thread::yield_now();
                        }
                    }
                    Stale::Quick if stale => {}
                    _ => std::thread::sleep(Duration::from_millis(2)),
                }
                active.fetch_sub(pool.threads(), Ordering::SeqCst);
                if cancelled(pool) {
                    if fakes.stale == Stale::PartialOnCancel {
                        return Ok((job.key, PARTIAL));
                    }
                    return Err("cancelled".to_string());
                }
                if fakes.panic == Some(job.key.input_start) {
                    panic_started.store(true, Ordering::SeqCst);
                    panic!("fake batch panic at {}", job.key.input_start);
                }
                if fakes.fail == Some(job.key.input_start) && fakes.fail_waits_for_panic {
                    let since = Instant::now();
                    while !panic_started.load(Ordering::SeqCst) {
                        assert!(
                            since.elapsed() < Duration::from_secs(10),
                            "no slot started the panicking batch within 10 s"
                        );
                        std::thread::yield_now();
                    }
                }
                if pool.threads() > 1 {
                    assert_eq!(rayon::current_num_threads(), pool.threads());
                }
                if fakes.fail == Some(job.key.input_start) {
                    return Err(format!("error at {}", job.key.input_start));
                }
                Ok((job.key, pool.threads()))
            };
            let shared: Shared<Fake> = Shared::new(switches);
            let (results, active_after_finish) = std::thread::scope(|scope| {
                let mut driver = Driver::new(
                    &shared,
                    &search,
                    scope,
                    input,
                    threads,
                    fakes.profile.unwrap_or(Profile {
                        allowance_mib: 1024,
                        ..Profile::of_task("blastn")
                    }),
                    Box::new(move |_, _| if fakes.lut == 0 { 8 } else { fakes.lut }),
                    None,
                );
                let mut mixer =
                    BatchSizeMixer::new(SplitSizes::default_for(false).mixer_max_batch_size());
                let mut size = batch_size.unwrap_or_else(|| mixer.batch_size(None) as u32);
                let mut start = 0;
                let mut results = Vec::new();
                while start < lengths.len() {
                    if fakes.loop_waits_for_stale && start > 0 && results.len() == 1 {
                        let since = Instant::now();
                        while !stale_started.load(Ordering::SeqCst) {
                            assert!(
                                since.elapsed() < Duration::from_secs(10),
                                "no slot started a batch of another size within 10 s"
                            );
                            std::thread::yield_now();
                        }
                    }
                    let (end, _) = next_ranged_batch(lengths, &query_input.skipped, start, size);
                    let key = BatchKey {
                        input_start: start,
                        input_end: end,
                        batch_size: size,
                    };
                    let result = match driver.take(key) {
                        Some(result) => result,
                        None => {
                            assert_eq!(rayon::current_thread_index(), None);
                            let job = Job {
                                key,
                                searched: start..end,
                                bytes: 0,
                            };
                            search(&job, &SearchPool::x_serial())
                        }
                    };
                    let stop = result.is_err();
                    results.push(result);
                    if stop {
                        break;
                    }
                    if let Some(delay) = fakes.loop_delay {
                        std::thread::sleep(delay);
                    }
                    let letters: u64 = lengths[start..end].iter().map(|&l| l as u64).sum();
                    let gie = letters / 10;
                    if batch_size.is_none() {
                        size = mixer.batch_size(Some(gie as u32 as i32)) as u32;
                    }
                    driver.commit(key, gie, batch_size.is_none().then_some(&mixer), size);
                    start = end;
                }
                driver.finish();
                (results, active.load(Ordering::SeqCst))
            });
            let started = started.into_inner().unwrap();
            Run {
                summary: shared.summary(),
                results,
                started,
                active_after_finish,
                line: shared.stats_line(),
            }
        }

        fn on() -> Switches {
            Switches {
                on: true,
                ..Switches::default()
            }
        }

        fn keys(results: &[Fake]) -> Vec<Result<BatchKey, String>> {
            results
                .iter()
                .map(|result| result.as_ref().map(|(key, _)| *key).map_err(Clone::clone))
                .collect()
        }

        fn off_run(lengths: &[usize], batch_size: Option<u32>) -> Run {
            run(
                lengths,
                batch_size,
                4,
                Switches::default(),
                Fakes::default(),
            )
        }

        #[test]
        fn slots_give_the_loops_batches_in_order() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(300));
            let run4 = run(&lengths, Some(300), 4, on(), Fakes::default());
            let stats = run4.summary.stats;
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(stats.batches as usize, run4.results.len());
            assert_eq!(stats.precomputed + stats.loop_searched, stats.batches);
            assert!(stats.precomputed >= stats.batches - 1, "{}", run4.line);
            assert_eq!((stats.k, stats.m, stats.idle), (4, 1, 0), "{}", run4.line);
            // 13 planned batches, the last of one query (no unit): 4 x 1 is predicted at 0.82 of
            // the loop on 4 threads.
            assert!(
                run4.line
                    .contains(" units=12 rule_k=4 rule_m=1 predicted=0.81"),
                "{}",
                run4.line
            );
            assert!(run4.summary.peak_busy <= 4);
            assert_eq!(run4.active_after_finish, 0);
            assert_eq!(run4.summary.panicked, 0);
            assert_eq!(off.summary.stats.precomputed, 0);
            assert_eq!(off.summary.stats.off, Some(Off::NotSet));
        }

        #[test]
        fn sub_pools_give_the_same_batches() {
            let lengths = vec![100; 30];
            let switches = Switches {
                slots: Some(2),
                ..on()
            };
            let run2 = run(&lengths, Some(200), 4, switches, Fakes::default());
            let stats = run2.summary.stats;
            assert_eq!((stats.k, stats.m), (2, 2), "{}", run2.line);
            assert!(stats.precomputed > 0, "{}", run2.line);
            assert!(run2
                .results
                .iter()
                .skip(1)
                .any(|result| matches!(result, Ok((_, 2)))));
        }

        #[test]
        fn the_mixer_batches_are_taken_from_the_slots() {
            let lengths: Vec<usize> = (0..400).map(|index| 3_000 + (index * 17) % 1_000).collect();
            let off = off_run(&lengths, None);
            // The mixer's second batch holds most of the input (one batch unit): two slots set.
            let switches = Switches {
                slots: Some(2),
                ..on()
            };
            let run4 = run(&lengths, None, 4, switches, Fakes::default());
            assert_eq!(keys(&run4.results), keys(&off.results));
            assert_eq!(run4.summary.stats.mispredicted, 0, "{}", run4.line);
            assert!(run4.summary.stats.precomputed > 0, "{}", run4.line);
        }

        // Corrupted predictions: the slots' batches under other keys are discarded, the loop's
        // batches do not change, and two misses within eight batches (or the waste bound) stop
        // the speculation.
        #[test]
        fn mispredicted_batches_are_discarded() {
            let lengths = vec![100; 60];
            let off = off_run(&lengths, Some(200));
            for every in [1, 2, 5] {
                let switches = Switches {
                    slots: Some(2),
                    miss_every: Some(every),
                    ..on()
                };
                let missed = run(&lengths, Some(200), 4, switches, Fakes::default());
                assert_eq!(keys(&missed.results), keys(&off.results), "every {every}");
                assert!(
                    missed.summary.stats.mispredicted > 0,
                    "every {every}: {}",
                    missed.line
                );
                // Every prediction corrupted: the second batch's and third's next keys are
                // both wrong, before the waste bound can be reached. With fewer corrupted, the
                // discards may reach the waste bound first.
                let stops: &[Off] = if every == 1 {
                    &[Off::Guard]
                } else {
                    &[Off::Guard, Off::Waste]
                };
                assert!(
                    missed
                        .summary
                        .stats
                        .off
                        .is_some_and(|off| stops.contains(&off)),
                    "every {every}: {}",
                    missed.line
                );
                assert_eq!(missed.summary.panicked, 0);
            }
        }

        // A running batch that a new plan drops is cancelled, and its result is never taken.
        #[test]
        fn a_running_batch_that_the_plan_drops_is_cancelled() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(2),
                miss_every: Some(1),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::WaitsForCancel,
                loop_waits_for_stale: true,
                ..Fakes::default()
            };
            let cancelled_run = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&cancelled_run.results), keys(&off.results));
            assert!(
                cancelled_run.summary.cancelled > 0,
                "{}",
                cancelled_run.line
            );
            assert_eq!(cancelled_run.active_after_finish, 0);
        }

        // Discarded batches beyond the waste bound stop the speculation, also when no next key
        // is mispredicted (four slots plan eight batches with seven predictions; every seventh
        // corrupted makes only the last planned batch wrong, and the loop waits after each batch
        // so that the slots search it before the next plan drops it).
        #[test]
        fn waste_beyond_the_bound_stops_the_speculation() {
            let lengths = vec![100; 180];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                miss_every: Some(7),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::Quick,
                loop_delay: Some(Duration::from_millis(5)),
                ..Fakes::default()
            };
            let wasted = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&wasted.results), keys(&off.results));
            assert_eq!(wasted.summary.stats.mispredicted, 0, "{}", wasted.line);
            assert_eq!(
                wasted.summary.stats.off,
                Some(Off::Waste),
                "{}",
                wasted.line
            );
        }

        // The memory budget: with none, a slot starts only the loop's next batch, so one batch
        // runs at a time.
        #[test]
        fn the_budget_limits_the_batches_in_flight() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                mem_mib: Some(0),
                ..on()
            };
            let tight = run(&lengths, Some(200), 4, switches, Fakes::default());
            assert_eq!(keys(&tight.results), keys(&off.results));
            assert_eq!(tight.summary.peak_busy, 1, "{}", tight.line);
            assert!(tight.summary.stats.precomputed > 0, "{}", tight.line);
        }

        // An error at batch j ends the loop there; batches the slots searched past j are
        // dropped, and so is a panic of one of them, which is never resumed (batch j's search
        // waits until a slot has started the batch that panics).
        #[test]
        fn an_error_ends_the_loop_and_drops_batches_past_it() {
            let lengths = vec![100; 40];
            // Four slots set: the batch past the error must be started by a slot.
            let switches = Switches {
                slots: Some(4),
                ..on()
            };
            let fakes = Fakes {
                fail: Some(10),
                panic: Some(12),
                fail_waits_for_panic: true,
                ..Fakes::default()
            };
            let failed = run(&lengths, Some(100), 4, switches, fakes);
            assert_eq!(failed.results.len(), 11);
            assert_eq!(failed.results[10], Err("error at 10".to_string()));
            assert!(failed.summary.discarded >= 1, "{}", failed.line);
            assert_eq!(failed.active_after_finish, 0);
            // The panic of batch 12, past the error, is counted though never resumed.
            assert_eq!(failed.summary.panicked, 1, "{}", failed.line);
        }

        // Once a slot's batch returns an error, the batches planned after it are not started:
        // one slot searches the batches in order and stops at the failing one.
        #[test]
        fn no_batch_after_a_failed_one_starts() {
            let lengths = vec![100; 40];
            let switches = Switches {
                slots: Some(1),
                ..on()
            };
            let fakes = Fakes {
                fail: Some(10),
                ..Fakes::default()
            };
            let failed = run(&lengths, Some(100), 4, switches, fakes);
            assert_eq!(failed.results.len(), 11);
            let last = failed
                .started
                .iter()
                .map(|key| key.input_start)
                .max()
                .unwrap();
            assert_eq!(last, 10, "{:?}", failed.started);
        }

        // A panic of a batch the loop takes is resumed on the loop's thread.
        #[test]
        fn a_panic_of_a_taken_batch_is_resumed() {
            let lengths = vec![100; 20];
            let caught = std::panic::catch_unwind(|| {
                let fakes = Fakes {
                    panic: Some(5),
                    ..Fakes::default()
                };
                run(&lengths, Some(100), 4, on(), fakes)
            });
            let payload = caught.expect_err("the panic of batch 5 is resumed");
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .unwrap_or_default();
            assert_eq!(message, "fake batch panic at 5");
        }

        // Review R2 on whole runs: stale searches that see their cancellation and return a value
        // as if complete are counted as cancelled and never taken.
        #[test]
        fn a_value_returned_after_a_cancel_is_never_taken() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(2),
                miss_every: Some(1),
                ..on()
            };
            let fakes = Fakes {
                stale: Stale::PartialOnCancel,
                loop_waits_for_stale: true,
                ..Fakes::default()
            };
            let partial = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&partial.results), keys(&off.results));
            assert!(partial
                .results
                .iter()
                .all(|result| !matches!(result, Ok((_, PARTIAL)))));
            assert!(partial.summary.cancelled > 0, "{}", partial.line);
            assert_eq!(partial.summary.panicked, 0);
        }

        // A budget between none and ample: the estimates (lookup width 10, 4 MiB each) and a
        // budget of 5 MiB admit the loop's next batch and one more, so two of the four slots
        // search at a time while the queued batches wait.
        #[test]
        fn a_budget_admits_some_batches_and_holds_the_rest() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(200));
            let switches = Switches {
                slots: Some(4),
                mem_mib: Some(5),
                ..on()
            };
            let fakes = Fakes {
                lut: 10,
                ..Fakes::default()
            };
            let bounded = run(&lengths, Some(200), 4, switches, fakes);
            assert_eq!(keys(&bounded.results), keys(&off.results));
            assert_eq!(bounded.summary.peak_busy, 2, "{}", bounded.line);
            let one = batch_estimate_bytes(10, 200, 20, QUERY_BYTES_PER_LETTER);
            assert!(
                bounded.summary.peak_bytes <= 5 * MIB + 2 * one,
                "{}",
                bounded.line
            );
            assert!(bounded.summary.stats.precomputed > 0, "{}", bounded.line);
        }

        // The model keeps the loop where slots would not be 3 % faster: 9 planned batches at 4
        // threads with megablast's efficiency (the loop: 9 / (4 x 0.81) = 2.78 batch times; the
        // best, 2 x 2: 5 x 1.05 / (2 x 0.89) = 2.95). The summary line states the units, the
        // layout judged and its predicted ratio.
        #[test]
        fn the_model_keeps_the_loop_where_slots_gain_little() {
            let lengths = vec![100; 40];
            let off = off_run(&lengths, Some(400));
            let fakes = Fakes {
                profile: Some(Profile {
                    allowance_mib: 1024,
                    ..Profile::of_task("megablast")
                }),
                ..Fakes::default()
            };
            let kept = run(&lengths, Some(400), 4, on(), fakes);
            let stats = kept.summary.stats;
            assert_eq!(keys(&kept.results), keys(&off.results));
            assert_eq!(stats.off, Some(Off::Model), "{}", kept.line);
            assert_eq!(
                (stats.precomputed, stats.k, stats.m),
                (0, 0, 0),
                "{}",
                kept.line
            );
            assert!(
                kept.started.len() == kept.results.len(),
                "{:?}",
                kept.started
            );
            assert!(
                kept.line
                    .contains(" units=9 rule_k=2 rule_m=2 predicted=1.062 "),
                "{}",
                kept.line
            );
            assert!(kept.line.ends_with(" off=model"), "{}", kept.line);
        }

        // LOSAT_X_DCTWOSLOTS on whole runs with dc-megablast's profile at 16 threads: for 9
        // planned batches (8 units, even) the model takes 8 x 2 (one round: the rule holds) and
        // two slots of 8 threads run; for 10 (9 units, odd) it takes 3 x 5 (predicted 0.806
        // of the loop) and the rule declines, so 3 x 5 runs. The loop's batches are the same
        // either way; the slots' pools have 8 threads where the rule acted; the summary line
        // states the model's layout and the slots that ran.
        #[test]
        fn the_dc_rule_runs_two_slots_of_half_the_threads_where_it_holds() {
            let dc = Fakes {
                profile: Some(Profile {
                    allowance_mib: 1024,
                    ..Profile::of_task("dc-megablast")
                }),
                ..Fakes::default()
            };
            let ruled = Switches {
                dc_two: true,
                ..on()
            };
            // Nine batches: even units, narrow slots.
            let nine = vec![100; 36];
            let off = off_run(&nine, Some(400));
            let model = run(&nine, Some(400), 16, on(), dc);
            let two = run(&nine, Some(400), 16, ruled, dc);
            for run in [&model, &two] {
                assert_eq!(keys(&run.results), keys(&off.results), "{}", run.line);
                assert!(run.summary.stats.precomputed > 0, "{}", run.line);
            }
            let (judged, ran) = (model.summary.stats, two.summary.stats);
            assert_eq!(
                (
                    judged.k,
                    judged.m,
                    judged.rule_k,
                    judged.rule_m,
                    judged.units
                ),
                (8, 2, 8, 2, 8),
                "{}",
                model.line
            );
            assert_eq!(
                (ran.k, ran.m, ran.idle, ran.rule_k, ran.rule_m, ran.units),
                (2, 8, 0, 8, 2, 8),
                "{}",
                two.line
            );
            assert_eq!(ran.predicted_milli, judged.predicted_milli, "{}", two.line);
            let pools: HashSet<usize> = two
                .results
                .iter()
                .filter_map(|result| result.as_ref().ok().map(|&(_, threads)| threads))
                .collect();
            assert_eq!(pools, HashSet::from([1, 8]), "{}", two.line);
            // Ten batches: odd units, the rule declines; the model's 3 x 5 runs with the switch.
            let ten = vec![100; 40];
            let off = off_run(&ten, Some(400));
            let declined = run(&ten, Some(400), 16, ruled, dc);
            assert_eq!(
                keys(&declined.results),
                keys(&off.results),
                "{}",
                declined.line
            );
            let stats = declined.summary.stats;
            assert_eq!(
                (stats.k, stats.m, stats.rule_k, stats.rule_m, stats.units),
                (3, 5, 3, 5, 9),
                "{}",
                declined.line
            );
            assert_eq!(stats.predicted_milli, Some(806), "{}", declined.line);
            // Megablast's profile has no rule: nine batches keep the model's 8 x 2.
            let megablast = Fakes {
                profile: Some(Profile {
                    allowance_mib: 1024,
                    ..Profile::of_task("megablast")
                }),
                ..Fakes::default()
            };
            let kept = run(&nine, Some(400), 16, ruled, megablast);
            let stats = kept.summary.stats;
            assert_eq!((stats.k, stats.m), (8, 2), "{}", kept.line);
        }

        // LOSAT_X_BNEIGHTSLOTS on whole runs with blastn's profile at 16 threads: for 31 planned
        // batches (30 units) the model takes 15 x 1 (predicted 0.371 of the loop) and the rule
        // runs its best of at most eight slots, 8 x 2; for 21 (20 units) it takes 4 x 4 itself,
        // which runs. The loop's batches are the same either way; the summary line states the
        // model's layout and the slots that ran.
        #[test]
        fn the_eight_slot_rule_runs_at_most_eight_slots_where_the_model_takes_more() {
            let blastn = Fakes {
                profile: Some(Profile {
                    allowance_mib: 1024,
                    ..Profile::of_task("blastn")
                }),
                ..Fakes::default()
            };
            let ruled = Switches {
                bn_eight: true,
                ..on()
            };
            let many = vec![100; 124];
            let off = off_run(&many, Some(400));
            let model = run(&many, Some(400), 16, on(), blastn);
            let eight = run(&many, Some(400), 16, ruled, blastn);
            for run in [&model, &eight] {
                assert_eq!(keys(&run.results), keys(&off.results), "{}", run.line);
                assert!(run.summary.stats.precomputed > 0, "{}", run.line);
            }
            let (judged, ran) = (model.summary.stats, eight.summary.stats);
            assert_eq!(
                (
                    judged.k,
                    judged.m,
                    judged.rule_k,
                    judged.rule_m,
                    judged.units
                ),
                (15, 1, 15, 1, 30),
                "{}",
                model.line
            );
            assert_eq!(judged.predicted_milli, Some(371), "{}", model.line);
            assert_eq!(
                (ran.k, ran.m, ran.idle, ran.rule_k, ran.rule_m, ran.units),
                (8, 2, 0, 15, 1, 30),
                "{}",
                eight.line
            );
            assert_eq!(
                ran.predicted_milli, judged.predicted_milli,
                "{}",
                eight.line
            );
            // Twenty-one batches: the model's 4 x 4 runs with the switch.
            let fewer = vec![100; 84];
            let off = off_run(&fewer, Some(400));
            let kept = run(&fewer, Some(400), 16, ruled, blastn);
            assert_eq!(keys(&kept.results), keys(&off.results), "{}", kept.line);
            let stats = kept.summary.stats;
            assert_eq!(
                (stats.k, stats.m, stats.rule_k, stats.rule_m, stats.units),
                (4, 4, 4, 4, 20),
                "{}",
                kept.line
            );
            // The dc-megablast profile does not take this rule: 31 batches keep the same layout
            // with the switch on and off, LOSAT_X_DCTWOSLOTS on in both arms.
            let dc = Fakes {
                profile: Some(Profile {
                    allowance_mib: 1024,
                    ..Profile::of_task("dc-megablast")
                }),
                ..Fakes::default()
            };
            let dc_two = Switches {
                dc_two: true,
                ..on()
            };
            let dc_model = run(&many, Some(400), 16, dc_two, dc);
            let dc_ruled = run(
                &many,
                Some(400),
                16,
                Switches {
                    bn_eight: true,
                    ..dc_two
                },
                dc,
            );
            let (a, b) = (dc_model.summary.stats, dc_ruled.summary.stats);
            assert_eq!((a.k, a.m), (b.k, b.m), "{}", dc_ruled.line);
        }

        // Without a budget for a second batch the model does not run: no layout is judged.
        #[test]
        fn no_budget_judges_no_layout() {
            let lengths = vec![100; 40];
            let switches = Switches {
                mem_mib: Some(0),
                ..on()
            };
            let none = run(&lengths, Some(400), 4, switches, Fakes::default());
            let stats = none.summary.stats;
            assert_eq!(stats.off, Some(Off::Budget), "{}", none.line);
            assert_eq!(
                (stats.units, stats.rule_k, stats.predicted_milli),
                (9, 0, None)
            );
            assert!(
                none.line.contains(" rule_k=0 rule_m=0 predicted=- "),
                "{}",
                none.line
            );
        }

        #[test]
        fn one_batch_after_the_first_does_not_speculate() {
            let lengths = vec![100; 2];
            let single = run(&lengths, Some(100), 4, on(), Fakes::default());
            assert_eq!(single.summary.stats.precomputed, 0);
            assert_eq!(
                single.summary.stats.off,
                Some(Off::Single),
                "{}",
                single.line
            );
        }
    }
}
