//! PROBE (line ST s2, branch probe/st-bnahead, never a PR): `LOSAT_PROBE_BNAHEAD=1`, with
//! `LOSAT_X_QBATCHPARSLOTS=2`: the slots of LOSAT_X_QBATCHPAR become one-thread set-up slots. A
//! slot sets its batch up (records, DUST, encoding, Karlin blocks, lookup table) on its own thread
//! with its one-thread pool, then waits at a gate and searches the batch on the loop's pool, one
//! batch at a time and in batch order, while the other slot sets the next batch up. The loop's
//! thread (the pool's thread 0) runs the pool's jobs while it waits for a slot. With
//! `LOSAT_PROBE_BNAHEAD=2` the gate opens at the end of a batch's preliminary search, so the next
//! batch's preliminary search runs beside its traceback. No NCBI counterpart: scheduling only.

use std::cell::Cell;
use std::collections::BTreeSet;
use std::sync::{Condvar, Mutex, MutexGuard, OnceLock};
use std::time::Duration;

use crate::utils::threading::SearchPool;

/// 0: off; 1: the gate is held for the whole search; 2: released at the end of the prelim.
pub(crate) fn mode() -> u8 {
    static MODE: OnceLock<u8> = OnceLock::new();
    *MODE.get_or_init(|| match std::env::var("LOSAT_PROBE_BNAHEAD").as_deref() {
        Ok("1") => 1,
        Ok("2") => 2,
        _ => 0,
    })
}

pub(crate) fn on() -> bool {
    mode() != 0
}

#[derive(Default)]
struct GateState {
    busy: bool,
    registered: BTreeSet<usize>,
}

/// One search on the loop's pool at a time, the registered batch of the lowest first record first.
#[derive(Default)]
pub(crate) struct Gate {
    state: Mutex<GateState>,
    wake: Condvar,
}

impl Gate {
    fn lock(&self) -> MutexGuard<'_, GateState> {
        self.state.lock().unwrap_or_else(|e| e.into_inner())
    }

    fn register(&self, start: usize) {
        self.lock().registered.insert(start);
        self.wake.notify_all();
    }

    fn deregister(&self, start: usize) {
        self.lock().registered.remove(&start);
        self.wake.notify_all();
    }

    /// Waits for the gate; `None` once `cancelled` says so.
    pub(crate) fn acquire(&self, start: usize, cancelled: &dyn Fn() -> bool) -> Option<Held<'_>> {
        let mut state = self.lock();
        loop {
            if cancelled() {
                return None;
            }
            if !state.busy && state.registered.first() == Some(&start) {
                state.busy = true;
                return Some(Held {
                    gate: self,
                    start,
                    open: true,
                });
            }
            state = self
                .wake
                .wait_timeout(state, Duration::from_millis(1))
                .unwrap_or_else(|e| e.into_inner())
                .0;
        }
    }
}

/// The gate held by one batch; released (and the batch deregistered) when dropped or `release`d.
pub(crate) struct Held<'a> {
    gate: &'a Gate,
    start: usize,
    open: bool,
}

impl Held<'_> {
    pub(crate) fn release(&mut self) {
        if self.open {
            self.open = false;
            let mut state = self.gate.lock();
            state.busy = false;
            state.registered.remove(&self.start);
            drop(state);
            self.gate.wake.notify_all();
        }
    }
}

impl Drop for Held<'_> {
    fn drop(&mut self) {
        self.release();
    }
}

/// What a slot's search needs: the loop's pool, the gate and the batch's first record.
pub(crate) struct Ctx<'a> {
    pub(crate) main: &'a SearchPool<'a>,
    pub(crate) gate: &'a Gate,
    pub(crate) start: usize,
}

thread_local! {
    static CTX: Cell<*const ()> = const { Cell::new(std::ptr::null()) };
}

/// Runs a slot's search `f` with `ctx` visible to `current` on this thread.
pub(crate) fn with<R>(ctx: &Ctx<'_>, f: impl FnOnce() -> R) -> R {
    struct Reset<'a>(&'a Gate, usize);
    impl Drop for Reset<'_> {
        fn drop(&mut self) {
            CTX.with(|c| c.set(std::ptr::null()));
            self.0.deregister(self.1);
        }
    }
    ctx.gate.register(ctx.start);
    let _reset = Reset(ctx.gate, ctx.start);
    CTX.with(|c| c.set(ctx as *const Ctx<'_> as *const ()));
    f()
}

/// The context of the slot's search running on this thread, if any.
pub(crate) fn current() -> Option<&'static Ctx<'static>> {
    let p = CTX.with(|c| c.get());
    // SAFETY (probe): set only inside `with`, which outlives the search that reads it.
    (!p.is_null()).then(|| unsafe { &*(p as *const Ctx<'static>) })
}

/// A closure sent to the loop's pool from the slot's thread, which waits for it (`install`).
pub(crate) struct AssertSend<F>(pub(crate) F);
// SAFETY (probe): the slot's thread blocks in `install` until the closure has run.
unsafe impl<F> Send for AssertSend<F> {}
