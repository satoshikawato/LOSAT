# Threaded Wasm: atomic instructions against a stale memory bound (2026-11-01)

Release v0.2.0 GO item (`release-v020` DECISIONS 54). A dc-megablast run at 16 threads on Node
26.8.2 ended with `RuntimeError: memory access out of bounds`, exit 143 and empty output
(`results/first_trap.dc.w39.r1.err`). It was 1 of the 483 threaded Wasm runs on record, and 1 of
the 42 dc-megablast runs at 16 threads. The host already had the `memory.grow(0)` before
`wasi_thread_start` (PR #513) and the bulk-memory guard.

## Cause

The trapping instruction is `i32.atomic.rmw.add` in `get_or_init_slow`, called from
`Stealer::steal`: the `Arc::clone` of crossbeam's `Collector`, right after a plain load of the
`Arc`'s address (`scripts/oob_locate.js`).

V8 14.6 (Node v26.8.2, `deps/v8/src/wasm/`):

- `baseline/liftoff-compiler.cc`: every atomic access calls `BoundsCheckMem(..., kDoForceCheck,
  kCheckAlignment)`, a compare with the instance's memory size. A plain load or store relies on
  the trap handler and the hardware page protection, which another thread's growth has already
  changed.
- `turboshaft-graph-interface.cc`: atomic loads, stores and read-modify-writes rely on the trap
  handler; `memory.atomic.wait` and `notify` compare with the size.
- The instance's size changes when its isolate handles the growth interrupt
  (`GROW_SHARED_MEMORY`), at a stack check or a tier-up check. `Liftoff::Loop` has no stack check
  under dynamic tiering. `FutexEmulation` returns "not equal" without handling interrupts.

So a thread in Liftoff code that reads an address in memory that another thread grew a moment
ago can trap on it, although the memory is there. Which of those windows the first trap took was
not determined; the fix does not depend on it.

## Experiments

| Experiment | Arm | Result |
|---|---|---|
| V8 only (`scripts/stale_atomic.js`), Node 26.8.2, ws and this machine | plain store | ok 20/20, 10/10 |
| | atomic `rmw.add` | trap 20/20, 10/10 |
| | atomic after `memory.grow(0)` | ok 20/20, 10/10 |
| The real module with the window widened (`scripts/amp_module.js`: 100,000 `atomic.fence` before the `Once` load), 9e519214, 16 threads, `results/amp1_runs.tsv` | as it is | trap 29/40 |
| | `memory.grow(0)` after the fences | trap 0/40 |
| The same, dev 74bbd3c7e, one run with the arms in turn, `results/amp2_runs.tsv` | dev's host | trap 7/40 |
| | the fixed guard | trap 0/40 |
| Chromium 156, V8 only (`scripts/browser/`), `results/browser_stale.log` | plain store | ok 5/5 |
| | atomic | trap 5/5 |
| | atomic after `memory.grow(0)` | ok 5/5 |
| | atomic through the fixed guard | ok 5/5 |
| Firefox 151, the same page | all four arms | ok 5/5 |
| Chromium 156, the real module widened, gbdraw 900992caf's host, 8 jobs x 6 passes, `results/chromium_real_module.jsonl` | fill and copy guarded only | 30 re-runs of a pair after a trap, 2 passes failed; all 48 recorded thread stacks are `get_or_init_slow` |
| | the fixed guard | 0 re-runs, 0 failed |

Widening the window makes the race frequent and changes nothing else: the trap has the first
trap's stack (`results/amp2_dev_trap.err`). Every run that did not trap wrote the same output
(SHA-256 prefix 5c357f81bdcadee4). The trap reports printed by `scripts/mk_host.py` give the
address and the true size: 44,337,472 < 49,610,752 and 60,756,160 < 61,276,160. The address is
inside the memory.

WebKit did not start on this machine, so it has no result.

## ARM64

`results/arm64_run_41308870798.txt`, `scripts/stale_atomic.js` on GitHub runners (20 rounds):

| Runner | Node | `i32.store` | `i32.store8` | `i32.load` | atomic | atomic after `grow(0)` |
|---|---|---|---|---|---|---|
| ubuntu-24.04-arm | 24.21.0 | trap 20/20 | trap 20/20 | ok | trap 6/20 | ok |
| ubuntu-24.04-arm | 26.8.2 | trap 20/20 | ok | ok | trap 5/20 | ok |
| macos-14 | 24.21.0 | ok | ok | ok | trap 3/20 | ok |
| macos-14 | 26.8.2 | ok | ok | ok | trap 3/20 | ok |
| ubuntu-24.04 | 24.21.0 | ok | ok | ok | trap 4/20 | ok |
| ubuntu-24.04 | 26.8.2 | ok | ok | ok | trap 5/20 | ok |

**The guard does not cover plain stores on ARM64 outside macOS.** V8 bounds checks them there
(`wasm/compilation-environment.h`, `kPartialOOBWritesAreNoops = false`; both compiler tiers).
The real module's rate there has not been measured.

The atomic counts below 20/20 on the runners are read as the function tiering up to Turboshaft
during the rounds, where atomics have no compare; this was not checked.

## The fix

`guardSharedMemory` rewrites every atomic memory instruction into a call of a helper, one per
opcode and memarg: `if address + offset + size > memory.size << 16: memory.grow(0)`, then the
original instruction. dev 74bbd3c7e: 28,725 sites, 1,412 helpers, 8,298 bytes more, 0.8 s to
guard (cached by the runner).

On V8 a helper's own entry stack check usually handles the pending interrupt first, so a test
of the guard must run with `--no-wasm-stack-checks` to reach the compare. With the flag, a helper
without its compare traps in all 4 tests of `test_wasi_shared_memory.js` on Node 24.21.0 and
26.8.2, and the real helper passes.

A differential test covers all 66 atomic opcodes against original Wasm, with results, memory and trap
messages compared at addresses in bounds, at the edge, past it, misaligned and wrapped.

## Natural rate

RATE_PENDING (`results/rate1_runs.tsv`).

## Time and memory

One thread, instructions of this branch's guard against dev's host, 2 pairs per case
(`results/cost1_runs.tsv`; output hashes equal in all):

| Case | blastp | blastp-fast | megablast | blastn | dc-megablast | blastn-short | tblastn | blastx | tblastx |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| instructions | x1.0015 | x1.0005 | x1.0008 | x1.0010 | x1.0003 | x0.9978 | x0.9934 | x1.0414 | x1.0034 |

TIMING_PENDING (`results/ab2_runs.tsv`, 9 cases at 8 threads, a null arm, 8 pairs; base = #693's host).
