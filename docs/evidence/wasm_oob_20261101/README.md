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

V8 bounds checks plain stores on ARM64 outside macOS (`wasm/compilation-environment.h`,
`kPartialOOBWritesAreNoops = false`; both compiler tiers). By that source this covers Android and
Windows on ARM64 too; only Linux was run.

The atomic counts below 20/20 on the runners are read as the function tiering up to Turboshaft
during the rounds, where atomics have no compare; this was not checked.

The real module, unwidened, on GitHub runners, blastp and blastn, AP027131 pairs, 16 threads, 100
runs per cell with the two arms in turn (`results/arm_v2_v3/`, `scripts/arm_rate.py`). v2 = the
guard with the compare before each atomic and nothing else; v3 = this branch. P(0/100 | 0.11) is
about 9e-6; 0/100 alone bounds the rate at 3 % (95 %).

| Runner, Node | blastp v2 | blastp v3 | blastn v2 | blastn v3 |
|---|---:|---:|---:|---:|
| ubuntu-24.04-arm, 24.21.0 | 29/100 | **0/100** | 0/100 | 0/100 |
| ubuntu-24.04-arm, 26.8.2 | 11/100 | **0/100** | 0/100 | 0/100 |
| ubuntu-24.04, 26.8.2 | 0/100 | 0/100 | 0/100 | 0/100 |

The first frame of the v2 traps is `dlmalloc` (one `dlfree`). Medians of wall time and max RSS of the
runs that finished are the same to 0.01 s and 1 MiB in each cell; every run that did not trap wrote the same output.
An earlier run with 40 runs per cell (run 41314393510) trapped 10/40 and 3/40 at v2 with the same
frames; megablast and dc-megablast 0/40; blastn traps nowhere. Declaring a 1 GiB initial memory
did not help (7/40, 6/40): wasi-libc's `sbrk` calls `memory.grow` whatever the host's initial size.

Chromium 156 on `ubuntu-24.04-arm`, V8 only (`scripts/browser/refresh.html`, a plain store to an
address in grown memory after an atomic read of another address, 2 jobs x 20 rounds per arm,
`results/arm_check_run_41324130443.txt`):

| Flags | Original | Through the guard |
|---|---:|---:|
| stack checks off | trap 40/40 | trap 0/40 |
| default | trap 36/40 | trap 0/40 |

Not run: a real search in a browser on ARM64, and any Android device.

## The fix

`guardSharedMemory` rewrites, per module:

- **Every atomic memory instruction** into a call of a helper, one per opcode and memarg: `if
  address + offset + size > memory.size << 16: memory.grow(0)`, then the original instruction.
- **Every `memory.grow`** into a call of a helper: the original `memory.grow`, then, when it
  succeeded, `GROWN_PAGES = max(GROWN_PAGES, memory.size)`, with a compare-exchange loop.
- **The helpers of atomics that read memory** (all but stores and notify), after the instruction:
  `if memory.size < GROWN_PAGES: memory.grow(0)`.

Plain stores stay as they are: a module has about 200,000 multi-byte ones. They do not need
wrapping, because a thread can only have an address in pages another thread grew through
something that thread wrote, and the atomic read that orders the two (a lock, or the published
address) now brings the bound up to date. `memory.size` reads the same per-instance bound as the
compare (Turboshaft does not cache the size of a shared growable memory, `InstanceCache` in
`turboshaft-graph-interface.cc`), so with all growth going through the helper, `memory.size <
GROWN_PAGES` means that this instance is behind. The hosts grow only with `memory.grow(0)`, and by
one page when they create the memory, before any instance exists. Should a compiler reuse an
earlier read of the size, that read is smaller or equal, which costs one more `memory.grow(0)`.

dev 74bbd3c7e: 28,725 atomic sites, 1,412 helpers, 1 `memory.grow` site; 39,016 bytes more than
the unguarded module with the bulk sites, 0.8 s to guard (cached by the runner).

`GROWN_PAGES` is address 16. rustc links this target with `--stack-first`, so the low addresses
are the deep end of the main thread's stack: the stack pointer starts at 1,048,576. 7 real runs (blastp twice, blastn, megablast, dc-megablast,
tblastn, blastx; `results/low_memory_after_runs.txt`): the deepest stack was 207,504 bytes, bytes
0-15 stayed zero, and the cell held the final page count in each. A main stack that came within
20 bytes of overflowing would overlap the cell.

On V8 a helper's own entry stack check usually handles the pending interrupt first, so a test
of the guard must run with `--no-wasm-stack-checks` and `--liftoff-only` to reach the compare.
With the flags, these mutants of the guard fail the tests of `test_wasi_shared_memory.js` on Node 24.21.0 and 26.8.2
(`scripts/mutants_v3.sh`):

| Mutant | Fails |
|---|---|
| no compare before the atomic; before a store of one width; before the instructions with two operands; offset left out of the compare | the every-opcode test |
| no refresh after the atomic; after wait64; after i64 results; after one cmpxchg | the every-opcode test |
| `memory.grow` does not raise the cell | the every-opcode test, the `memory.grow` test |
| the cell is overwritten whatever its value | the `memory.grow` test |
| no retry when the compare-exchange of the raise lost | **no test** |

The last mutant survives: losing the race needs another thread between two instructions of the
helper, which a test cannot force, and a stress test of 4 threads over 300 rounds did not hit it.
A lost raise leaves the cell too small until the next `memory.grow`, so threads go without the
refresh for that long.

A differential test covers all 66 atomic opcodes against original Wasm, with results, memory and trap
messages compared at addresses in bounds, at the edge, past it, misaligned and wrapped.

## Natural rate

RATE_PENDING (`results/rate1_runs.tsv`).

## Time and memory

One thread, instructions, output hashes equal in all, round 0 not used: 2 pairs per case
(`results/cost1_runs.tsv`: v2 against dev's host; `results/cost2_runs.tsv`: v3 against v2):

| Case | blastp | blastp-fast | megablast | blastn | dc-megablast | blastn-short | tblastn | blastx | tblastx |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| v2 / dev, instructions | x1.0015 | x1.0005 | x1.0008 | x1.0010 | x1.0003 | x0.9978 | x0.9934 | x1.0414 | x1.0034 |
| v3 / v2, instructions | x0.9990 | x1.0002 | x1.0011 | x1.0003 | x1.0002 | x1.0001 | x0.9984 | x1.0007 | x1.0031 |

Blastx: cycles x1.009-1.017 for its x1.041 instructions; 72 % of its cycles are in
`x_word_finder_direct_sent` (`LOSAT/src/algorithm/blastx/x_seed_scan.rs`), where an atomic load
sits in the hot path (read as a `Once` check, not confirmed).

TIMING_PENDING (`results/ab2_runs.tsv`, 9 cases at 8 threads, a null arm, 8 pairs; base = #693's host).
