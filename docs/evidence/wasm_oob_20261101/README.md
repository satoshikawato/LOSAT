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

The real module, unwidened, on GitHub runners, AP027131 pairs, 16 threads, 100 runs per cell with
the arms in turn (`scripts/arm_rate.py`), this branch against dev's guard, which covers `memory.fill`
and `memory.copy` only (`results/arm_dev_fix/`, run 41326262707, the fix branch at 7b4fe6132, dev at
d2ef1913f):

| Runner, Node | blastp dev | blastp this branch | blastn dev | blastn this branch |
|---|---:|---:|---:|---:|
| ubuntu-24.04-arm, 24.21.0 | **100/100** | **0/100** | 0/100 | 0/100 |
| ubuntu-24.04-arm, 26.8.2 | **64/100** | **0/100** | 0/100 | 0/100 |

P(0/100 | 0.64) is about 1e-44; 0/100 alone bounds the rate at 3 % (95 %). So on ARM64 Linux
threaded blastp at 16 threads mostly did not finish on dev: the runs that trapped wrote no output.

The same with two guards (`results/arm_v2_v3/`, run 41318924971). v2 = the compare before each
atomic and nothing else (more function entries, so more stack checks); v3 = this branch's guard.

| Runner, Node | blastp v2 | blastp v3 | blastn v2 | blastn v3 |
|---|---:|---:|---:|---:|
| ubuntu-24.04-arm, 24.21.0 | 29/100 | **0/100** | 0/100 | 0/100 |
| ubuntu-24.04-arm, 26.8.2 | 11/100 | **0/100** | 0/100 | 0/100 |
| ubuntu-24.04, 26.8.2 | 0/100 | 0/100 | 0/100 | 0/100 |

The first frame of the traps is `dlmalloc` (one `dlfree` at v2). Medians of wall time and max RSS of the
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

A real search in Chromium on ARM64 is under "After merge". No Android device has been run.

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
0-15 stayed zero, and the cell held the final page count in each. The threaded reactor, the same
way (`results/low_memory_reactor.txt`, `scripts/reactor_low_memory.js`; blastn, megablast, blastp,
blastx; its Web API has no dc-megablast or tblastn): deepest stack 166,016 bytes, the same two
readings. A main stack that came within 20 bytes of overflowing would overlap the cell.

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

x64, Node 26.8.2, the first trap's form (dc-megablast, 16 threads, AP027131 x AP027132, 16 CPUs of
an otherwise idle node), dev's host without this guard, modules 9e519214 and 74bbd3c7e in turn
(`results/rate1_runs.tsv`):

| Module | Runs | Trapped | Output hashes |
|---|---:|---:|---:|
| 9e519214 | 1,000 | 0 | 1 |
| 74bbd3c7e | 1,000 | 0 | 1 |

Another line's process filled the machine's memory around 01:00-01:04, during this run; only traps
are read here. So on this machine, in these conditions, the rate is below 0.15 % per run (95 %); the first trap
was 1 of 42 such runs on record, and P(0/2000 | 1/42) is about 1e-21, so that sighting's conditions
were not these. It is too rare for a natural run to show a fix on x64: 0/2000 after the fix could
not be told from 0/2000 before it. The tests that do separate the arms are the widened module
(7/40 against 0/40) and ARM64's real blastp (above).

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

Time and max RSS at 8 threads, in one run with the arms in turn, x64, Node 26.8.2 (`results/ab2_runs.tsv`,
`scripts/wo_ab2.sh`): this branch's guard and a null arm (the same host again) against the host of
#693, 8 CPUs of an otherwise idle node, each arm with its own guard cache, round 0 not used, output
hashes equal in all. Another line's process filled the machine's memory around 01:59-02:02, inside
round 1, so that round is left out: 7 triples per case. Keeping it changes the wall means by 0.0004
at most, and moves the null arm's RSS mean from x1.007 to x1.001.

| Case | dev s | wall, fix/dev | wall, null/dev | RSS, fix/dev | RSS, null/dev |
|---|---:|---:|---:|---:|---:|
| blastp | 11.10 | x0.998 | x1.005 | x1.001 | x1.007 |
| blastp-fast | 14.29 | x1.006 | x0.999 | x1.018 | x1.003 |
| megablast | 11.08 | x1.007 | x1.001 | x1.018 | x1.035 |
| blastn | 14.77 | x1.008 | x1.000 | x1.010 | x0.996 |
| dc-megablast | 3.99 | x1.003 | x1.018 | x1.025 | x1.008 |
| blastn-short | 4.60 | x1.000 | x1.000 | x1.009 | x1.007 |
| tblastn | 3.42 | x0.994 | x1.006 | x0.990 | x0.999 |
| blastx | 2.25 | x1.009 | x0.996 | x1.005 | x0.998 |
| tblastx | 14.41 | x1.010 | x0.999 | x0.983 | x1.013 |
| geometric mean | | x1.004 | x1.003 | x1.007 | x1.007 |

The geometric means match the null arm. Single cases differ by up to about 1 % in time; megablast's
7 fix/dev ratios are all above 1 (x1.002-1.018), so that one is likely a real small cost.

## Every mode, window widened

The widened module of 74bbd3c7e, 16 threads, 16 logical CPUs of node 0 shared with other lines' runs,
10 runs per arm per mode, arms in turn; no time read (`results/modes_runs.tsv`, `scripts/wo_modes.sh`):

| Case | blastp | blastp-fast | megablast | blastn | dc-megablast | blastn-short | tblastn | blastx | tblastx | all |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| #693's host, trapped | 8 | 0 | 0 | 1 | 0 | 4 | 0 | 0 | 0 | 13/90 |
| this branch, trapped | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0/90 |

One-sided Fisher p = 8e-5. Every run that finished wrote one output per mode. The widening sits in one
function, on the path that the first trap took; a mode at 0 does not show that mode is free of
the race without it, since dc-megablast trapped 7/40 in an earlier run on other CPUs.

## After merge

The clearing runs of the release's GO item (release DECISIONS 54, line DECISIONS 15). dev
f62da7f6b, with #703; "before" is dev just before #703 (7aa5329e7), so its guard covers `memory.fill`
and `memory.copy` only. Arms in turn, in one run.

| Case | Before | After | P(0/N at the rate before) |
|---|---:|---:|---:|
| ubuntu-24.04-arm, Node 24.21.0, real blastp, 16 threads, N=100 | 100/100 | **0/100** | 0 |
| ubuntu-24.04-arm, Node 26.8.2, real blastp, 16 threads, N=100 | 76/100 | **0/100** | 3e-62 |
| ubuntu-24.04-arm, both Nodes, real blastn, 16 threads, N=100 | 0/100, 0/100 | 0/100, 0/100 | no information |
| x64, Node 26.8.2, widened module, blastp, 16 threads, N=40 | 34/40 | **0/40** | 7e-34 |
| x64, Node 26.8.2, widened module, dc-megablast, 16 threads, N=40 | 1/40 | 0/40 | 0.36 (at the earlier 7/40: 5e-4) |

(`results/clear_arm/`, run 41358045153; `results/clear_x64_runs.tsv`, `scripts/wo_clear.sh`, 8 cores of
node 0.) The widened module is 74bbd3c7e's with merged dev's host, since the guard is what changes.
The widened dc-megablast rate depends on the CPUs: 7/40 and 0/10 earlier, 1/40 here, so that row
says little alone; blastp is the row that separates the arms on x64.

The artifacts of dev f62da7f6b (wasm-threading run 41358038188; the reverse build gives the same
bytes): `losat-threaded-command.wasm` sha256 2f4cab62...304263, `losat-threaded-reactor.wasm`
a8b77914...630a2d, record version 3.

**Real searches in Chromium on ARM64**, run by gbdraw (gbdraw Actions 41348561507, its host, GitHub
`ubuntu-24.04-arm`, 4 vCPUs, Chromium 156.0.8078.4, blastp AP027131 x AP027132 CDS translations,
5 searches per pass, 56 passes per cell, arms in turn; 16 threads means `hardwareConcurrency` spoofed
on the 4 vCPUs). Control: LOSAT dde737f18 with gbdraw's own guard. Trial: this branch at 7b4fe6132,
which has the same guard code as dev f62da7f6b but older Rust.

| Threads | Arm | Passes with a trap | Failed passes | Searches completed |
|---:|---|---:|---:|---:|
| 4 | control | 11/56 | 3 | 265/280 |
| 4 | trial | **0/56** | 0 | 280/280 |
| 16 | control | 56/56 | 56 | 0/280 |
| 16 | trial | **0/56** | 0 | 280/280 |

Every completed search in both arms gave the same 7,532 rows. One pair and one program, on few CPUs.

## What users see

- **Command on Node:** `RuntimeError: memory access out of bounds` with a Wasm stack on stderr, exit
  143, no output. On dev: on ARM64 Linux, most blastp runs at 16 threads; on x64, rarely (above).
- **Chromium:** the same trap. In gbdraw's host on x64 the widened real module printed "running the
  pair again" 30 times in 48 jobs, and 2 of 6 passes failed. On ARM64 the real search at 16 threads never
  completed, and at 4 threads 3 of 56 passes failed ("After merge").
- **LOSAT Web:** read from the code, not run on a trap: a trap in a thread ends the engine and the
  search fails with "The engine stopped: memory access out of bounds"
  (`web/app/src/infra/engine-worker/`, `reactor/abi.ts`).
- **Firefox 151:** no trap on the V8-level page. **WebKit:** did not start here, so unknown.
- **With this fix:** none of these in the runs above. No Android device has been run.

