# Threaded Wasm: `memory access out of bounds` at thread start (2026-10-25)

The v0.2.0 GO blocker of `release-v020` (Owner 2026-10-25): nightly run 40157413236 (attempt 1,
LOSAT source = `dev` 97cf1fde5) failed `check_wasm_threading.py` at
`format/blastx-blastx-query-2/threaded-n16`:

```
RuntimeError: memory access out of bounds
    at losat-62fea96b929aa2a3.wasm.__wasi_thread_start_C (...wasm-function[7234]:0x70b775)
    at losat-62fea96b929aa2a3.wasm.wasi_thread_start (...wasm-function[7257]:0x70c9fc)
    at .../LOSAT/tests/wasi_thread_host.js:712:26
```

The same source passed in run 40157422136 and on the rerun of the failed job.

## The trapping instruction

V8 reports offsets in the module the host compiles, the output of `guardSharedMemory`
(`LOSAT/tests/wasi_shared_memory.js`), not in the artifact. In the guarded module of that run's
artifact (`losat-threaded-command.wasm`, sha256 `6e299aea...`), 0x70c9fc is the call of
`__wasi_thread_start_C` in `wasi_thread_start`, and 0x70b775 is the first instruction of
`__wasi_thread_start_C` that uses the start arguments:

```
global.get __tls_base; i32.const 1364; i32.add; local.tee 2; local.get 0
i32.atomic.store offset=20          <- trap: self->tid = tid (wasi-libc pthread_create.c)
```

`wasi_thread_start` had just read `args->stack` and `args->tls_base` with plain `i32.load`s
from the same block of memory (the start arguments sit right below the TLS block, both in the
`malloc` that `pthread_create` made for the thread).

## Repro on GitHub runners (Node 24.21.0, 4 vCPUs)

Branch `audit/wasm-oob-repro` (scratch, not merged), workflow `.github/workflows/oob-repro.yml`
there: each job builds the threaded command as `wasm-threading.yml` does, patches the host on
the runner with `probe_patch.py` (the file here is the run-3 version) and runs the failing
command with `repro.py`. `results/run<k>-<job>.json` are the jobs' summaries.

| Arm | Host | Runs | Failures |
| --- | --- | ---: | ---: |
| default | as on `dev` (`LOSAT_WASI_PRESPAWN=auto`, standby workers) | 4,800 | 14 |
| default, commit a8faa7225 (before #488) | as at that commit | 800 | 2 |
| prespawn0 | `LOSAT_WASI_PRESPAWN=0` (a new worker per thread) | 800 | 0 |
| refresh | `memory.grow(0)` right before `wasi_thread_start` | 2,800 | 0 |
| fixed host | this branch's `wasi_thread_host.js` (dafbbf0d0), no probe arm | 2,000 | 0 |

- Run 1, 40163586559: arms interleaved in each job (default, prespawn0, refresh). Its probe
  read the start arguments through `DataView(memory.buffer)` before printing and threw
  `RangeError: Offset is outside the bounds of the DataView`, which replaced the trap's
  message; these 2 failures on `dev` and the 2 before #488 are counted as failures.
- Run 2, 40165233595: the default arm only; all 6 failures are this trap. The probe's output
  was lost (a worker's `process.stderr.write` is asynchronous and the host ends the process).
- Run 3, 40166384902: one arm per job (2 default, 2 refresh); the probe writes with
  `fs.writeSync`. 6 default failures, all this trap.
- Run 4, 40168085788: the host of this branch (dafbbf0d0, cherry-picked onto the scratch
  branch) in its default arm, 2 jobs of 1,000 runs: no failure.

Runs 1 and 3 ran the default and refresh arms side by side, 2,800 runs each: 8 failures
against 0. If both arms failed at the same rate, all 8 would fall in the default arm with
probability (1/2)^8 = 0.0039. prespawn0's 0 of 800 alone is weak (2.3 failures expected at the
default rate, P = 0.10). Median run time in run 1: default 0.46 s, refresh 0.445 s (means
0.443 and 0.442 s), prespawn0 0.775 s.

Every probe line of run 3 (6 of 6) came from a standby worker (with `auto`, all 15 threads of
the command start in standby workers, so this says nothing by itself). In each, the memory
length that the worker's isolate saw after the trap ended below the start arguments, and after
`memory.grow(0)` it covered them and the TLS store:

```
{"tid":11,"startArg":13363016,"standby":10,"before":{"length":12320768},
 "after":{"length":13369344,"stack":13363008,"tls":13363040,"tidStoreInBounds":true}}
```

(the five others alike: the view ended 0.91-1.00 MiB below the start arguments, 15-32 pages
short of the length after `memory.grow(0)`).

## What the data shows, and what it does not

Shown: a standby worker (`LOSAT_WASI_PRESPAWN`, `auto` since 2026-10-18) instantiates the
module on the shared memory before its thread exists; the guest's `pthread_create` then
allocates the thread's stack and TLS, growing the memory, and hands the start arguments to that
worker. In the failing runs the worker's view of the memory was still older than that growth
after the trap: plain loads of the start arguments beyond the view's end passed, and the first
atomic store, to the new TLS, trapped. This fits V8 checking atomics against the per-isolate
size (and, an inference, plain accesses relying on the guard region). `memory.grow(0)` updates
the calling isolate's view (the header of `wasi_shared_memory.js`), and with it no run failed.

Not examined: why V8 had not updated the worker's view (the start order reaches the worker
through `Atomics`, after the growth), and whether that is a V8 defect. #488 (65,535-page
maximum) did not introduce it: 2 of 800 before it.

The browser's thread host (`web/app/src/infra/engine-worker/thread-host.ts`) starts every
thread in a worker that instantiated the module beforehand; it was not measured in a browser.

## Fix

`memory.grow(0)` right before `wasi_thread_start`: `LOSAT/tests/wasi_thread_host.js` (worker
loop: new, standby and recycled workers) and `web/app/src/infra/reactor/instance.ts`
(`instantiateThread`'s `start`, the browser's only path).

Not covered: a running thread whose atomic access reaches memory that another thread grew
before this thread's isolate updates its view. Not seen in these runs. No test checks the new
call in either host.
