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
the runner with `probe_patch.py` and runs the failing command with `repro.py`, arms
interleaved. `results/run<k>-<job>.json` are the jobs' summaries.

| Arm | Host | Runs | Failures |
| --- | --- | ---: | ---: |
| default | as on `dev` (`LOSAT_WASI_PRESPAWN=auto`, standby workers) | 4,800 | 14 |
| default, commit a8faa7225 (before #488) | as at that commit | 800 | 2 |
| prespawn0 | `LOSAT_WASI_PRESPAWN=0` (a new worker per thread) | 800 | 0 |
| refresh | `memory.grow(0)` right before `wasi_thread_start` | 2,800 | 0 |

Runs 40163586559 (run 1: the probe itself threw on its first read, see below), 40165233595
(run 2: every failure is this trap; the probe's output was lost because a worker's
`process.stderr.write` is asynchronous and the host ends the process) and 40166384902 (run 3:
the probe writes with `fs.writeSync`). With 14 failures among 7,600 default and refresh runs,
the chance that none falls in the 2,800 refresh runs is (4,800 / 7,600)^14 = 0.0016.

Every probe line of run 3 (6 of 6) came from a standby worker, and the memory length that the
worker's isolate saw after the trap was below the start arguments; after `memory.grow(0)` it
covered them and the TLS store:

```
{"tid":11,"startArg":13363016,"standby":10,"before":{"length":12320768},
 "after":{"length":13369344,"stack":13363008,"tls":13363040,"tidStoreInBounds":true}}
```

(the five others alike, 1-2 MiB short). In run 1 the probe read the start arguments through
`DataView(memory.buffer)` before printing and threw `RangeError: Offset is outside the bounds
of the DataView` for the same reason.

## Cause

A standby worker (`LOSAT_WASI_PRESPAWN`, `auto` since 2026-10-18) instantiates the module on
the shared memory before its thread exists. The guest's `pthread_create` then allocates the
thread's stack and TLS, which grows the memory, and hands the start arguments to that worker.
V8 updates another isolate's view of a grown shared memory later, and checks an atomic access
against that view, while plain accesses rely on the guard region. So the plain loads of the
start arguments pass and the first atomic store, to the new TLS, traps. It is not an engine
bug, and #488 (65,535-page maximum) did not cause it (2 of 800 before it).

The browser's thread host (`web/app/src/infra/engine-worker/thread-host.ts`) starts every
thread in a worker that instantiated the module beforehand; it was not measured in a browser.

## Fix

`memory.grow(0)` right before `wasi_thread_start`, which updates the calling isolate's view at
once: `LOSAT/tests/wasi_thread_host.js` (worker loop) and `web/app/src/infra/reactor/instance.ts`
(`instantiateThread`'s `start`).

Not covered: a running thread whose atomic access reaches memory that another thread grew
before this thread's isolate updates its view. Not seen in these runs.
