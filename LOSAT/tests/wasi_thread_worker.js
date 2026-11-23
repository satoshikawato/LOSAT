"use strict";
// The thread worker of the threaded WASI host (wasi_thread_host.js), and what a guest thread
// needs from the host: the WASI context, the imports and the instance on the shared memory. A
// Node worker loads only this file (and wasi_artifact.js and the file table, wasi_files.js), not
// the host: every worker of a 16-thread command boots at the same time, and what each one loads
// and compiles adds to the time until the command's threads can start.
const { isMainThread, parentPort, workerData } = require("worker_threads");
const { WASI } = require("wasi");
const fs = require("fs");
const { unsignedWasiImports, blockingStdioWrite } = require("./wasi_artifact");
// The table of the files the command opens (wasi_files.js).
const { sharedFileImports, directoryReads } = require("./wasi_files");

const DEBUG = process.env.LOSAT_WASI_THREADS_DEBUG === "1";

function debug(message) {
  if (DEBUG) {
    // A failing diagnostic write must never end a worker or stop the SIGTERM of failWorker.
    try { fs.writeSync(2, `[losat-wasi-threads] ${message}\n`); } catch { /* diagnostics only */ }
  }
}

function writeError(err) {
  const message = err && err.stack ? err.stack : String(err);
  fs.writeSync(2, `${message}\n`);
}

// NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-266
// A Wasm trap cannot unwind Rust/Rayon's shared state. The main thread may be
// blocked inside Wasm, unable to receive a Worker error/exit event. Terminate
// this command process from the failing worker after flushing its diagnostic;
// process.exit() in a Node Worker would terminate only that worker.
function failWorker(err) {
  try { writeError(err); } finally { process.kill(process.pid, "SIGTERM"); }
}

// Control words shared with a thread worker: the start confirmation (every worker); for a
// worker that waits for its thread (a LOSAT_WASI_PRESPAWN standby worker, or a reactor's
// worker whose earlier thread returned, LOSAT_WASI_RECYCLE), the start order with the tid and
// startArg; and for a standby or recycled worker, whether it waits for a thread (SLOT_STATE).
const START_STATE = 0, START_ORDER = 1, START_TID = 2, START_ARG = 3, SLOT_STATE = 4, START_WORDS = 5;
// SLOT_STATE: a thread runs (or the worker boots for its first thread), the worker waits for a
// thread, the host gave it one, or its thread returned and it instantiates the module again.
const SLOT_RUNNING = 0, SLOT_READY = 1, SLOT_TAKEN = 2, SLOT_PREPARING = 3;

function markWorkerStart(control, state) {
  const view = new Int32Array(control);
  Atomics.store(view, START_STATE, state);
  Atomics.notify(view, START_STATE, 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
// A standby or recycled worker blocks here, after instantiation, until spawnThread assigns
// it a thread. The order is consumed, so a recycled worker waits again after this thread.
function waitForStartOrder(control) {
  const view = new Int32Array(control);
  while (Atomics.load(view, START_ORDER) === 0) Atomics.wait(view, START_ORDER, 0);
  const order = { tid: Atomics.load(view, START_TID), startArg: Atomics.load(view, START_ARG) };
  Atomics.store(view, START_ORDER, 0);
  return order;
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:177-188
// Node 18/24 lib/wasi.js initialize binds memory and optionally calls _initialize.
// Child entry has its own thread-pointer initialization; do not rerun main CRT.
function initializeWasiThreadInstance(wasi, instance) {
  wasi.initialize({ exports: { memory: instance.exports.memory } });
}

// `directory`: a LOSAT Web reactor's run directory, preopened as "/work" instead of "/" (docs/web/abi_v2.md
// §3, web/adapter/src/run.rs RUN_DIRECTORY, which the reactor makes its current directory; a search of -db
// opens its database there by the argv's names, line WQA D12a, 2026-11-20).
function createWasi(args, directory) {
  return new WASI({
    version: "preview1",
    // NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:172-176
    // Handle the command's exit explicitly on every supported Node version.
    returnOnExit: true,
    args,
    env: { ...process.env, PWD: process.cwd() },
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729
// Preserve lexical filenames at the recorded host cwd.
// Root preopen permits absolute paths; PWD initializes the guest cwd for
// relative filenames through std::env::set_current_dir, without rewriting argv.
    preopens: directory === undefined ? { "/": "/" } : { "/work": directory },
  });
}

// NCBI reference (598d8ae6): c++/src/algo/blast/format/blast_format.cpp:68-93
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...) : ..., m_Outfile(outfile),
// A module may import host functions besides WASI, such as the output stream of the
// LOSAT Web adapter (docs/web/abi_v2.md). The main instance receives the functions the
// caller provides; a worker instance imports the same names but never calls them, so
// each one there traps if it is called.
function hostImports(module, provided = {}) {
  const imports = {};
  for (const { module: name, name: field, kind } of WebAssembly.Module.imports(module)) {
    if (kind !== "function" || ["env", "wasi", "wasi_snapshot_preview1"].includes(name)) continue;
    imports[name] = imports[name] || {};
    imports[name][field] = (provided[name] && provided[name][field]) || (() => {
      throw new Error(`${name}.${field} is not available on this thread`);
    });
  }
  return imports;
}

function createImports(wasi, memory, spawnThread, extra, files, directory) {
  const wasiImport = unsignedWasiImports(wasi.wasiImport);
  // A folder opened for reading fails its reads with EISDIR, as on Linux (directoryReads, wasi_files.js).
  const local = directoryReads({ ...wasiImport, fd_write: blockingStdioWrite(wasiImport.fd_write) }, memory, directory ?? "/");
  return {
    ...extra,
    env: { memory },
    wasi: { "thread-spawn": spawnThread },
    wasi_snapshot_preview1: files ? { ...local, ...sharedFileImports(files, memory, local, directory) } : local,
  };
}

async function instantiateWithMemory(module, memory, args, spawnThread, provided, wasi, files = null, directory = undefined) {
  wasi ??= createWasi(args, directory);
  const instance = await WebAssembly.instantiate(
    module,
    createImports(wasi, memory, spawnThread, hostImports(module, provided), files, directory),
  );
  return { instance, wasi };
}

module.exports = {
  DEBUG, debug, writeError, START_STATE, START_ORDER, START_TID, START_ARG, SLOT_STATE, START_WORDS,
  SLOT_RUNNING, SLOT_READY, SLOT_TAKEN, SLOT_PREPARING, instantiateWithMemory,
};

// Only as a worker's entry: the host loads this file for the functions above.
if (!isMainThread && require.main === module) {
  const { control, fault, standby, recycle, early } = workerData;
  // An early worker (createThreadHost) boots before the module and the memory exist and gets
  // them in one message (`init`) from prespawnWorkers.
  let { module: workerModule, wasmPath, args, memory, files, directory } = workerData;
  // A standby worker (LOSAT_WASI_PRESPAWN) learns its tid and startArg from its start order,
  // and so does a recycled worker (LOSAT_WASI_RECYCLE) for every thread after its first.
  let { startArg, tid } = workerData;
  const name = standby === undefined ? `tid=${tid}` : `standby=${standby}`;
  // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-266
  process.on("uncaughtException", failWorker);
  // NCBI reference (598d8ae6): c++/src/app/blast/blast_app_util.hpp:260-266
  // This handler runs in the worker even while the caller is blocked in Wasm.
  let completed = false;
  process.on("exit", code => {
    if (!completed) failWorker(new Error(`WASI worker ${tid ?? `standby ${standby}`} exited unexpectedly with code ${code}`));
  });

  Promise.resolve()
    .then(async () => {
      debug(`worker boot ${name}`);
      if (fault === "startup-timeout") Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0);
      if (fault === "startup-error") throw new Error("injected worker startup error");
      if (early) {
        ({ module: workerModule, wasmPath, args, memory, files, directory } =
          await new Promise(resolve => parentPort.once("message", resolve)));
      }
      const module = workerModule || await require("./wasi_thread_host").compileThreadedModule(wasmPath);
      const spawnThread = () => -1;
      const view = new Int32Array(control);
      let waits = standby !== undefined, wasi = null;
      for (;;) {
        debug(waits ? `worker instantiate ${name}` : `worker instantiate tid=${tid} start_arg=${startArg}`);
        // A new instance on the shared memory for every guest thread, recycled worker or not. A
        // recycled worker keeps its one WASI context (bound once to the same shared memory): a
        // new context per thread would hold one more open preopen descriptor for each thread.
        const created = wasi === null;
        const instantiated = await instantiateWithMemory(module, memory, args, spawnThread, undefined, wasi ?? undefined, files, directory);
        const { instance } = instantiated;
        wasi = instantiated.wasi;
        if (created) initializeWasiThreadInstance(wasi, instance);
        if (typeof instance.exports.wasi_thread_start !== "function") {
          throw new Error(`${wasmPath} does not export wasi_thread_start`);
        }
        if (waits) {
          // A recycled worker that waits can be claimed, and a standby worker that waits needs no
          // start confirmation (the host's startStandby).
          Atomics.store(view, SLOT_STATE, SLOT_READY); Atomics.notify(view, SLOT_STATE);
          debug(`worker standby ${name}`);
          ({ tid, startArg } = waitForStartOrder(control));
          debug(`worker ${name} takes tid=${tid} start_arg=${startArg}`);
        }
        debug(`worker start tid=${tid}`);
        markWorkerStart(control, 1);
        if (fault === "trap") throw new WebAssembly.RuntimeError("injected worker trap");
        if (fault === "abnormal-exit") process.exit(23);
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
        // Host memory only. This instance may be older than the guest's last growth of the
        // shared memory: a standby worker (or a recycled one, not measured) instantiated
        // before pthread_create allocated this thread's stack and TLS and grew the memory for
        // them. In 12 of 4,800 runs of one 16-thread check (2026-10-25, Node 24.21.0) this
        // isolate's view of the memory still ended below those after the thread trapped
        // "memory access out of bounds" at its first atomic store (`self->tid` in its TLS),
        // although plain loads of its start arguments there had passed; 2 more failed the
        // same way with the message hidden. memory.grow(0) updates this isolate's view (the
        // header of wasi_shared_memory.js): 0 of 2,800 (docs/evidence/wasm_oob_20261025/).
        memory.grow(0);
        instance.exports.wasi_thread_start(tid, startArg);
        debug(`worker done tid=${tid}`);
        if (!recycle) break;
        // Claimable when the next instance is ready; until then the host may wait for it.
        Atomics.store(view, SLOT_STATE, SLOT_PREPARING);
        parentPort.postMessage({ type: "returned", tid });
        waits = true;
      }
      completed = true;
      if (parentPort) {
        parentPort.close();
      }
    })
    .catch(failWorker);
}
