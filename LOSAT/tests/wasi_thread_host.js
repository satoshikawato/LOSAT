const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const { WASI } = require("wasi");
const fs = require("fs");
const { performance } = require("node:perf_hooks");
const { inspectArtifact } = require("./wasi_artifact");
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Keep the parallel command's memory operations valid across V8 isolates.
const { guardSharedMemory } = require("./wasi_shared_memory");

const DEBUG = process.env.LOSAT_WASI_THREADS_DEBUG === "1";
const WORKER_START_TIMEOUT_MS = parsePositiveIntEnv(
  "LOSAT_WASI_THREADS_START_TIMEOUT_MS",
  30000,
);

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// ```c
// int num_threads = args[kArgNumThreads].AsInteger();
// if (num_threads > kMaxValue) {
//     m_NumThreads = kMaxValue;
// } else {
//     m_NumThreads = num_threads;
// }
// ```
//
// This runner provides the WASI `thread-spawn` import needed by Rust's
// wasm32-wasip1-threads target so LOSAT can execute its NCBI-style
// num_threads paths instead of silently falling back to serial command-WASI.

function parsePositiveIntEnv(name, fallback) {
  const raw = process.env[name];
  if (raw === undefined) {
    return fallback;
  }
  const parsed = /^\d+$/.test(raw) ? Number(raw) : NaN;
  if (!Number.isSafeInteger(parsed) || parsed <= 0) {
    throw new Error(`${name} must be a positive integer`);
  }
  return parsed;
}

function debug(message) {
  if (DEBUG) {
    fs.writeSync(2, `[losat-wasi-threads] ${message}\n`);
  }
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Apply the same shared-memory compatibility guard before main or worker
// compilation. Workers normally receive this already compiled module.
async function compileThreadedModule(wasmPath) {
  const guarded = guardSharedMemory(fs.readFileSync(wasmPath));
  debug(`shared-memory guards: ${guarded.fill} fills, ${guarded.copy} copies`);
  return WebAssembly.compile(guarded.bytes);
}

function writeError(err) {
  const message = err && err.stack ? err.stack : String(err);
  fs.writeSync(2, `${message}\n`);
}

// NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
// catch (const std::exception& e) {
//     LOG_POST(Error << "Error: " << e.what());
//     exit_code = BLAST_UNKNOWN_ERROR;
// }
// A Wasm trap cannot unwind Rust/Rayon's shared state. The main thread may be
// blocked inside Wasm, unable to receive a Worker error/exit event. Terminate
// this command process from the failing worker after flushing its diagnostic;
// process.exit() in a Node Worker would terminate only that worker.
function failWorker(err) {
  writeError(err);
  process.kill(process.pid, "SIGTERM");
}

function markWorkerStart(control, state) {
  const view = new Int32Array(control);
  Atomics.store(view, 0, state);
  Atomics.notify(view, 0, 1);
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
// Node 18/24 lib/wasi.js initialize binds memory and optionally calls _initialize.
// Child entry has its own thread-pointer initialization; do not rerun main CRT.
function initializeWasiThreadInstance(wasi, instance) {
  wasi.initialize({ exports: { memory: instance.exports.memory } });
}

function createWasi(args) {
  return new WASI({
    version: "preview1",
    // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    // CATCH_ALL(status)
    // return status;
    // Handle the command's exit explicitly on every supported Node version.
    returnOnExit: true,
    args,
    env: process.env,
    preopens: { "/": "/" },
  });
}

function createImports(wasi, memory, spawnThread) {
  return {
    env: { memory },
    wasi: { "thread-spawn": spawnThread },
    wasi_snapshot_preview1: wasi.wasiImport,
  };
}

async function instantiateWithMemory(module, memory, args, spawnThread) {
  const wasi = createWasi(args);
  const instance = await WebAssembly.instantiate(
    module,
    createImports(wasi, memory, spawnThread),
  );
  return { instance, wasi };
}

async function instantiateMain(module, args, spawnThread, limits) {
  const initial = parsePositiveIntEnv("LOSAT_WASM_MEMORY_INITIAL_PAGES", limits.initial);
  const maximum = parsePositiveIntEnv("LOSAT_WASM_MEMORY_MAXIMUM_PAGES", limits.maximum);
  const memory = new WebAssembly.Memory({ initial, maximum, shared: true });
  return { ...await instantiateWithMemory(module, memory, args, spawnThread), memory };
}

async function terminateWorkers(workers) {
  const pending = [];
  for (const worker of workers) {
    pending.push(worker.terminate());
  }
  workers.clear();
  await Promise.allSettled(pending);
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Run(); (*thread)->Join(&result);
// Host-owned lifecycle records count attempts, ready confirmations, and real exit
// events separately. Test-only failure injection never changes search behavior.
async function createThreadHost(wasmPath, kind, argv = [], options = {}) {
  return (await prepareThreadHost(wasmPath, kind)).create(argv, options);
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// (*thread)->Run(); (*thread)->Join(&result);
// Explicit module reuse for host measurements: each create still owns a fresh
// instance/memory/worker set. There is no global cache or retained search pool.
async function prepareThreadHost(wasmPath, kind) {
  const start = performance.now();
  const bytes = fs.readFileSync(wasmPath);
  const identity = inspectArtifact(bytes, kind);
  const inspected = performance.now();
  const guarded = guardSharedMemory(bytes);
  const guardedAt = performance.now();
  const module = await WebAssembly.compile(guarded.bytes);
  const timings = { validation_including_raw_compile_seconds: (inspected-start)/1000, guard_seconds: (guardedAt-inspected)/1000,
    guarded_compile_seconds: (performance.now()-guardedAt)/1000 };
  debug(`shared-memory guards: ${guarded.fill} fills, ${guarded.copy} copies`);
  return { identity, timings, create: (argv = [], options = {}) =>
    instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) };
}

async function instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) {
  const args = [wasmPath, ...argv];
  const workers = new Set(), exits = new Map(), events = [];
  let nextTid = 1, attempt = 0, mainMemory, closing = false, failure = null;
  const event = (type, fields = {}) => {
    const record = { event: type, at_ms: performance.now(), ...fields }; events.push(record);
    if (DEBUG) fs.writeSync(2, `[losat-wasi-event] ${JSON.stringify(record)}\n`);
  };
  const spawnThread = (startArg) => {
    attempt++; const tid = nextTid++; event("spawn_attempt", { tid, attempt });
    if (options.rejectSpawn && options.rejectSpawn(attempt)) { event("spawn_rejected", { tid }); return -1; }
    const control = new SharedArrayBuffer(4), view = new Int32Array(control);
    let worker;
    try {
      worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
        startArg, tid, control, fault: options.fault || null } });
    } catch (err) { event("spawn_rejected", { tid, error: String(err) }); return -1; }
    workers.add(worker); event("spawned", { tid });
    exits.set(worker, new Promise(resolve => {
      worker.on("exit", code => {
        workers.delete(worker); event("exited", { tid, code });
        if (!closing && code !== 0) failure = new Error(`worker ${tid} exited with code ${code}`);
        resolve();
      });
    }));
    worker.on("error", err => { failure = err; writeError(err); });
    const wait = Atomics.wait(view, 0, 0, options.startTimeoutMs || WORKER_START_TIMEOUT_MS);
    if (wait === "timed-out" || Atomics.load(view, 0) !== 1) {
      // The guest has not received a successful tid and may release startArg.
      // Abort the host rather than allowing a delayed worker to use freed state.
      writeError(new Error(`timed out waiting for WASI thread ${tid} to start`));
      process.kill(process.pid, "SIGTERM");
      return -1;
    }
    event("ready", { tid }); return tid;
  };
  const instantiated = await instantiateMain(module, args, spawnThread, identity.memory);
  mainMemory = instantiated.memory;
  const { instance, wasi } = instantiated;
  if (kind === "threaded-reactor") wasi.initialize(instance);
  return {
    instance, memory: mainMemory, module, identity, events,
    start: () => wasi.start(instance),
    setRejectSpawn: predicate => { options.rejectSpawn = predicate; },
    async waitForWorkers(timeoutMs = 30000) {
      let timer;
      try {
        await Promise.race([Promise.all([...exits.values()]), new Promise((_, reject) => {
          timer = setTimeout(() => reject(new Error("WASI worker exit timeout")), timeoutMs);
        })]);
        exits.clear();
        if (failure) throw failure;
        if (workers.size !== 0) throw new Error("unreaped WASI workers");
        // Return this search's evidence to its owner instead of retaining a
        // growing history on a repeatedly used reactor host.
        return events.splice(0);
      } finally { clearTimeout(timer); }
    },
    async close() { closing = true; await terminateWorkers(workers); exits.clear(); },
  };
}
module.exports = { createThreadHost, prepareThreadHost, compileThreadedModule };

if (!isMainThread) {
  const { module: workerModule, wasmPath, args, memory, startArg, tid, control, fault } = workerData;
  // NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
  // LOG_POST(Error << "Error: " << e.what());
  // exit_code = BLAST_UNKNOWN_ERROR;
  process.on("uncaughtException", failWorker);
  // NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
  // LOG_POST(Error << "Error: " << e.what()); exit_code = BLAST_UNKNOWN_ERROR;
  // This handler runs in the worker even while the caller is blocked in Wasm.
  let completed = false;
  process.on("exit", code => {
    if (!completed) failWorker(new Error(`WASI worker ${tid} exited unexpectedly with code ${code}`));
  });

  Promise.resolve()
    .then(async () => {
      debug(`worker boot tid=${tid}`);
      if (fault === "startup-timeout") Atomics.wait(new Int32Array(new SharedArrayBuffer(4)), 0, 0);
      if (fault === "startup-error") throw new Error("injected worker startup error");
      const module = workerModule || await compileThreadedModule(wasmPath);
      const spawnThread = () => -1;
      debug(`worker instantiate tid=${tid} start_arg=${startArg}`);
      const { instance, wasi } = await instantiateWithMemory(module, memory, args, spawnThread);
      initializeWasiThreadInstance(wasi, instance);
      if (typeof instance.exports.wasi_thread_start !== "function") {
        throw new Error(`${wasmPath} does not export wasi_thread_start`);
      }
      debug(`worker start tid=${tid}`);
      markWorkerStart(control, 1);
      if (fault === "trap") throw new WebAssembly.RuntimeError("injected worker trap");
      if (fault === "abnormal-exit") process.exit(23);
      instance.exports.wasi_thread_start(tid, startArg);
      completed = true;
      debug(`worker done tid=${tid}`);
      if (parentPort) {
        parentPort.close();
      }
    })
    .catch(failWorker);
}
