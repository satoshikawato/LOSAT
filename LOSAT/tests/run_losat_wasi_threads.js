const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const { WASI } = require("wasi");
const fs = require("fs");
const os = require("os");

const DEFAULT_INITIAL_MEMORY_PAGES = 21;
const DEFAULT_MAXIMUM_MEMORY_PAGES = 16384;
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
  if (raw === undefined || raw === "") {
    return fallback;
  }
  const parsed = Number.parseInt(raw, 10);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    throw new Error(`${name} must be a positive integer`);
  }
  return parsed;
}

// NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// ```c
// const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
// int num_threads = args[kArgNumThreads].AsInteger();
// if (num_threads > kMaxValue) {
//     m_NumThreads = kMaxValue;
// } else {
//     m_NumThreads = num_threads;
// }
// ```
function defaultRunnerThreadCap() {
  if (typeof os.availableParallelism === "function") {
    return Math.max(1, os.availableParallelism());
  }
  return Math.max(1, os.cpus().length || 1);
}

function debug(message) {
  if (DEBUG) {
    fs.writeSync(2, `[losat-wasi-threads] ${message}\n`);
  }
}

function writeError(err) {
  const message = err && err.stack ? err.stack : String(err);
  fs.writeSync(2, `${message}\n`);
}

function markWorkerStart(control, state) {
  const view = new Int32Array(control);
  Atomics.store(view, 0, state);
  Atomics.notify(view, 0, 1);
}

function findWasiSymbol(wasi, name) {
  return Object.getOwnPropertySymbols(wasi).find((symbol) => symbol.description === name);
}

function initializeWasiThreadInstance(wasi, instance) {
  try {
    wasi.initialize(instance);
    return;
  } catch (err) {
    if (!String(err && err.message).includes("_start")) {
      throw err;
    }
  }

  const setMemorySymbol = findWasiSymbol(wasi, "kSetMemory");
  const startedSymbol = findWasiSymbol(wasi, "kStarted");
  const instanceSymbol = findWasiSymbol(wasi, "kInstance");
  if (!setMemorySymbol || !startedSymbol || !instanceSymbol) {
    throw new Error("Node WASI internals are unavailable; cannot initialize a WASI thread instance");
  }

  wasi[setMemorySymbol](instance.exports.memory);
  wasi[instanceSymbol] = instance;
  wasi[startedSymbol] = true;
}

function createWasi(args) {
  return new WASI({
    version: "preview1",
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

async function instantiateMain(module, args, spawnThread) {
  let initial = parsePositiveIntEnv(
    "LOSAT_WASM_MEMORY_INITIAL_PAGES",
    DEFAULT_INITIAL_MEMORY_PAGES,
  );
  let maximum = parsePositiveIntEnv(
    "LOSAT_WASM_MEMORY_MAXIMUM_PAGES",
    DEFAULT_MAXIMUM_MEMORY_PAGES,
  );

  for (let attempt = 0; attempt < 3; attempt += 1) {
    const memory = new WebAssembly.Memory({ initial, maximum, shared: true });
    try {
      const instantiated = await instantiateWithMemory(module, memory, args, spawnThread);
      return { ...instantiated, memory };
    } catch (err) {
      const message = err && err.message ? err.message : String(err);
      const initialMatch = message.match(/smaller than initial (\d+)/);
      if (initialMatch) {
        initial = Number.parseInt(initialMatch[1], 10);
        continue;
      }
      const maximumMatch = message.match(/larger maximum size \d+ than the module's declared maximum (\d+)/);
      if (maximumMatch) {
        maximum = Number.parseInt(maximumMatch[1], 10);
        continue;
      }
      throw err;
    }
  }

  throw new Error("failed to instantiate threaded Wasm memory after limit retries");
}

async function terminateWorkers(workers) {
  const pending = [];
  for (const worker of workers) {
    pending.push(worker.terminate());
  }
  workers.clear();
  await Promise.allSettled(pending);
}

if (isMainThread) {
  const wasmPath = process.argv[2];
  if (!wasmPath) {
    console.error("usage: node tests/run_losat_wasi_threads.js <LOSAT.wasm> [args...]");
    process.exit(2);
  }

  const args = [wasmPath, ...process.argv.slice(3)];
  const bytes = fs.readFileSync(wasmPath);
  const workers = new Set();
  let nextTid = 1;
  let terminatingWorkers = false;
  // NCBI reference: ncbi-blast/c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
  // ```c
  // const int kMaxValue = static_cast<int>(CSystemInfo::GetCpuCount());
  // int num_threads = args[kArgNumThreads].AsInteger();
  // if (num_threads > kMaxValue) {
  //     m_NumThreads = kMaxValue;
  // } else {
  //     m_NumThreads = num_threads;
  // }
  // ```
  const runnerThreadCap = parsePositiveIntEnv(
    "LOSAT_WASI_THREAD_CAP",
    defaultRunnerThreadCap(),
  );
  process.env.LOSAT_WASI_THREAD_CAP = String(runnerThreadCap);
  debug(`runner thread cap=${runnerThreadCap}`);

  WebAssembly.compile(bytes)
    .then(async (module) => {
      const spawnThread = (startArg) => {
        const tid = nextTid;
        nextTid += 1;
        debug(`spawn tid=${tid} start_arg=${startArg}`);
        const control = new SharedArrayBuffer(Int32Array.BYTES_PER_ELEMENT);
        const controlView = new Int32Array(control);
        let worker;
        try {
          worker = new Worker(__filename, {
            workerData: {
              module,
              wasmPath,
              args,
              memory: mainMemory,
              startArg,
              tid,
              control,
            },
          });
        } catch (err) {
          writeError(err);
          return -1;
        }
        workers.add(worker);
        worker.on("exit", (code) => {
          debug(`worker exit code=${code}`);
          workers.delete(worker);
          if (!terminatingWorkers && code !== 0 && process.exitCode === undefined) {
            process.exitCode = code;
          }
        });
        worker.on("error", (err) => {
          writeError(err);
          if (process.exitCode === undefined) {
            process.exitCode = 1;
          }
        });
        const wait = Atomics.wait(controlView, 0, 0, WORKER_START_TIMEOUT_MS);
        if (wait === "timed-out") {
          console.error(`timed out waiting for WASI thread ${tid} to start`);
          worker.terminate();
          return -1;
        }
        if (Atomics.load(controlView, 0) !== 1) {
          return -1;
        }
        return tid;
      };

      let mainMemory = null;
      const instantiated = await instantiateMain(module, args, spawnThread);
      mainMemory = instantiated.memory;
      const { instance, wasi } = instantiated;
      if (typeof instance.exports._start !== "function") {
        throw new Error(`${wasmPath} does not export _start; use the LOSAT command Wasm artifact.`);
      }

      try {
        debug("starting main _start");
        wasi.start(instance);
      } finally {
        debug("main _start returned; terminating workers");
        terminatingWorkers = true;
        await terminateWorkers(workers);
      }
    })
    .catch(async (err) => {
      terminatingWorkers = true;
      await terminateWorkers(workers);
      writeError(err);
      process.exitCode = 1;
    });
} else {
  const { module: workerModule, wasmPath, args, memory, startArg, tid, control } = workerData;

  Promise.resolve()
    .then(async () => {
      debug(`worker boot tid=${tid}`);
      const module = workerModule || await WebAssembly.compile(fs.readFileSync(wasmPath));
      const spawnThread = () => -1;
      debug(`worker instantiate tid=${tid} start_arg=${startArg}`);
      const { instance, wasi } = await instantiateWithMemory(module, memory, args, spawnThread);
      initializeWasiThreadInstance(wasi, instance);
      if (typeof instance.exports.wasi_thread_start !== "function") {
        throw new Error(`${wasmPath} does not export wasi_thread_start`);
      }
      debug(`worker start tid=${tid}`);
      markWorkerStart(control, 1);
      instance.exports.wasi_thread_start(tid, startArg);
      debug(`worker done tid=${tid}`);
      if (parentPort) {
        parentPort.close();
      }
    })
    .catch((err) => {
      writeError(err);
      markWorkerStart(control, -1);
      process.exitCode = 1;
    });
}
