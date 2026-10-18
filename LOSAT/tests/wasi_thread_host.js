const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const { WASI } = require("wasi");
const fs = require("fs");
const crypto = require("node:crypto");
const os = require("node:os");
const path = require("node:path");
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

// NCBI reference (598d8ae6): c++/src/algo/blast/blastinput/blast_args.cpp:3207-3222
// ```c++
//     if (args.Exist(kArgNumThreads) &&
//         args[kArgNumThreads].HasValue()) {  // could be cancelled by the exclusion in CRemoteArgs
//         ...
//         int num_threads = args[kArgNumThreads].AsInteger();
//         ...
//             m_NumThreads = num_threads;
// ```
// Host-only (LOSAT_WASI_PRESPAWN; `auto` when unset since 2026-10-18, `0` boots none): the
// number of thread workers to boot ahead of `_start`. `auto` reads the requested
// -num_threads from argv (the caller is one of the N threads, so N-1 workers); a value the
// guest will reject gives 0. The limit is the
// guest's: rayon-core 1.13.0 sleep/counters.rs:57-77 (THREADS_BITS = 8 on 32-bit targets,
// THREADS_MAX = 255) bounds the pool that LOSAT/src/utils/threading.rs validate_threads
// accepts, so at most 254 workers. Booting ahead changes when a worker boots, not which
// workers the guest gets or what they run.
const PRESPAWN_MAX = 254;
function prespawnCount(argv) {
  const raw = process.env.LOSAT_WASI_PRESPAWN ?? "auto";
  if (raw === "auto") {
    let value;
    for (let i = 0; i < argv.length; i++) {
      const option = /^--?num_threads(?:=(.*))?$/s.exec(String(argv[i]));
      if (option) value = option[1] !== undefined ? option[1] : argv[++i];
    }
    const requested = /^\d+$/.test(String(value)) ? Number(value) : 0;
    return requested >= 2 && requested - 1 <= PRESPAWN_MAX ? requested - 1 : 0;
  }
  const count = /^\d+$/.test(raw) ? Number(raw) : NaN;
  if (!(count <= PRESPAWN_MAX)) {
    throw new Error(`LOSAT_WASI_PRESPAWN must be auto or an integer from 0 to ${PRESPAWN_MAX}`);
  }
  return count;
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host-only (LOSAT_WASI_GUARD_CACHE; on by default since 2026-10-19, `0` turns it off): reuse
// the guarded bytes of an earlier run. The key is the SHA-256 of the input bytes and of the
// guard's own source, so a changed artifact or guard never reuses an entry. A hit is used
// only when its sidecar names this key, the sidecar's own SHA-256 covers its fields (the
// fill / copy counts included), and the cached bytes have the recorded SHA-256 and length;
// anything else is recomputed by guardSharedMemory and replaced atomically (temporary file,
// then rename). The guarded bytes themselves are never changed here. Unset, the directory is
// the user's cache directory (defaultGuardCacheDir), and a failure to use it only means that
// the guard is computed as without a cache; a directory that the variable names must work.
// A write keeps the GUARD_CACHE_KEEP entries used last (a hit marks its entry used).
const GUARD_CACHE_FORMAT = 1;
const GUARD_CACHE_KEEP = 8;
const sha256 = bytes => crypto.createHash("sha256").update(bytes).digest("hex");
let guardVersionKey = null;
function guardVersion() {
  if (guardVersionKey === null) guardVersionKey = sha256(fs.readFileSync(require.resolve("./wasi_shared_memory")));
  return guardVersionKey;
}

function guardCacheRecord(inputSha256, guarded) {
  const fields = { format: GUARD_CACHE_FORMAT, input_sha256: inputSha256, guard_version: guardVersion(),
    output_sha256: sha256(guarded.bytes), output_bytes: guarded.bytes.length, fill: guarded.fill, copy: guarded.copy };
  return { ...fields, record_sha256: sha256(JSON.stringify(fields)) };
}

function readGuardCache(stem, inputSha256) {
  let record, bytes;
  try {
    record = JSON.parse(fs.readFileSync(`${stem}.json`, "utf8"));
    bytes = fs.readFileSync(`${stem}.wasm`);
  } catch (err) {
    if (err.code === "ENOENT") return { cache: "miss" };
    if (err instanceof SyntaxError) return { cache: "invalid" };
    throw err;
  }
  if (record === null || typeof record !== "object") return { cache: "invalid" };
  const { record_sha256: recordSha256, ...fields } = record;
  const count = value => Number.isSafeInteger(value) && value >= 0;
  const valid = recordSha256 === sha256(JSON.stringify(fields)) && fields.format === GUARD_CACHE_FORMAT
    && fields.input_sha256 === inputSha256 && fields.guard_version === guardVersion()
    && count(fields.fill) && count(fields.copy) && fields.output_bytes === bytes.length
    && fields.output_sha256 === sha256(bytes)
    // guardSharedMemory returns its input unchanged exactly when it guards nothing.
    && (fields.fill + fields.copy === 0) === (fields.output_sha256 === inputSha256);
  return valid ? { cache: "hit", bytes, fill: fields.fill, copy: fields.copy } : { cache: "invalid" };
}

function writeAtomic(file, data) {
  const temp = `${file}.${process.pid}.${crypto.randomBytes(6).toString("hex")}.tmp`;
  try {
    fs.writeFileSync(temp, data, { flag: "wx" });
    fs.renameSync(temp, file);
  } catch (err) {
    fs.rmSync(temp, { force: true });
    throw err;
  }
}

// The user's cache directory: $XDG_CACHE_HOME (absolute) or ~/.cache, ~/Library/Caches on
// macOS, %LOCALAPPDATA% on Windows; null without a home directory.
function defaultGuardCacheDir(env = process.env, platform = process.platform) {
  let home;
  try { home = os.homedir(); } catch { home = ""; }
  let base = null;
  if (platform === "win32") base = env.LOCALAPPDATA || (home && path.join(home, "AppData", "Local"));
  else if (platform === "darwin") base = home && path.join(home, "Library", "Caches");
  else base = env.XDG_CACHE_HOME && path.isAbsolute(env.XDG_CACHE_HOME) ? env.XDG_CACHE_HOME : home && path.join(home, ".cache");
  return base ? path.join(base, "losat", "wasi-guard") : null;
}

const GUARD_ENTRY = /^guarded-[0-9a-f]{64}-[0-9a-f]{64}\.json$/;
// Removes the entries (sidecar and bytes) beyond the GUARD_CACHE_KEEP used last.
function pruneGuardCache(dir) {
  const entries = [];
  for (const name of fs.readdirSync(dir)) {
    if (!GUARD_ENTRY.test(name)) continue;
    try {
      entries.push({ stem: path.join(dir, name.slice(0, -5)), used: fs.statSync(path.join(dir, name)).mtimeMs });
    } catch { /* removed by another run meanwhile */ }
  }
  entries.sort((a, b) => b.used - a.used);
  for (const { stem } of entries.slice(GUARD_CACHE_KEEP)) {
    // The sidecar first: bytes without a sidecar are never used.
    fs.rmSync(`${stem}.json`, { force: true });
    fs.rmSync(`${stem}.wasm`, { force: true });
  }
}

function guardSharedMemoryCached(bytes, inputSha256 = sha256(bytes)) {
  const named = process.env.LOSAT_WASI_GUARD_CACHE;
  if (named === "0") return { ...guardSharedMemory(bytes), cache: "off" };
  if (named === "") throw new Error("LOSAT_WASI_GUARD_CACHE must name a directory or be 0");
  const dir = named ?? defaultGuardCacheDir();
  if (dir === null) return { ...guardSharedMemory(bytes), cache: "off" };
  // The default directory may be unusable (no home, read-only, a file in the way).
  const unusable = (err) => {
    if (named !== undefined) throw err;
    debug(`guard cache unusable (${err.message}): ${dir}`);
  };
  const stem = path.join(dir, `guarded-${inputSha256}-${guardVersion()}`);
  let cached;
  try {
    cached = readGuardCache(stem, inputSha256);
  } catch (err) {
    unusable(err);
    return { ...guardSharedMemory(bytes), cache: "off" };
  }
  debug(`guard cache ${cached.cache}: ${stem}.wasm`);
  if (cached.cache === "hit") {
    const now = new Date();
    try { fs.utimesSync(`${stem}.json`, now, now); } catch { /* only the order of pruning */ }
    return cached;
  }
  const guarded = guardSharedMemory(bytes);
  try {
    fs.mkdirSync(dir, { recursive: true });
    // The sidecar is written last; a reader that pairs it with other bytes rejects the pair.
    writeAtomic(`${stem}.wasm`, guarded.bytes);
    writeAtomic(`${stem}.json`, `${JSON.stringify(guardCacheRecord(inputSha256, guarded))}\n`);
  } catch (err) {
    unusable(err);
    return { ...guarded, cache: "off" };
  }
  try {
    pruneGuardCache(dir);
  } catch (err) {
    debug(`guard cache not pruned (${err.message}): ${dir}`);
  }
  return { ...guarded, cache: cached.cache };
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

// Control words shared with a thread worker: the start confirmation (every worker); for a
// worker that waits for its thread (a LOSAT_WASI_PRESPAWN standby worker, or a reactor's
// worker whose earlier thread returned, LOSAT_WASI_RECYCLE), the start order with the tid and
// startArg; and for a recycled worker, whether it waits for a thread (SLOT_STATE).
const START_STATE = 0, START_ORDER = 1, START_TID = 2, START_ARG = 3, SLOT_STATE = 4, START_WORDS = 5;
// SLOT_STATE: a thread runs (or the worker boots for its first thread), the worker waits for a
// thread, the host gave it one, or its thread returned and it instantiates the module again.
const SLOT_RUNNING = 0, SLOT_READY = 1, SLOT_TAKEN = 2, SLOT_PREPARING = 3;

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
// ```c++
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     }
// ```
// Host-only (LOSAT_WASI_RECYCLE; `1` when unset, `0` turns it off): in a reactor, a thread
// worker outlives its guest thread. When the thread returns, the worker instantiates the
// module again on the shared memory (a new instance per guest thread, as before) and waits
// for a start order, so a later search's `thread-spawn` takes it instead of booting a new
// Node worker (about 20 ms each, one after another while the guest waits). This is what the
// browser's thread host (web/app/src/infra/engine-worker/thread-host.ts) does. A command
// runs one search, so its workers still exit with their threads. The guest gets the same
// threads in the same order; only the worker that runs each one changes.
function recycleWorkers(kind) {
  const raw = process.env.LOSAT_WASI_RECYCLE ?? "1";
  if (raw !== "0" && raw !== "1") throw new Error("LOSAT_WASI_RECYCLE must be 0 or 1");
  return kind === "threaded-reactor" && raw === "1";
}

function markWorkerStart(control, state) {
  const view = new Int32Array(control);
  Atomics.store(view, START_STATE, state);
  Atomics.notify(view, START_STATE, 1);
}

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
// ```c++
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     }
// ```
// A standby or recycled worker blocks here, after instantiation, until spawnThread assigns
// it a thread. The order is consumed, so a recycled worker waits again after this thread.
function waitForStartOrder(control) {
  const view = new Int32Array(control);
  while (Atomics.load(view, START_ORDER) === 0) Atomics.wait(view, START_ORDER, 0);
  const order = { tid: Atomics.load(view, START_TID), startArg: Atomics.load(view, START_ARG) };
  Atomics.store(view, START_ORDER, 0);
  return order;
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
    env: { ...process.env, PWD: process.cwd() },
    // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729
// ```c++
//             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
// ```
// Preserve lexical filenames at the recorded host cwd.
// Root preopen permits absolute paths; PWD initializes the guest cwd for
// relative filenames through std::env::set_current_dir, without rewriting argv.
    preopens: { "/": "/" },
  });
}

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:68-93
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

function createImports(wasi, memory, spawnThread, extra) {
  return {
    ...extra,
    env: { memory },
    wasi: { "thread-spawn": spawnThread },
    wasi_snapshot_preview1: wasi.wasiImport,
  };
}

async function instantiateWithMemory(module, memory, args, spawnThread, provided, wasi = createWasi(args)) {
  const instance = await WebAssembly.instantiate(
    module,
    createImports(wasi, memory, spawnThread, hostImports(module, provided)),
  );
  return { instance, wasi };
}

async function instantiateMain(module, args, spawnThread, limits, provided) {
  const initial = parsePositiveIntEnv("LOSAT_WASM_MEMORY_INITIAL_PAGES", limits.initial);
  const maximum = parsePositiveIntEnv("LOSAT_WASM_MEMORY_MAXIMUM_PAGES", limits.maximum);
  const memory = new WebAssembly.Memory({ initial, maximum, shared: true });
  return { ...await instantiateWithMemory(module, memory, args, spawnThread, provided), memory };
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
// instance/memory/worker set. There is no global cache or retained search pool
// (LOSAT_WASI_GUARD_CACHE keeps only guarded bytes on disk); a reactor's
// thread workers stay with its host between searches (LOSAT_WASI_RECYCLE) until close().
async function prepareThreadHost(wasmPath, kind) {
  const start = performance.now();
  const bytes = fs.readFileSync(wasmPath);
  const identity = inspectArtifact(bytes, kind);
  const inspected = performance.now();
  const guarded = guardSharedMemoryCached(bytes, identity.sha256);
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
  const prespawn = prespawnCount(argv), recycle = recycleWorkers(kind);
  const workers = new Set(), exits = new Map(), events = [], standby = [];
  // LOSAT_WASI_RECYCLE: the workers that run one guest thread after another, the guest threads
  // they run now (tids), and how each of those threads' exit promise is resolved.
  const recycled = [], running = new Set(), returned = new Map();
  let nextTid = 1, attempt = 0, nextSlot = 0, mainMemory, closing = false, failure = null;
  const event = (type, fields = {}) => {
    const record = { event: type, at_ms: performance.now(), ...fields }; events.push(record);
    if (DEBUG) fs.writeSync(2, `[losat-wasi-event] ${JSON.stringify(record)}\n`);
  };
  const threadExited = (worker, tid, code) => {
    workers.delete(worker); event("exited", { tid, code });
    if (!closing && code !== 0) failure = new Error(`worker ${tid} exited with code ${code}`);
  };
  const awaitThreadStart = (view, tid) => {
    const wait = Atomics.wait(view, START_STATE, 0, options.startTimeoutMs || WORKER_START_TIMEOUT_MS);
    if (wait === "timed-out" || Atomics.load(view, START_STATE) !== 1) {
      // The guest has not received a successful tid and may release startArg.
      // Abort the host rather than allowing a delayed worker to use freed state.
      writeError(new Error(`timed out waiting for WASI thread ${tid} to start`));
      process.kill(process.pid, "SIGTERM");
      return -1;
    }
    event("ready", { tid }); return tid;
  };
  // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
  // ```c++
  //     // ... launch the threads ...
  //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
  //         (*thread)->Run();
  //     }
  // ```
  // A standby worker becomes the guest's thread only here: it receives this tid and
  // startArg, and the guest gets the tid after the same start confirmation as a new worker.
  const startStandby = (slot, tid, startArg) => {
    const view = new Int32Array(slot.control);
    slot.tid = tid; slot.worker.ref();
    workers.add(slot.worker); exits.set(slot.worker, slot.exited); event("spawned", { tid, standby: slot.index });
    Atomics.store(view, START_TID, tid); Atomics.store(view, START_ARG, startArg);
    Atomics.store(view, START_ORDER, 1); Atomics.notify(view, START_ORDER, 1);
    return awaitThreadStart(view, tid);
  };
  // LOSAT_WASI_RECYCLE: a recycled worker's guest thread ends with the worker's `returned`
  // message (the worker lives on). While the guest runs, the main thread cannot receive
  // messages, so whether a worker waits for a thread is read from its control words.
  const threadStarted = (tid) => {
    running.add(tid); exits.set(tid, new Promise(resolve => returned.set(tid, resolve)));
  };
  const threadEnded = (tid, code) => {
    if (!running.delete(tid)) return;
    event("exited", { tid, code }); returned.get(tid)(); returned.delete(tid);
  };
  const recycledWorker = (data) => {
    const control = new SharedArrayBuffer(4 * START_WORDS);
    const worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
      control, fault: options.fault || null, recycle: true, ...data } });
    const slot = { index: nextSlot++, worker, control, state: new Int32Array(control), tid: data.tid ?? null };
    if (data.tid === undefined) Atomics.store(slot.state, SLOT_STATE, SLOT_PREPARING);
    worker.on("message", message => {
      if (!message || message.type !== "returned") return;
      threadEnded(message.tid, 0);
      // An idle worker never keeps the process alive; spawnThread refs it again.
      if (Atomics.load(slot.state, SLOT_STATE) !== SLOT_TAKEN) worker.unref();
    });
    worker.on("exit", code => {
      const at = recycled.indexOf(slot);
      if (at >= 0) recycled.splice(at, 1);
      if (closing) return;
      failure = new Error(`WASI worker ${slot.tid ?? `slot ${slot.index}`} exited with code ${code}`);
      if (slot.tid !== null) threadEnded(slot.tid, code);
    });
    worker.on("error", err => { failure = err; writeError(err); });
    recycled.push(slot);
    return slot;
  };
  // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
  // ```c++
  //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
  //         (*thread)->Run();
  //     }
  // ```
  // A recycled worker that waits receives this tid and startArg, and the guest gets the tid
  // after the same start confirmation as from a new worker.
  // A worker that waits, else one whose thread has returned and that instantiates the module
  // again (a few ms; a new worker takes about 20 ms to boot). A worker whose thread still runs
  // is not waited for: it may be a thread of the previous pool, which returns only after this
  // spawn (the guest holds a lock across thread-spawn that an ending thread needs).
  const claimRecycled = () => {
    const deadline = performance.now() + (options.startTimeoutMs || WORKER_START_TIMEOUT_MS);
    for (;;) {
      const waiting = recycled.find(slot =>
        Atomics.compareExchange(slot.state, SLOT_STATE, SLOT_READY, SLOT_TAKEN) === SLOT_READY);
      if (waiting) return waiting;
      const preparing = recycled.find(slot => Atomics.load(slot.state, SLOT_STATE) === SLOT_PREPARING);
      if (!preparing || performance.now() > deadline) return undefined;
      Atomics.wait(preparing.state, SLOT_STATE, SLOT_PREPARING, 50);
    }
  };
  const startRecycled = (slot, tid, startArg) => {
    const view = slot.state;
    Atomics.store(view, START_STATE, 0);
    slot.tid = tid; slot.worker.ref(); threadStarted(tid); event("spawned", { tid, recycled: slot.index });
    Atomics.store(view, START_TID, tid); Atomics.store(view, START_ARG, startArg);
    Atomics.store(view, START_ORDER, 1); Atomics.notify(view, START_ORDER, 1);
    return awaitThreadStart(view, tid);
  };
  const spawnThread = (startArg) => {
    attempt++; const tid = nextTid++; event("spawn_attempt", { tid, attempt });
    if (options.rejectSpawn && options.rejectSpawn(attempt)) { event("spawn_rejected", { tid }); return -1; }
    if (standby.length) return startStandby(standby.shift(), tid, startArg);
    if (recycle) {
      const waiting = claimRecycled();
      if (waiting) return startRecycled(waiting, tid, startArg);
      // None waits or is about to (the first search, or the threads of the previous pool have
      // not returned yet): boot a worker for this thread, which then stays.
      let slot;
      try { slot = recycledWorker({ tid, startArg }); } catch (err) { event("spawn_rejected", { tid, error: String(err) }); return -1; }
      threadStarted(tid); event("spawned", { tid });
      return awaitThreadStart(slot.state, tid);
    }
    const control = new SharedArrayBuffer(4 * START_WORDS), view = new Int32Array(control);
    let worker;
    try {
      worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
        startArg, tid, control, fault: options.fault || null } });
    } catch (err) { event("spawn_rejected", { tid, error: String(err) }); return -1; }
    workers.add(worker); event("spawned", { tid });
    exits.set(worker, new Promise(resolve => {
      worker.on("exit", code => { threadExited(worker, tid, code); resolve(); });
    }));
    worker.on("error", err => { failure = err; writeError(err); });
    return awaitThreadStart(view, tid);
  };
  // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:156-163
  // ```c++
  //     // Create the threads ...
  //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
  //         thread->Reset(new CPrelimSearchThread(internal_data,
  //                                               opts_memento.get()));
  // ```
  // LOSAT_WASI_PRESPAWN: boot the workers together, before the guest asks for them. Each
  // instantiates the module with the shared memory and waits for its start order; spawnThread
  // hands them out in order and creates new workers when none is left. A standby worker never
  // keeps the process alive (unref until taken), and close() terminates the unused ones.
  const prespawnWorkers = (count) => {
    for (let index = 0; index < count; index++) {
      if (recycle) {
        // A recycled worker that waits from the start; spawnThread takes it like the others.
        try { recycledWorker({ standby: index }).worker.unref(); }
        catch (err) { event("standby_rejected", { standby: index, error: String(err) }); return; }
        event("standby_spawned", { standby: index });
        continue;
      }
      const control = new SharedArrayBuffer(4 * START_WORDS);
      let worker;
      try {
        worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
          control, standby: index, fault: options.fault || null } });
      } catch (err) { event("standby_rejected", { standby: index, error: String(err) }); return; }
      worker.unref();
      const slot = { index, worker, control, tid: null };
      slot.exited = new Promise(resolve => {
        worker.on("exit", code => {
          if (slot.tid !== null) threadExited(worker, slot.tid, code);
          else {
            const at = standby.indexOf(slot);
            if (at >= 0) standby.splice(at, 1);
            event("standby_exited", { standby: index, code });
            if (!closing) failure = new Error(`standby WASI worker ${index} exited with code ${code}`);
          }
          resolve();
        });
      });
      worker.on("error", err => { failure = err; writeError(err); });
      standby.push(slot); event("standby_spawned", { standby: index });
    }
  };
  const instantiated = await instantiateMain(module, args, spawnThread, identity.memory, options.imports);
  mainMemory = instantiated.memory;
  const { instance, wasi } = instantiated;
  if (kind === "threaded-reactor") wasi.initialize(instance);
  if (prespawn) prespawnWorkers(prespawn);
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
        if (workers.size !== 0 || running.size !== 0) throw new Error("unreaped WASI workers");
        // Return this search's evidence to its owner instead of retaining a
        // growing history on a repeatedly used reactor host.
        return events.splice(0);
      } finally { clearTimeout(timer); }
    },
    async close() {
      closing = true;
      for (const slot of standby) event("standby_unused", { standby: slot.index });
      await Promise.all([terminateWorkers(workers), terminateWorkers(new Set(standby.splice(0).map(slot => slot.worker))),
        terminateWorkers(new Set(recycled.splice(0).map(slot => slot.worker)))]);
      exits.clear();
    },
  };
}
module.exports = { createThreadHost, prepareThreadHost, compileThreadedModule, guardSharedMemoryCached, defaultGuardCacheDir };

if (!isMainThread) {
  const { module: workerModule, wasmPath, args, memory, control, fault, standby, recycle } = workerData;
  // A standby worker (LOSAT_WASI_PRESPAWN) learns its tid and startArg from its start order,
  // and so does a recycled worker (LOSAT_WASI_RECYCLE) for every thread after its first.
  let { startArg, tid } = workerData;
  const name = standby === undefined ? `tid=${tid}` : `standby=${standby}`;
  // NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
  // LOG_POST(Error << "Error: " << e.what());
  // exit_code = BLAST_UNKNOWN_ERROR;
  process.on("uncaughtException", failWorker);
  // NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
  // LOG_POST(Error << "Error: " << e.what()); exit_code = BLAST_UNKNOWN_ERROR;
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
      const module = workerModule || await compileThreadedModule(wasmPath);
      const spawnThread = () => -1;
      const view = new Int32Array(control);
      let waits = standby !== undefined, wasi = null;
      for (;;) {
        debug(waits ? `worker instantiate ${name}` : `worker instantiate tid=${tid} start_arg=${startArg}`);
        // A new instance on the shared memory for every guest thread, recycled worker or not. A
        // recycled worker keeps its one WASI context (bound once to the same shared memory): a
        // new context per thread would hold one more open preopen descriptor for each thread.
        const created = wasi === null;
        const instantiated = await instantiateWithMemory(module, memory, args, spawnThread, undefined, wasi ?? undefined);
        const { instance } = instantiated;
        wasi = instantiated.wasi;
        if (created) initializeWasiThreadInstance(wasi, instance);
        if (typeof instance.exports.wasi_thread_start !== "function") {
          throw new Error(`${wasmPath} does not export wasi_thread_start`);
        }
        if (waits) {
          if (recycle) { Atomics.store(view, SLOT_STATE, SLOT_READY); Atomics.notify(view, SLOT_STATE); }
          debug(`worker standby ${name}`);
          ({ tid, startArg } = waitForStartOrder(control));
          debug(`worker ${name} takes tid=${tid} start_arg=${startArg}`);
        }
        debug(`worker start tid=${tid}`);
        markWorkerStart(control, 1);
        if (fault === "trap") throw new WebAssembly.RuntimeError("injected worker trap");
        if (fault === "abnormal-exit") process.exit(23);
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
