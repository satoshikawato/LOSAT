const { Worker, SHARE_ENV } = require("worker_threads");
const fs = require("fs");
const path = require("node:path");
const { performance } = require("node:perf_hooks");
const { inspectArtifact, prepareArtifact, createSharedMemory } = require("./wasi_artifact");
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Keep the parallel command's memory operations valid across V8 isolates (wasi_shared_memory.js,
// loaded when a module is guarded or checked for the build's guard here, not when the guard cache
// has the guarded bytes).
const guardSharedMemory = bytes => require("./wasi_shared_memory").guardSharedMemory(bytes);
const guardRecord = bytes => require("./wasi_shared_memory").guardRecord(bytes);
// The table of the files the command opens (wasi_files.js).
const { sharedFileTable, sharedFileViews } = require("./wasi_files");
// What every thread of the command runs on, and the file that a thread worker loads.
const {
  DEBUG, debug, writeError, START_STATE, START_ORDER, START_TID, START_ARG, SLOT_STATE, START_WORDS,
  SLOT_READY, SLOT_TAKEN, SLOT_PREPARING, instantiateWithMemory,
} = require("./wasi_thread_worker");
const WORKER_FILE = require.resolve("./wasi_thread_worker");
// Thread workers share this process's environment (SHARE_ENV) instead of each copying it at its
// creation, which this thread does for every worker (15 workers: 14 instead of 18 ms on ws). They
// only read it (createWasi).
// Their V8 heaps may grow to just under 2 GiB, not to V8's default: an isolate whose old
// generation may reach 2 GiB gets its own copy of V8's embedded builtins in its code range
// (--short-builtin-calls), 1.8 MiB of RSS per worker: a process with the 15 idle workers of a
// 16-thread command held 323 MiB on ws with no limit or 2048 MiB, 296 MiB with 2047 or 1024 MiB.
// Their young generations stay at 1 MiB, which V8 would grow from 2.25 to 4 MiB once a worker
// allocates (up to 3 MiB of RSS per worker); after its start a worker allocates only the small
// views of its file calls (sharedFileImports). At 16 threads on ws the peak RSS of 9 commands
// fell by 13 to 32 MiB with the old limit alone and by 40 to 76 MiB with both; none was slower.
// A worker's heap holds its WASI context and its instance; the guest's memory is the shared
// memory, outside it.
const WORKER_LIMITS = Object.freeze({ maxOldGenerationSizeMb: 2047, maxYoungGenerationSizeMb: 1 });

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
// The -num_threads that argv requests (its last value), or 0 when none or not a number.
function requestedThreads(argv) {
  let value;
  for (let i = 0; i < argv.length; i++) {
    const option = /^--?num_threads(?:=(.*))?$/s.exec(String(argv[i]));
    if (option) value = option[1] !== undefined ? option[1] : argv[++i];
  }
  return /^\d+$/.test(String(value)) ? Number(value) : 0;
}

function prespawnCount(argv) {
  const raw = process.env.LOSAT_WASI_PRESPAWN ?? "auto";
  if (raw === "auto") {
    const requested = requestedThreads(argv);
    return requested >= 2 && requested - 1 <= PRESPAWN_MAX ? requested - 1 : 0;
  }
  const count = /^\d+$/.test(raw) ? Number(raw) : NaN;
  if (!(count <= PRESPAWN_MAX)) {
    throw new Error(`LOSAT_WASI_PRESPAWN must be auto or an integer from 0 to ${PRESPAWN_MAX}`);
  }
  return count;
}

// LOSAT_WASI_PRESPAWN=auto only: when a spawn finds no standby worker left, the workers of the
// next set of threads boot together, as many as the -num_threads that argv requests. That is
// the most one set asks for: a pool whose threads are all new (LOSAT/src/utils/threading.rs
// x_spawn_search_pool), as blastn's query-batch slots start after the search's pool, k pools
// of m threads with k x m = -num_threads. One new worker per spawn instead, each spawn waits
// for its own worker's boot (16 one after another, about 40 ms each, for the slots of a
// 16-thread search on the workstation). Unused workers are terminated at close, as the first
// set's. 0 (a new worker per spawn, as before) for a number of workers and without standby.
function refillCount(argv) {
  if ((process.env.LOSAT_WASI_PRESPAWN ?? "auto") !== "auto" || !prespawnCount(argv)) return 0;
  return requestedThreads(argv);
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
// node:crypto is loaded at its first use: with the guard cache's index (prepareThreadHost), the
// cached module is compiling before anything is hashed.
let cryptoModule = null;
const crypto = () => cryptoModule ??= require("node:crypto");
const sha256 = bytes => crypto().createHash("sha256").update(bytes).digest("hex");
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

// The checks of a cached entry: its sidecar (`record`) names this input and guard version and
// covers its own fields, and the guarded bytes have the recorded length and SHA-256.
function validGuardEntry(record, bytesLength, bytesSha256, inputSha256) {
  if (record === null || typeof record !== "object") return false;
  const { record_sha256: recordSha256, ...fields } = record;
  const count = value => Number.isSafeInteger(value) && value >= 0;
  return recordSha256 === sha256(JSON.stringify(fields)) && fields.format === GUARD_CACHE_FORMAT
    && fields.input_sha256 === inputSha256 && fields.guard_version === guardVersion()
    && count(fields.fill) && count(fields.copy) && fields.output_bytes === bytesLength
    && fields.output_sha256 === bytesSha256
    // guardSharedMemory returns its input unchanged exactly when it guards nothing.
    && (fields.fill + fields.copy === 0) === (fields.output_sha256 === inputSha256);
}

// The sidecar and the guarded bytes of an entry, unchecked.
function readGuardEntry(stem) {
  try {
    return { record: JSON.parse(fs.readFileSync(`${stem}.json`, "utf8")), bytes: fs.readFileSync(`${stem}.wasm`) };
  } catch (err) {
    if (err.code === "ENOENT") return { cache: "miss" };
    if (err instanceof SyntaxError) return { cache: "invalid" };
    throw err;
  }
}

function readGuardCache(stem, inputSha256) {
  const entry = readGuardEntry(stem);
  if (entry.cache) return entry;
  const { record, bytes } = entry;
  return validGuardEntry(record, bytes.length, sha256(bytes), inputSha256)
    ? { cache: "hit", bytes, fill: record.fill, copy: record.copy } : { cache: "invalid" };
}

// A whole file read into the start of `into`; the view of the bytes read. A file longer than
// `into` is an error (RangeError), as one that changes while it is read may be.
function readInto(file, into) {
  const fd = fs.openSync(file, "r");
  try {
    const size = fs.fstatSync(fd).size;
    if (size > into.length) throw new RangeError(`${file} is larger than its buffer`);
    let done = 0;
    while (done < size) {
      const read = fs.readSync(fd, into, done, size - done, done);
      if (read === 0) throw new RangeError(`${file} ended while it was read`);
      done += read;
    }
    return into.subarray(0, size);
  } finally {
    fs.closeSync(fd);
  }
}

// Frees the memory of `view` at once, not at a garbage collection, which the main thread of a
// command may not run before it exits. Only a view of a whole ArrayBuffer is freed (never a part
// of Node's Buffer pool); the ArrayBuffer is detached (ArrayBuffer.prototype.transfer), and a
// later read of the view sees 0 bytes.
function freeBytes(view) {
  const buffer = view.buffer;
  if (view.byteOffset === 0 && view.byteLength === buffer.byteLength && buffer.byteLength > 0
    && typeof buffer.transfer === "function") buffer.transfer(0);
}

function writeAtomic(file, data) {
  const temp = `${file}.${process.pid}.${crypto().randomBytes(6).toString("hex")}.tmp`;
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
function defaultGuardCacheDir(env = process.env, platform = process.platform, home = homeDir()) {
  let base = null;
  if (platform === "win32") base = env.LOCALAPPDATA || (home && path.join(home, "AppData", "Local"));
  else if (platform === "darwin") base = home && path.join(home, "Library", "Caches");
  else base = env.XDG_CACHE_HOME && path.isAbsolute(env.XDG_CACHE_HOME) ? env.XDG_CACHE_HOME : home && path.join(home, ".cache");
  return base ? path.join(base, "losat", "wasi-guard") : null;
}

function homeDir() {
  try { return require("node:os").homedir(); } catch { return ""; }
}

const GUARD_ENTRY = /^guarded-[0-9a-f]{64}-[0-9a-f]{64}\.json$/;
const GUARD_LEFTOVER = /^guarded-[0-9a-f]{64}-[0-9a-f]{64}\.wasm(?:\..*\.tmp)?$|^(?:guarded|index)-.*\.json\..*\.tmp$/;
// Older than this, bytes without a sidecar and temporary files are left over from a run that
// stopped between its writes; a younger one may belong to a run that is writing now.
const GUARD_LEFTOVER_MS = 3_600_000;
// Removes the entries (sidecar and bytes) beyond the GUARD_CACHE_KEEP used last, the index files
// beyond the GUARD_CACHE_KEEP used last, and leftovers.
function pruneGuardCache(dir) {
  const entries = [], indexes = [];
  const names = fs.readdirSync(dir);
  const sidecars = new Set(names.filter((name) => GUARD_ENTRY.test(name)));
  for (const name of names) {
    try {
      if (sidecars.has(name)) {
        entries.push({ stem: path.join(dir, name.slice(0, -5)), used: fs.statSync(path.join(dir, name)).mtimeMs });
      } else if (GUARD_INDEX.test(name)) {
        indexes.push({ file: path.join(dir, name), used: fs.statSync(path.join(dir, name)).mtimeMs });
      } else if (GUARD_LEFTOVER.test(name) && !sidecars.has(name.replace(/\.wasm$/, ".json"))) {
        const file = path.join(dir, name);
        if (Date.now() - fs.statSync(file).mtimeMs > GUARD_LEFTOVER_MS) fs.rmSync(file, { force: true });
      }
    } catch { /* removed by another run meanwhile */ }
  }
  entries.sort((a, b) => b.used - a.used);
  for (const { stem } of entries.slice(GUARD_CACHE_KEEP)) {
    // The sidecar first: bytes without a sidecar are never used.
    fs.rmSync(`${stem}.json`, { force: true });
    fs.rmSync(`${stem}.wasm`, { force: true });
  }
  indexes.sort((a, b) => b.used - a.used);
  for (const { file } of indexes.slice(GUARD_CACHE_KEEP)) fs.rmSync(file, { force: true });
}

// The cache's directory, or null when the cache is off (LOSAT_WASI_GUARD_CACHE=0, or no default).
function guardCacheDir() {
  const named = process.env.LOSAT_WASI_GUARD_CACHE;
  if (named === "0") return null;
  if (named === "") throw new Error("LOSAT_WASI_GUARD_CACHE must name a directory or be 0");
  return named ?? defaultGuardCacheDir();
}

function guardSharedMemoryCached(bytes, inputSha256 = sha256(bytes)) {
  const named = process.env.LOSAT_WASI_GUARD_CACHE;
  const dir = guardCacheDir();
  // An artifact guarded at build time (`wasi_shared_memory.js guard`) runs as it is: there is
  // nothing to compute, and nothing is written to the cache.
  if (guardRecord(bytes) !== null) return { ...guardSharedMemory(bytes), cache: "guarded" };
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
// Host diagnostics only (LOSAT_WASI_THREADS_DEBUG): the guard counts of this run, or of the build.
function guardSummary(guarded) {
  if (guarded.fill + guarded.copy === 0 && guarded.record) {
    return `shared-memory guards: ${guarded.record.fill} fills, ${guarded.record.copy} copies (guarded at build time)`;
  }
  return `shared-memory guards: ${guarded.fill} fills, ${guarded.copy} copies`;
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host-only (LOSAT_WASI_GUARD_CACHE): the cache's index. An entry is named by the SHA-256 of the
// input, so finding it reads and hashes the whole input before the compile can start, and the
// hit then hashes the cached bytes too (the 9-10 MB command: about 40 ms on ws, where the compile
// takes about 30). An index file named by the identity of the input file (device, inode, size,
// and the modification and change times in ns) records the input's SHA-256 and the guard version
// that an earlier run found. The cached guarded bytes are then compiled at once, and every check
// of a hit (readGuardCache: the input's SHA-256, from the input read again, the guard version,
// the sidecar, the length and SHA-256 of the cached bytes) is made while V8 compiles, before the
// module is used. A check that fails discards that module, and the input is read and guarded as
// without the index. A write or a rename of the input usually changes its change time and so its
// index name. Where the timestamps are coarse the name can stay the same after a rewrite, and
// then the SHA-256 check of the input rejects the index. The index saves only the time before
// the compile, never a check.
const GUARD_INDEX = /^index-\d+-\d+-\d+-\d+-\d+\.json$/;
const HEX64 = /^[0-9a-f]{64}$/;
const guardIndexFile = (dir, stat) =>
  path.join(dir, `index-${stat.dev}-${stat.ino}-${stat.size}-${stat.mtimeNs}-${stat.ctimeNs}.json`);

// After a run that found or wrote the entry of this input (`stat`, taken after its read).
function writeGuardIndex(stat, inputSha256) {
  try {
    const dir = guardCacheDir();
    if (dir === null || stat === null) return;
    const index = { format: GUARD_CACHE_FORMAT, input_sha256: inputSha256, guard_version: guardVersion() };
    writeAtomic(guardIndexFile(dir, stat), `${JSON.stringify(index)}\n`);
    // Each new copy or build of a module adds an index file, also when its entry exists.
    pruneGuardCache(dir);
  } catch (err) {
    debug(`guard cache index not written (${err.message})`);
  }
}

// The buffer that compileIndexedGuard reads the module into is freed once the compile has copied
// the bytes (V8 keeps its own copy with the module) and the identity is read; it used to stay to
// the end of a command (10 MB of RSS). It is at least 33 MiB, of which only the module's pages
// are touched: glibc's free() of a mapped block of up to 32 MiB raises the mmap threshold to that
// size, and later blocks below it come from the arenas and stay there once freed. On ws (16
// threads on 4 CPUs, 2 runs each) the peak RSS of a blastx, blastp, megablast and tblastx command
// rose by 7, 12, 22 and 30 MiB with a 10 MB buffer freed, and fell by 9, 11 and 13 MiB (tblastx:
// within its spread) with the buffer at 33 MiB.
const INDEX_BUFFER_BYTES = 33 << 20;

// The module compiled from the cached guarded bytes that the input's index names, with the input,
// the buffer it was read into (INDEX_BUFFER_BYTES) and its SHA-256, or null when there is no
// index, no entry, or a check fails. `compileStarted` runs once the compile has started
// (createThreadHost's early workers).
async function compileIndexedGuard(wasmPath, compileStarted) {
  let dir, stat, index, stem, record, buffer, cached;
  try {
    dir = guardCacheDir();
    if (dir === null) return null;
    stat = fs.statSync(wasmPath, { bigint: true });
    index = JSON.parse(fs.readFileSync(guardIndexFile(dir, stat), "utf8"));
    if (!index || index.format !== GUARD_CACHE_FORMAT || !HEX64.test(index.input_sha256)
      || !HEX64.test(index.guard_version)) return null;
    stem = path.join(dir, `guarded-${index.input_sha256}-${index.guard_version}`);
    record = JSON.parse(fs.readFileSync(`${stem}.json`, "utf8"));
    // The guard replaces instructions and adds a few functions: its output is never twice the input.
    if (!record || !Number.isSafeInteger(record.output_bytes) || record.output_bytes < 0
      || record.output_bytes > 2 * Number(stat.size) + 65536) return null;
    // One buffer for the cached bytes and then the input: V8 copies the bytes when the compile
    // starts, and a new 10 MB buffer costs its page faults (about 6 ms on ws).
    buffer = Buffer.allocUnsafeSlow(Math.max(record.output_bytes, Number(stat.size), INDEX_BUFFER_BYTES));
    cached = readInto(`${stem}.wasm`, buffer);
  } catch (err) {
    debug(`guard cache index not used (${err.message})`);
    return null;
  }
  const compiled = WebAssembly.compile(cached);
  // Bytes that do not compile fail the checks below like any other stale entry.
  compiled.catch(() => {});
  compileStarted();
  let input, inputSha256;
  try {
    const cachedLength = cached.length, cachedSha256 = sha256(cached);
    input = readInto(wasmPath, buffer);
    inputSha256 = sha256(input);
    if (inputSha256 !== index.input_sha256 || index.guard_version !== guardVersion()
      || !validGuardEntry(record, cachedLength, cachedSha256, inputSha256)) {
      debug(`guard cache index stale: ${stem}.wasm`);
      return null;
    }
  } catch (err) {
    debug(`guard cache index not used (${err.message})`);
    return null;
  }
  let module;
  try {
    module = await compiled;
  } catch (err) {
    debug(`guard cache index: cached bytes did not compile (${err.message})`);
    return null;
  }
  const now = new Date();
  try {
    fs.utimesSync(`${stem}.json`, now, now);
    fs.utimesSync(guardIndexFile(dir, stat), now, now);
  } catch { /* only the order of pruning */ }
  debug(`guard cache index hit: ${stem}.wasm`);
  return { module, buffer, input, inputSha256, fill: record.fill, copy: record.copy };
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Apply the same shared-memory compatibility guard before main or worker
// compilation. Workers normally receive this already compiled module.
async function compileThreadedModule(wasmPath) {
  const guarded = guardSharedMemory(fs.readFileSync(wasmPath));
  debug(guardSummary(guarded));
  return WebAssembly.compile(guarded.bytes);
}

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

async function instantiateMain(module, args, spawnThread, limits, provided, files) {
  const initial = parsePositiveIntEnv("LOSAT_WASM_MEMORY_INITIAL_PAGES", limits.initial);
  // A maximum given in the environment is used as it is, without stepping down.
  const requested = parsePositiveIntEnv("LOSAT_WASM_MEMORY_MAXIMUM_PAGES", limits.maximum);
  const { memory, maximum } = createSharedMemory(initial, requested,
    process.env.LOSAT_WASM_MEMORY_MAXIMUM_PAGES === undefined);
  debug(`shared memory: initial ${initial} pages, maximum ${maximum} pages (module declares ${limits.maximum})`);
  return { ...await instantiateWithMemory(module, memory, args, spawnThread, provided, undefined, files), memory };
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
  // LOSAT_WASI_PRESPAWN: a command's standby workers start booting while the module is compiled,
  // and get the module and the memory once they exist (`early`). Booting one is most of a
  // thread's start (Node's worker bootstrap: 0.1-0.3 s for 15 workers at once), so they do not
  // start only after the compile. They are created once the compile has started: with the guard
  // cache's index right after the cached guarded bytes are read and handed to V8
  // (compileIndexedGuard; the input read and both hashes then run after the workers exist), else
  // after the module is read and guarded. The reads of the module come before the workers on the
  // path without the index: a read of the 10 MB module right after 15 workers were created took
  // 5-15 ms instead of 2 (ws node 0), and with the workers created first the guest started at 121
  // instead of 78 ms.
  const count = kind === "threaded-command" && !recycleWorkers(kind) ? prespawnCount(argv) : 0;
  let early = [];
  try {
    const prepared = await prepareThreadHost(wasmPath, kind, () => { early = bootEarlyWorkers(count, options.fault); });
    return await prepared.create(argv, early.length ? { ...options, early } : options);
  } catch (err) {
    await terminateWorkers(new Set(early.map(slot => slot.worker)));
    throw err;
  }
}

// LOSAT_WASI_PRESPAWN (createThreadHost): `count` standby workers created before the module is
// compiled; each boots and waits for its `init` message from prespawnWorkers. Fewer when a worker
// cannot be created (prespawnWorkers creates the rest as before).
function bootEarlyWorkers(count, fault) {
  const early = [];
  for (let index = 0; index < count; index++) {
    const control = new SharedArrayBuffer(4 * START_WORDS);
    let worker;
    try {
      worker = new Worker(WORKER_FILE, { env: SHARE_ENV, resourceLimits: WORKER_LIMITS, workerData: { control, standby: index, early: true, fault: fault || null } });
    } catch { break; }
    worker.unref();
    early.push({ worker, control });
  }
  return early;
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// (*thread)->Run(); (*thread)->Join(&result);
// Explicit module reuse for host measurements: each create still owns a fresh
// instance/memory/worker set. There is no global cache or retained search pool
// (LOSAT_WASI_GUARD_CACHE keeps only guarded bytes on disk); a reactor's
// thread workers stay with its host between searches (LOSAT_WASI_RECYCLE) until close().
async function prepareThreadHost(wasmPath, kind, compiling = () => {}) {
  const start = performance.now();
  // createThreadHost boots its early workers at the first compile that starts, while V8 compiles.
  let started = false;
  const compileStarted = () => { if (!started) { started = true; compiling(); } };
  const indexed = await compileIndexedGuard(wasmPath, compileStarted);
  if (indexed) {
    const identity = prepareArtifact(indexed.input, kind, indexed.module, indexed.inputSha256).identity;
    freeBytes(indexed.buffer);
    const timings = { read_and_hash_seconds: 0, guard_seconds: 0, guarded_compile_seconds: (performance.now()-start)/1000 };
    debug(`shared-memory guards: ${indexed.fill} fills, ${indexed.copy} copies`);
    const { module } = indexed;
    return { identity, timings, create: (argv = [], options = {}) =>
      instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) };
  }
  const bytes = fs.readFileSync(wasmPath);
  let stat = null;
  try { stat = fs.statSync(wasmPath, { bigint: true }); } catch { /* no index for this run */ }
  // One compile: the identity (imports, exports) is read from the guarded module, whose imports
  // and exports are the input's (the guard rewrites only memory.fill and memory.copy), instead of
  // compiling the input too (about 40 ms for the 9 MB command). An input that cannot be guarded
  // or compiled is inspected first, so that an artifact of another kind is refused as before.
  const inputSha256 = sha256(bytes);
  const inspected = performance.now();
  let guarded, module, guardedAt;
  try {
    guarded = guardSharedMemoryCached(bytes, inputSha256);
    guardedAt = performance.now();
    const compiled = WebAssembly.compile(guarded.bytes);
    compileStarted();
    module = await compiled;
  } catch (err) {
    inspectArtifact(bytes, kind);
    throw err;
  }
  const identity = prepareArtifact(bytes, kind, module, inputSha256).identity;
  // An artifact guarded at build time has no cache entry for an index to name.
  if (guarded.cache !== "off" && guarded.cache !== "guarded") writeGuardIndex(stat, inputSha256);
  const timings = { read_and_hash_seconds: (inspected-start)/1000, guard_seconds: (guardedAt-inspected)/1000,
    guarded_compile_seconds: (performance.now()-guardedAt)/1000 };
  debug(guardSummary(guarded));
  return { identity, timings, create: (argv = [], options = {}) =>
    instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) };
}

async function instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) {
  const args = [wasmPath, ...argv];
  const prespawn = prespawnCount(argv), recycle = recycleWorkers(kind), refill = recycle ? 0 : refillCount(argv);
  // The table of the files that every thread reads and writes (sharedFileImports), in a command
  // at every thread count (a reactor opens no files).
  const files = kind === "threaded-command" ? sharedFileTable() : null;
  const workers = new Set(), exits = new Map(), events = [], standby = [];
  // LOSAT_WASI_RECYCLE: the workers that run one guest thread after another, the guest threads
  // they run now (tids), and how each of those threads' exit promise is resolved.
  const recycled = [], running = new Set(), returned = new Map();
  let nextTid = 1, attempt = 0, nextSlot = 0, nextStandby = 0, mainMemory, closing = false, failure = null;
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
  // startArg. The guest gets the tid after the start confirmation only from a worker that still
  // boots (the paragraph below).
  // A worker that already waits for its order (SLOT_READY) has instantiated the module and only
  // takes the order and starts the thread: the spawn returns the tid without waiting for its start
  // confirmation, as pthread_create returns before the new thread runs (each confirmation is a
  // wake-up of that worker, ~0.5 ms, one after another for the 15 threads of a 16-thread pool on
  // ws). From here a failure of that thread ends the command as it does after a confirmation (the
  // worker's failWorker). A worker that still boots or instantiates is waited for as before, so a
  // worker that never starts still ends the command with the start timeout.
  const startStandby = (slot, tid, startArg) => {
    const view = new Int32Array(slot.control);
    const waiting = Atomics.compareExchange(view, SLOT_STATE, SLOT_READY, SLOT_TAKEN) === SLOT_READY;
    slot.tid = tid; slot.worker.ref();
    workers.add(slot.worker); exits.set(slot.worker, slot.exited); event("spawned", { tid, standby: slot.index });
    Atomics.store(view, START_TID, tid); Atomics.store(view, START_ARG, startArg);
    Atomics.store(view, START_ORDER, 1); Atomics.notify(view, START_ORDER, 1);
    if (waiting) { event("ready", { tid, waiting: true }); return tid; }
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
    const worker = new Worker(WORKER_FILE, { env: SHARE_ENV, resourceLimits: WORKER_LIMITS, workerData: { module, wasmPath, args, memory: mainMemory,
      control, files, fault: options.fault || null, recycle: true, ...data } });
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
    // refillCount: the next set boots together; this spawn takes the first of it.
    if (!standby.length && refill) { event("standby_refill", { tid, count: refill }); prespawnWorkers(refill); }
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
      worker = new Worker(WORKER_FILE, { env: SHARE_ENV, resourceLimits: WORKER_LIMITS, workerData: { module, wasmPath, args, memory: mainMemory,
        startArg, tid, control, files, fault: options.fault || null } });
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
  // hands them out in order and, when none is left, boots the next set (refillCount) or creates
  // a new worker per spawn. A standby worker never keeps the process alive (unref until taken),
  // and close() terminates the unused ones.
  const prespawnWorkers = (count) => {
    for (let n = 0; n < count; n++) {
      const index = nextStandby++;
      if (recycle) {
        // A recycled worker that waits from the start; spawnThread takes it like the others.
        try { recycledWorker({ standby: index }).worker.unref(); }
        catch (err) { event("standby_rejected", { standby: index, error: String(err) }); return; }
        event("standby_spawned", { standby: index });
        continue;
      }
      // An early worker (createThreadHost) has booted meanwhile and gets what it instantiates now.
      const early = options.early ? options.early[index] : undefined;
      const control = early ? early.control : new SharedArrayBuffer(4 * START_WORDS);
      let worker;
      try {
        if (early) {
          worker = early.worker;
          worker.postMessage({ module, wasmPath, args, memory: mainMemory, files });
        } else {
          worker = new Worker(WORKER_FILE, { env: SHARE_ENV, resourceLimits: WORKER_LIMITS, workerData: { module, wasmPath, args, memory: mainMemory,
            control, files, standby: index, fault: options.fault || null } });
        }
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
  const instantiated = await instantiateMain(module, args, spawnThread, identity.memory, options.imports, files);
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
      // The shared files that the guest left open.
      if (files) for (const descriptor of sharedFileViews(files).os) if (descriptor >= 0) fs.closeSync(descriptor);
      exits.clear();
    },
  };
}
module.exports = { createThreadHost, prepareThreadHost, compileThreadedModule, guardSharedMemoryCached, defaultGuardCacheDir, freeBytes };
