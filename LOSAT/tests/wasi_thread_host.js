const { Worker, isMainThread, parentPort, workerData } = require("worker_threads");
const { WASI } = require("wasi");
const fs = require("fs");
const crypto = require("node:crypto");
const os = require("node:os");
const path = require("node:path");
const { performance } = require("node:perf_hooks");
const { inspectArtifact, prepareArtifact, unsignedWasiImports, createSharedMemory } = require("./wasi_artifact");
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
function defaultGuardCacheDir(env = process.env, platform = process.platform, home = homeDir()) {
  let base = null;
  if (platform === "win32") base = env.LOCALAPPDATA || (home && path.join(home, "AppData", "Local"));
  else if (platform === "darwin") base = home && path.join(home, "Library", "Caches");
  else base = env.XDG_CACHE_HOME && path.isAbsolute(env.XDG_CACHE_HOME) ? env.XDG_CACHE_HOME : home && path.join(home, ".cache");
  return base ? path.join(base, "losat", "wasi-guard") : null;
}

function homeDir() {
  try { return os.homedir(); } catch { return ""; }
}

const GUARD_ENTRY = /^guarded-[0-9a-f]{64}-[0-9a-f]{64}\.json$/;
const GUARD_LEFTOVER = /^guarded-[0-9a-f]{64}-[0-9a-f]{64}\.wasm(?:\..*\.tmp)?$|^guarded-.*\.json\..*\.tmp$/;
// Older than this, bytes without a sidecar and temporary files are left over from a run that
// stopped between its writes; a younger one may belong to a run that is writing now.
const GUARD_LEFTOVER_MS = 3_600_000;
// Removes the entries (sidecar and bytes) beyond the GUARD_CACHE_KEEP used last, and leftovers.
function pruneGuardCache(dir) {
  const entries = [];
  const names = fs.readdirSync(dir);
  const sidecars = new Set(names.filter((name) => GUARD_ENTRY.test(name)));
  for (const name of names) {
    try {
      if (sidecars.has(name)) {
        entries.push({ stem: path.join(dir, name.slice(0, -5)), used: fs.statSync(path.join(dir, name)).mtimeMs });
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

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:68-93
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...) : ..., m_Outfile(outfile),
// The native command's write to standard output or error blocks while the reader drains a
// full pipe. A Worker started with Node's default stdio opens process.stdout and
// process.stderr, and Node sets a pipe there non-blocking, so the guest's fd_write to a full
// pipe returns EAGAIN (WASI errno 6). Node's uv_fs_write returns EAGAIN only when it wrote
// nothing (a partial write returns its count), so the host waits 1 ms (Atomics.wait, which
// Node allows on every thread) and repeats the same write until the pipe takes it.
const WASI_ERRNO_AGAIN = 6;
const STDIO_WAIT = new Int32Array(new SharedArrayBuffer(4));

function blockingStdioWrite(fdWrite) {
  return (fd, iovs, iovsLen, nwritten) => {
    for (;;) {
      const errno = fdWrite(fd, iovs, iovsLen, nwritten);
      if (errno !== WASI_ERRNO_AGAIN || (fd !== 1 && fd !== 2)) return errno;
      Atomics.wait(STDIO_WAIT, 0, 0, 1);
    }
  };
}

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:114-127
// ```c++
// CMemoryFile* CSeqDBAtlas::GetMemoryFile(const string& fileName)
// {
//     std::lock_guard<std::mutex> guard(m_FileMemMapMutex);
//     auto it = m_FileMemMap.find(fileName);
//     ...
//     CAtlasMappedFile* file(new CAtlasMappedFile(fileName));
//     m_FileMemMap[fileName].reset(file);
// ```
// NCBI opens a database file once for the process, and every thread reads it there. LOSAT does
// the same with one descriptor (LOSAT/src/seqdb/file.rs: a seek and a read under a lock), so a
// descriptor must name the same open file on every thread, as in the native process's table.
// Node's WASI keeps its table in the instance (uvwasi), and each thread worker has its own, so a
// descriptor that one thread opened was EBADF on the others (-db with 2 or more threads,
// 2026-10-26). In a command that may run more than one thread, the host therefore opens itself
// (Node's fs, whose descriptors belong to the process) each regular file that path_open opens
// from the preopened "/" for reading only (no open or descriptor flags, no write right), enters
// it in a table that every thread shares (the OS descriptor, the file position, the rights and
// the kind) and gives the guest SHARED_FD_BASE + its slot; each thread then reads, seeks, stats
// and closes such a descriptor itself, through the table. Every other open for reading (a
// directory, a FIFO, one that fails, with uvwasi's errno) stays with the calling thread's uvwasi
// as before.
// (One I/O worker that ran the file functions of every thread was tried first: each call waited
// until that worker was scheduled on the 16 busy CPUs, and BLASTN at 16 threads took 22.6 s
// instead of 9.3 s on ws, 2026-10-27.)
//
// NCBI reference (598d8ae6): c++/src/app/blast/blastp_app.cpp:225-228
// ```c++
//         CBlastFormat formatter(opt, *db_adapter,
//                                fmt_args->GetFormattedOutputChoice(),
//                                query_opts->GetParseDeflines(),
//                                m_CmdLineArgs->GetOutputStream(),
// ```
// The output (-out) is one stream of the process too. LOSAT creates the file at the report's
// first write, and blastp-fast's query blocks searched two at a time (LOSAT_X_BPQSLOTS,
// LOSAT/src/algorithm/blastp/x_query_slots.rs) are written by whichever pool thread finished
// the block that completes the order: the file opened by one thread was written by another,
// EBADF (blastp-fast -subject -out with 2 or more threads, 2026-10-28). Each open for writing
// therefore enters the table too: from the preopened "/", with no open flag but create,
// exclusive and truncate and no descriptor flag but append, of a regular file (a new one with
// create) or of a character device or a FIFO (/dev/null; /dev/stdout when standard output is one
// of these; a named pipe). A regular file not opened to append is written at the table's
// position, the bytes of one write reserved before they are written as one write() moves the
// native offset; an append and a device or FIFO are written at the OS's position (the end of the
// file, the stream). A device or FIFO is not seeked through the table (ESPIPE, where native
// /dev/null seeks; LOSAT does not seek an output). Every other open for writing (from another
// directory descriptor, with other flags, of a directory, a socket or a block device, a symbolic
// link not followed, one that fails) stays with the calling thread's uvwasi.
const SHARED_FD_BASE = 0x40000000, SHARED_FDS = 1024;
const WASI_EBADF = 8, WASI_EINVAL = 28, WASI_EIO = 29, WASI_ENOTSUP = 58, WASI_ESPIPE = 70, WASI_ENOTCAPABLE = 76;
const WASI_ERRNO = { EACCES: 2, EAGAIN: 6, EBADF: WASI_EBADF, EDQUOT: 19, EEXIST: 20, EFBIG: 22, EINTR: 27, EINVAL: WASI_EINVAL,
  EIO: WASI_EIO, EISDIR: 31, ELOOP: 32, ENOENT: 44, ENOSPC: 51, ENOTDIR: 54, EOVERFLOW: 61, EPERM: 63, EPIPE: 64,
  EROFS: 69, ESPIPE: WASI_ESPIPE };
const RIGHT_FD_READ = 1n << 1n, RIGHT_FD_WRITE = 1n << 6n, RIGHT_FD_ALLOCATE = 1n << 8n, RIGHT_FD_FILESTAT_SET_SIZE = 1n << 22n;
const LOOKUP_SYMLINK_FOLLOW = 1, OFLAGS_CREAT = 1, OFLAGS_EXCL = 4, OFLAGS_TRUNC = 8, FDFLAGS_APPEND = 1, PREOPEN_ROOT_FD = 3;
// The filetypes of the table's files, as uvwasi gives them (a FIFO is a socket stream there).
const FILETYPE_CHARACTER_DEVICE = 2, FILETYPE_REGULAR_FILE = 4, FILETYPE_SOCKET_STREAM = 6;
// The table: per slot the OS descriptor (-1 free, -2 being entered), the position, the rights and
// the kind (the filetype, with KIND_APPEND for a file opened to append).
const KIND_APPEND = 0x100;
const sharedFileViews = table => ({ os: new Int32Array(table, 0, SHARED_FDS),
  position: new BigInt64Array(table, 4 * SHARED_FDS, SHARED_FDS),
  rights: new BigUint64Array(table, 12 * SHARED_FDS, 2 * SHARED_FDS),
  kind: new Int32Array(table, 28 * SHARED_FDS, SHARED_FDS) });
function sharedFileTable() {
  const table = new SharedArrayBuffer(32 * SHARED_FDS);
  sharedFileViews(table).os.fill(-1);
  return table;
}

// The file functions of one thread with the shared table in front of its own uvwasi functions
// (`local`, which take the guest's signed i32 values themselves).
function sharedFileImports(table, memory, local) {
  const { os, position, rights, kind } = sharedFileViews(table);
  const slotOf = fd => (fd >>> 0) - SHARED_FD_BASE;
  const shared = fd => slotOf(fd) >= 0 && slotOf(fd) < SHARED_FDS;
  // A regular file has the table's position; it is written there unless it was opened to append.
  const regular = slot => (kind[slot] & 0xff) === FILETYPE_REGULAR_FILE;
  // The guest's memory as this isolate sees it; memory.grow(0) updates the view when the guest
  // grew the memory on another thread (see the start of a thread below).
  const buffer = end => {
    if (end > memory.buffer.byteLength) memory.grow(0);
    return memory.buffer;
  };
  const view = (pointer, length) => new DataView(buffer(pointer + length), pointer, length);
  const lengthOf = (iovs, count) => {
    let total = 0n;
    const vectors = view(iovs, 8 * count);
    for (let i = 0; i < count; i++) total += BigInt(vectors.getUint32(8 * i + 4, true));
    return total;
  };
  // fs.readSync or fs.writeSync over the guest's vectors at `at` (null: the OS's position); an
  // error after some bytes ends the call with their count, as read() and write() do. writeSync
  // takes a Number (it writes at the OS's position when given a BigInt); readSync a BigInt.
  const transfer = (move, descriptor, iovs, count, at) => {
    let total = 0n;
    const vectors = view(iovs, 8 * count);
    for (let i = 0; i < count; i++) {
      const pointer = vectors.getUint32(8 * i, true), length = vectors.getUint32(8 * i + 4, true);
      if (length === 0) continue;
      let done;
      try {
        done = move(descriptor, new Uint8Array(buffer(pointer + length), pointer, length), 0, length,
          at === null ? null : move === fs.writeSync ? Number(at + total) : at + total);
      } catch (err) {
        if (total > 0n) break;
        throw err;
      }
      total += BigInt(done);
      if (done < length) break;
    }
    return total;
  };
  const imports = {};
  // Any function given a shared descriptor where that is not one of the files: no right to it.
  for (const [name, call] of Object.entries(local)) {
    if (!/^(fd|path|sock)_/.test(name)) continue;
    const fds = { path_symlink: [2], path_rename: [0, 3], path_link: [0, 4], fd_renumber: [0, 1] }[name] || [0];
    imports[name] = (...args) => (fds.some(i => shared(args[i])) ? WASI_ENOTCAPABLE : call(...args));
  }
  const own = (name, op) => {
    const call = imports[name];
    imports[name] = (fd, ...args) => {
      if (!shared(fd)) return call(fd, ...args);
      const slot = slotOf(fd), descriptor = Atomics.load(os, slot);
      if (descriptor < 0) return WASI_EBADF;
      try {
        return op(slot, descriptor, ...args.map(arg => (typeof arg === "number" ? arg >>> 0 : arg))) ?? 0;
      } catch (err) {
        return WASI_ERRNO[err && err.code] ?? WASI_EIO;
      }
    };
  };
  const allowed = (slot, right) => (rights[2 * slot] & right) !== 0n;
  own("fd_read", (slot, descriptor, iovs, count, nread) => {
    if (!allowed(slot, RIGHT_FD_READ)) return WASI_ENOTCAPABLE;
    const got = transfer(fs.readSync, descriptor, iovs, count, regular(slot) ? Atomics.load(position, slot) : null);
    if (regular(slot)) Atomics.add(position, slot, got);
    view(nread, 4).setUint32(0, Number(got), true);
  });
  own("fd_pread", (slot, descriptor, iovs, count, offset, nread) => {
    if (!allowed(slot, RIGHT_FD_READ)) return WASI_ENOTCAPABLE;
    if (!regular(slot)) return WASI_ESPIPE;
    view(nread, 4).setUint32(0, Number(transfer(fs.readSync, descriptor, iovs, count, BigInt.asUintN(64, offset))), true);
  });
  own("fd_write", (slot, descriptor, iovs, count, nwritten) => {
    if (!allowed(slot, RIGHT_FD_WRITE)) return WASI_ENOTCAPABLE;
    let wrote;
    if (kind[slot] === FILETYPE_REGULAR_FILE) {
      const length = lengthOf(iovs, count), at = Atomics.add(position, slot, length);
      try {
        wrote = transfer(fs.writeSync, descriptor, iovs, count, at);
      } finally {
        // The bytes not written given back, unless another write came after them.
        if (wrote !== length) Atomics.compareExchange(position, slot, at + length, at + (wrote ?? 0n));
      }
    } else {
      wrote = transfer(fs.writeSync, descriptor, iovs, count, null);
      // An append leaves the position at the end of the file, as it leaves the native offset.
      if (regular(slot)) Atomics.store(position, slot, fs.fstatSync(descriptor, { bigint: true }).size);
    }
    view(nwritten, 4).setUint32(0, Number(wrote), true);
  });
  own("fd_pwrite", (slot, descriptor, iovs, count, offset, nwritten) => {
    if (!allowed(slot, RIGHT_FD_WRITE)) return WASI_ENOTCAPABLE;
    if (!regular(slot)) return WASI_ESPIPE;
    view(nwritten, 4).setUint32(0, Number(transfer(fs.writeSync, descriptor, iovs, count, BigInt.asUintN(64, offset))), true);
  });
  own("fd_seek", (slot, descriptor, offset, whence, result) => {
    if (!regular(slot)) return WASI_ESPIPE;
    const from = [0n, Atomics.load(position, slot), fs.fstatSync(descriptor, { bigint: true }).size][whence];
    if (from === undefined || from + offset < 0n) return WASI_EINVAL;
    Atomics.store(position, slot, from + offset);
    view(result, 8).setBigUint64(0, from + offset, true);
  });
  own("fd_tell", (slot, descriptor, result) => {
    if (!regular(slot)) return WASI_ESPIPE;
    view(result, 8).setBigUint64(0, Atomics.load(position, slot), true);
  });
  own("fd_close", (slot, descriptor) => {
    if (Atomics.compareExchange(os, slot, descriptor, -1) !== descriptor) return WASI_EBADF;
    fs.closeSync(descriptor);
  });
  own("fd_filestat_get", (slot, descriptor, result) => {
    const stat = fs.fstatSync(descriptor, { bigint: true }), out = view(result, 64);
    for (const [at, value] of [[0, stat.dev], [8, stat.ino], [16, 0n], [24, stat.nlink], [32, stat.size],
      [40, stat.atimeNs], [48, stat.mtimeNs], [56, stat.ctimeNs]]) out.setBigUint64(at, value, true);
    out.setUint8(16, kind[slot] & 0xff);
  });
  own("fd_fdstat_get", (slot, descriptor, result) => {
    const out = view(result, 24);
    out.setBigUint64(0, 0n, true);
    out.setUint8(0, kind[slot] & 0xff);
    out.setUint16(2, kind[slot] & KIND_APPEND ? FDFLAGS_APPEND : 0, true);
    out.setBigUint64(8, rights[2 * slot], true);
    out.setBigUint64(16, rights[2 * slot + 1], true);
  });
  own("fd_fdstat_set_rights", (slot, descriptor, base, inheriting) => {
    [base, inheriting] = [BigInt.asUintN(64, base), BigInt.asUintN(64, inheriting)];
    if ((base & ~rights[2 * slot]) || (inheriting & ~rights[2 * slot + 1])) return WASI_ENOTCAPABLE;
    [rights[2 * slot], rights[2 * slot + 1]] = [base, inheriting];
  });
  own("fd_fdstat_set_flags", (slot, descriptor, flags) => (flags === (kind[slot] & KIND_APPEND ? FDFLAGS_APPEND : 0) ? 0 : WASI_ENOTSUP));
  own("fd_filestat_set_size", (slot, descriptor, size) => {
    if (!allowed(slot, RIGHT_FD_FILESTAT_SET_SIZE)) return WASI_ENOTCAPABLE;
    fs.ftruncateSync(descriptor, Number(BigInt.asUintN(64, size)));
  });
  // Grows the file to offset + length, as uvwasi does (Node has no posix_fallocate).
  own("fd_allocate", (slot, descriptor, offset, length) => {
    if (!allowed(slot, RIGHT_FD_ALLOCATE)) return WASI_ENOTCAPABLE;
    const end = BigInt.asUintN(64, offset) + BigInt.asUintN(64, length);
    if (fs.fstatSync(descriptor, { bigint: true }).size < end) fs.ftruncateSync(descriptor, Number(end));
  });
  own("fd_advise", () => 0);
  own("fd_sync", (slot, descriptor) => { fs.fsyncSync(descriptor); });
  own("fd_datasync", (slot, descriptor) => { fs.fdatasyncSync(descriptor); });
  own("fd_prestat_get", () => WASI_EBADF);
  own("fd_prestat_dir_name", () => WASI_EBADF);
  // A free slot marked as being entered, or -1 when the table is full.
  const reserve = () => {
    for (let slot = 0; slot < SHARED_FDS; slot++) {
      if (Atomics.load(os, slot) === -1 && Atomics.compareExchange(os, slot, -1, -2) === -1) return slot;
    }
    return -1;
  };
  const pathOpen = imports.path_open;
  imports.path_open = (...args) => {
    const [dirfd, dirflags, path, pathLength, oflags] = args.slice(0, 5).map(arg => arg >>> 0);
    const [base, inheriting, fdflags, opened] = [BigInt.asUintN(64, args[5]), BigInt.asUintN(64, args[6]), args[7] >>> 0, args[8] >>> 0];
    const write = (base & RIGHT_FD_WRITE) !== 0n;
    if (dirfd !== PREOPEN_ROOT_FD || (write
      ? (oflags & ~(OFLAGS_CREAT | OFLAGS_EXCL | OFLAGS_TRUNC)) || (fdflags & ~FDFLAGS_APPEND)
      : oflags !== 0 || fdflags !== 0)) return pathOpen(...args);
    let slot = -1, descriptor = -1;
    try {
      // The guest's bytes as they are (the path may create a file).
      const host = Buffer.concat([Buffer.from("/"), new Uint8Array(buffer(path + pathLength), path, pathLength)]);
      const follow = (dirflags & LOOKUP_SYMLINK_FOLLOW) !== 0;
      // Not opened for reading unless it is a regular file now (opening a FIFO waits for its
      // writer); for writing, a regular file, a new one, a character device or a FIFO.
      const stat = (follow ? fs.statSync : fs.lstatSync)(host, { throwIfNoEntry: false });
      const kept = stat => stat.isFile() || (write && (stat.isCharacterDevice() || stat.isFIFO()));
      if (stat ? !kept(stat) : !(write && (oflags & OFLAGS_CREAT))) return pathOpen(...args);
      // A full table: the thread's own uvwasi opens it, as before (reserved before the open,
      // which may create the file).
      slot = reserve();
      if (slot < 0) return pathOpen(...args);
      const { O_RDONLY, O_WRONLY, O_RDWR, O_CREAT, O_EXCL, O_TRUNC, O_APPEND, O_NOFOLLOW } = fs.constants;
      descriptor = fs.openSync(host, (write ? (base & RIGHT_FD_READ ? O_RDWR : O_WRONLY) : O_RDONLY)
        | (oflags & OFLAGS_CREAT ? O_CREAT : 0) | (oflags & OFLAGS_EXCL ? O_EXCL : 0) | (oflags & OFLAGS_TRUNC ? O_TRUNC : 0)
        | (fdflags & FDFLAGS_APPEND ? O_APPEND : 0) | (follow ? 0 : O_NOFOLLOW));
      const now = fs.fstatSync(descriptor);
      if (!kept(now)) throw new Error("not a file of the table");
      Atomics.store(position, slot, 0n);
      [rights[2 * slot], rights[2 * slot + 1]] = [base, inheriting];
      kind[slot] = (now.isFile() ? FILETYPE_REGULAR_FILE : now.isFIFO() ? FILETYPE_SOCKET_STREAM : FILETYPE_CHARACTER_DEVICE)
        | (fdflags & FDFLAGS_APPEND ? KIND_APPEND : 0);
      Atomics.store(os, slot, descriptor);
      view(opened, 4).setUint32(0, SHARED_FD_BASE + slot, true);
      return 0;
    } catch {
      if (descriptor >= 0) fs.closeSync(descriptor);
      if (slot >= 0) Atomics.store(os, slot, -1);
      return pathOpen(...args);
    }
  };
  return imports;
}

function createImports(wasi, memory, spawnThread, extra, files) {
  const wasiImport = unsignedWasiImports(wasi.wasiImport);
  const local = { ...wasiImport, fd_write: blockingStdioWrite(wasiImport.fd_write) };
  return {
    ...extra,
    env: { memory },
    wasi: { "thread-spawn": spawnThread },
    wasi_snapshot_preview1: files ? { ...local, ...sharedFileImports(files, memory, local) } : local,
  };
}

async function instantiateWithMemory(module, memory, args, spawnThread, provided, wasi = createWasi(args), files = null) {
  const instance = await WebAssembly.instantiate(
    module,
    createImports(wasi, memory, spawnThread, hostImports(module, provided), files),
  );
  return { instance, wasi };
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
  // LOSAT_WASI_PRESPAWN: a command's standby workers start booting here, while the module is
  // read, guarded and compiled, and get the module and the memory once they exist (`early`).
  // Booting one is most of a thread's start (Node's worker bootstrap: 0.1-0.3 s for 15 workers at
  // once), so they no longer start only after the compile.
  const early = kind === "threaded-command" && !recycleWorkers(kind)
    ? bootEarlyWorkers(prespawnCount(argv), options.fault) : [];
  try {
    return await (await prepareThreadHost(wasmPath, kind)).create(argv, early.length ? { ...options, early } : options);
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
      worker = new Worker(__filename, { workerData: { control, standby: index, early: true, fault: fault || null } });
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
async function prepareThreadHost(wasmPath, kind) {
  const start = performance.now();
  const bytes = fs.readFileSync(wasmPath);
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
    module = await WebAssembly.compile(guarded.bytes);
  } catch (err) {
    inspectArtifact(bytes, kind);
    throw err;
  }
  const identity = prepareArtifact(bytes, kind, module, inputSha256).identity;
  const timings = { read_and_hash_seconds: (inspected-start)/1000, guard_seconds: (guardedAt-inspected)/1000,
    guarded_compile_seconds: (performance.now()-guardedAt)/1000 };
  debug(`shared-memory guards: ${guarded.fill} fills, ${guarded.copy} copies`);
  return { identity, timings, create: (argv = [], options = {}) =>
    instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) };
}

async function instantiatePreparedHost(wasmPath, kind, module, identity, argv, options) {
  const args = [wasmPath, ...argv];
  const prespawn = prespawnCount(argv), recycle = recycleWorkers(kind), refill = recycle ? 0 : refillCount(argv);
  // The table of the files that every thread reads (sharedFileImports), in a command that
  // requests 2 or more threads.
  const files = kind === "threaded-command" && requestedThreads(argv) >= 2 ? sharedFileTable() : null;
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
      worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
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
          worker = new Worker(__filename, { workerData: { module, wasmPath, args, memory: mainMemory,
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
module.exports = { createThreadHost, prepareThreadHost, compileThreadedModule, guardSharedMemoryCached, defaultGuardCacheDir };

if (!isMainThread) {
  const { control, fault, standby, recycle, early } = workerData;
  // An early worker (createThreadHost) boots before the module and the memory exist and gets
  // them in one message (`init`) from prespawnWorkers.
  let { module: workerModule, wasmPath, args, memory, files } = workerData;
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
      if (early) {
        ({ module: workerModule, wasmPath, args, memory, files } =
          await new Promise(resolve => parentPort.once("message", resolve)));
      }
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
        const instantiated = await instantiateWithMemory(module, memory, args, spawnThread, undefined, wasi ?? undefined, files);
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
        // NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:171-174
        // ```c++
        //     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
        //         (*thread)->Run();
        //     }
        // ```
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
