"use strict";

const assert = require("node:assert/strict");
const { spawnSync } = require("node:child_process");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { test } = require("node:test");

// The runs of this file keep guarded bytes in a cache of their own, not in the user's cache
// directory (the host's default); the tests of the default name its environment themselves.
const GUARD_CACHE = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-guard-tests-"));
process.env.LOSAT_WASI_GUARD_CACHE = GUARD_CACHE;
process.on("exit", () => fs.rmSync(GUARD_CACHE, { recursive: true, force: true }));

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status)
// return status;
// NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
// LOG_POST(Error << "Error: " << e.what());
// exit_code = BLAST_UNKNOWN_ERROR;
// Exercise the host failure boundary with tiny real Wasm commands, including
// a worker trap while its parent waits in Wasm. No BLAST behavior is simulated.
// Encode fixtures directly so this test needs only Node, not a Wasm compiler.
function uleb(value) {
  const bytes = [];
  do {
    const byte = value & 127;
    value >>>= 7;
    bytes.push(byte | (value ? 128 : 0));
  } while (value);
  return bytes;
}

function i32(value) {
  const bytes = [0x41];
  while (true) {
    const byte = value & 127;
    value >>= 7;
    const done = (value === 0 && !(byte & 64)) || (value === -1 && (byte & 64));
    bytes.push(byte | (done ? 0 : 128));
    if (done) return bytes;
  }
}

const vector = (items) => [...uleb(items.length), ...items.flat()];
const name = (value) => [...uleb(Buffer.byteLength(value)), ...Buffer.from(value)];
const section = (id, bytes) => [id, ...uleb(bytes.length), ...bytes];
const body = (instructions, locals = 0) => {
  const bytes = [...(locals ? [1, locals, 0x7f] : [0]), ...instructions, 0x0b];
  return [...uleb(bytes.length), ...bytes];
};

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Entry/export/memory variations exercise Wasm ABI inspection only. `writes` (threaded only)
// appends that many fd_write calls of 4,096 bytes to standard output after the worker part,
// exiting with the first nonzero WASI errno.
function command({ threaded, status = 0, worker = null, entry = "_start",
  aliases = [], exportMemory = true, memoryLimits = [1, 0, 1], writes = 0 }) {
  const types = vector([
    [0x60, 1, 0x7f, 0],                 // proc_exit(i32)
    [0x60, 0, 0],                       // _start()
    [0x60, 1, 0x7f, 1, 0x7f],          // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],          // wasi_thread_start(i32, i32)
    [0x60, 4, 0x7f, 0x7f, 0x7f, 0x7f, 1, 0x7f], // fd_write(fd, iovs, iovs_len, nwritten) -> errno
  ]);
  const imports = [[...name("wasi_snapshot_preview1"), ...name("proc_exit"), 0, 0]];
  if (threaded) {
    imports.push([...name("env"), ...name("memory"), 2, 3, 21, ...uleb(16384)]);
    imports.push([...name("wasi"), ...name("thread-spawn"), 0, 2]);
  }
  if (writes) imports.push([...name("wasi_snapshot_preview1"), ...name("fd_write"), 0, 4]);
  // Imported functions come first in the function index space (the memory is not one).
  const imported = (threaded ? 2 : 1) + (writes ? 1 : 0);
  const exports = exportMemory ? [[...name("memory"), 2, 0]] : [];
  for (const exported of [entry, ...aliases].filter(Boolean)) {
    exports.push([...name(exported), 0, imported]);
  }
  const functions = [1];
  let main = [...i32(status), 0x10, 0];
  const bodies = [];
  if (worker) {
    // Spawn, then block in memory.atomic.wait32. A trapped worker cannot
    // release this wait, so a callback-only error handler hangs this command.
    main = [...i32(0), 0x10, 1, 0x1a,
      ...i32(0), ...i32(0), 0x42, 0x7f, 0xfe, 1, 2, 0, 0x1a];
    functions.push(3);
    if (worker !== "missing-export") exports.push([...name("wasi_thread_start"), 0, imported + 1]);
    bodies.push(body(worker === "trap"
      ? [...i32(-1), ...i32(0), ...i32(1), 0xfc, 11, 0] // OOB memory.fill
      : [...i32(0), ...i32(1), 0xfe, 0x17, 2, 0,      // atomic store + notify
        ...i32(0), ...i32(1), 0xfe, 0, 2, 0, 0x1a]));
  }
  if (threaded && !worker) {
    functions.push(3);
    exports.push([...name("wasi_thread_start"), 0, imported + 1]);
    bodies.push(body([]));
  }
  if (writes) {
    // The iovec at 16 names the 4,096 bytes at 4,096; nwritten goes to 24. Local 0 holds
    // the errno, local 1 counts the writes.
    main = [...main, ...i32(16), ...i32(4096), 0x36, 2, 0, ...i32(20), ...i32(4096), 0x36, 2, 0,
      0x03, 0x40,
      ...i32(1), ...i32(16), ...i32(1), ...i32(24), 0x10, imported - 1, 0x22, 0,
      0x04, 0x40, 0x20, 0, 0x10, 0, 0x0b,
      0x20, 1, ...i32(1), 0x6a, 0x22, 1, ...i32(writes), 0x49, 0x0d, 0,
      0x0b];
  }
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector(functions.map((type) => [type]))),
    ...(threaded ? [] : section(5, memoryLimits)),
    ...section(7, vector(exports)),
    ...section(10, vector([body(main, writes ? 2 : 0), ...bodies])),
  ]);
}

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Host memory boundary only: one WASI call whose pointer argument is above 2 GiB.
// clock_time_get writes the time to 0x80000008; the command exits with the call's
// errno, or with 99 if no time arrived there. The serial memory starts at 2 GiB and
// one page; the threaded one declares LOSAT's 65,535 pages and grows to that size.
function highPointerCommand(threaded) {
  const pages = 32769;
  const types = vector([
    [0x60, 1, 0x7f, 0],                    // proc_exit(i32)
    [0x60, 0, 0],                          // _start()
    [0x60, 1, 0x7f, 1, 0x7f],              // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],              // wasi_thread_start(i32, i32)
    [0x60, 3, 0x7f, 0x7e, 0x7f, 1, 0x7f],  // clock_time_get(i32, i64, i32) -> i32
  ]);
  const imports = [
    [...name("wasi_snapshot_preview1"), ...name("proc_exit"), 0, 0],
    [...name("wasi_snapshot_preview1"), ...name("clock_time_get"), 0, 4],
  ];
  if (threaded) {
    imports.push([...name("env"), ...name("memory"), 2, 3, 1, ...uleb(65535)]);
    imports.push([...name("wasi"), ...name("thread-spawn"), 0, 2]);
  }
  const entry = imports.length - (threaded ? 1 : 0);
  const address = i32(0x80000008 | 0);
  const start = [1, 1, 0x7f,                                  // (local i32)
    ...(threaded ? [...i32(pages - 1), 0x40, 0, 0x1a] : []),  // memory.grow; drop
    ...i32(1), 0x42, 1, ...address, 0x10, 1, 0x22, 0,         // errno = clock_time_get(monotonic, 1, address)
    0x04, 0x40, 0x20, 0, 0x10, 0, 0x0b,                       // if (errno) proc_exit(errno)
    ...address, 0x29, 3, 0, 0x50,                             // i64.load(address) == 0
    0x04, 0x40, ...i32(99), 0x10, 0, 0x0b,                    // if so: proc_exit(99)
    0x0b];
  const exports = [[...name("memory"), 2, 0], [...name("_start"), 0, entry]];
  if (threaded) exports.push([...name("wasi_thread_start"), 0, entry + 1]);
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector(threaded ? [[1], [3]] : [[1]])),
    ...(threaded ? [] : section(5, [1, 0, ...uleb(pages)])),
    ...section(7, vector(exports)),
    ...section(10, vector([[...uleb(start.length), ...start], ...(threaded ? [body([])] : [])])),
  ]);
}

function run(options) {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-runner-"));
  try {
    const fixture = path.join(directory, "command.wasm");
    const bytes = options.bytes || command(options);
    assert.ok(WebAssembly.validate(bytes));
    fs.writeFileSync(fixture, bytes);
    const runner = options.threaded ? "run_losat_wasi_threads.js" : "run_losat_wasi.js";
    const result = spawnSync(process.execPath, [path.join(__dirname, runner), fixture], {
      env: hostEnv(), encoding: "utf8", timeout: 10000, killSignal: "SIGKILL",
    });
    assert.ifError(result.error);
    return result;
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
}

for (const threaded of [false, true]) {
  for (const status of [0, 23, 70]) {
    test(`${threaded ? "threaded" : "serial"} preserves WASI exit ${status}`, () => {
      const result = run({ threaded, status });
      assert.equal(result.signal, null, result.stderr);
      assert.equal(result.status, status, result.stderr);
    });
  }
}

test("worker memory trap terminates the command while its main thread waits", () => {
  const result = run({ threaded: true, worker: "trap" });
  assert.equal(result.signal, "SIGTERM", result.stderr);
  assert.match(result.stderr, /RuntimeError: memory access out of bounds/);
});

test("missing thread entry is rejected before starting the command", () => {
  const result = run({ threaded: true, worker: "missing-export" });
  assert.equal(result.signal, null, result.stderr);
  assert.equal(result.status, 1);
  assert.match(result.stderr, /requires wasi_thread_start/);
});

test("successful worker releases its parent and exits successfully", () => {
  const result = run({ threaded: true, worker: "success" });
  assert.equal(result.signal, null, result.stderr);
  assert.equal(result.status, 0, result.stderr);
});

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:68-93
// CBlastFormat::CBlastFormat(..., CNcbiOstream& outfile, ...) : ..., m_Outfile(outfile),
// The native command's output write blocks while a slow reader drains a full pipe. Starting a
// Worker makes Node open process.stdout, which sets a piped stdout non-blocking, so the
// guest's fd_write to a full pipe returned EAGAIN (WASI errno 6) and LOSAT stopped with
// "Resource temporarily unavailable (os error 6)". 64 writes of 4,096 bytes overfill the
// 64 KiB pipe of a reader that waits a second; every byte must arrive and the exit be 0.
test("threaded command output waits for a slow reader of a full pipe",
  { skip: process.platform === "win32" && "needs sh" }, () => {
    const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-runner-"));
    try {
      const fixture = path.join(directory, "command.wasm");
      const bytes = command({ threaded: true, worker: "success", writes: 64 });
      assert.ok(WebAssembly.validate(bytes));
      fs.writeFileSync(fixture, bytes);
      const result = spawnSync("sh", ["-c", '{ "$0" "$1" "$2"; echo "exit $?" >&2; } | { sleep 1; wc -c; }',
        process.execPath, path.join(__dirname, "run_losat_wasi_threads.js"), fixture], {
        env: hostEnv(), encoding: "utf8", timeout: 20000, killSignal: "SIGKILL",
      });
      assert.ifError(result.error);
      assert.match(result.stderr, /^exit 0$/m, result.stderr);
      assert.equal(Number(result.stdout.trim()), 64 * 4096, result.stderr);
    } finally {
      fs.rmSync(directory, { recursive: true, force: true });
    }
  });

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// node:wasi answers EINVAL (28) to a pointer argument that JavaScript sees as negative;
// both runners pass WASI arguments as unsigned values (wasi_artifact.js).
for (const threaded of [false, true]) {
  test(`${threaded ? "threaded" : "serial"} WASI call with a pointer above 2 GiB succeeds`, () => {
    const result = run({ threaded, bytes: highPointerCommand(threaded) });
    assert.equal(result.signal, null, result.stderr);
    assert.equal(result.status, 0, result.stderr);
  });
}

test("WASI imports receive i32 arguments as unsigned values", () => {
  const { unsignedWasiImports } = require("./wasi_artifact");
  const seen = [];
  const imports = unsignedWasiImports({
    clock_time_get: (...args) => { seen.push(args); return 0; },
    proc_exit: (...args) => { seen.push(args); },
  });
  assert.deepEqual(Object.keys(imports), ["clock_time_get", "proc_exit"]);
  assert.equal(imports.clock_time_get(1, 1n, -2147483640), 0);
  imports.proc_exit(-1);
  assert.deepEqual(seen, [[1, 1n, 2147483656], [4294967295]]);
});

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host memory only: when the declared maximum cannot be reserved, the host provides a
// shared memory with a smaller one. An engine either refuses to create the memory or
// creates one that cannot grow; an explicit maximum is never replaced.
test("shared memory steps down when the declared maximum cannot be reserved", () => {
  const { createSharedMemory } = require("./wasi_artifact");
  const tried = [];
  // `refuses`: the constructor throws above the limit. Otherwise the memory is created
  // with its initial pages only and cannot grow.
  const engine = (limit, refuses) => descriptor => {
    tried.push(descriptor.maximum);
    if (descriptor.maximum > limit && refuses) throw new RangeError("could not allocate memory");
    let pages = descriptor.initial;
    return { descriptor, pages: () => pages, grow(delta) {
      if (descriptor.maximum > limit) throw new RangeError("Unable to grow instance memory");
      pages += delta;
    } };
  };
  for (const refuses of [true, false]) {
    tried.length = 0;
    const whole = createSharedMemory(23, 65535, true, engine(65535, refuses));
    assert.deepEqual(whole.memory.descriptor, { initial: 22, maximum: 65535, shared: true });
    assert.deepEqual([whole.maximum, whole.memory.pages()], [65535, 23]);
    const less = createSharedMemory(23, 65535, true, engine(20000, refuses));
    assert.deepEqual([less.maximum, less.memory.pages()], [16384, 23]);
    assert.deepEqual(tried, [65535, 65535, 49152, 32768, 24576, 16384]);
    // Never below the 1 GiB of the earlier modules, nor below the initial size.
    assert.throws(() => createSharedMemory(23, 65535, true, engine(16383, refuses)), RangeError);
    assert.throws(() => createSharedMemory(23, 16384, true, engine(16383, refuses)), RangeError);
    tried.length = 0;
    assert.throws(() => createSharedMemory(30000, 65535, true, engine(20000, refuses)), RangeError);
    assert.deepEqual(tried, [65535, 49152, 32768]);
    // An explicit maximum is used as it is.
    tried.length = 0;
    assert.throws(() => createSharedMemory(23, 65535, false, engine(20000, refuses)), RangeError);
    assert.deepEqual(tried, [65535]);
  }
  assert.throws(() => createSharedMemory(23, 65535, true, () => { throw new TypeError("other"); }), TypeError);
  const real = createSharedMemory(3, 4, true);
  assert.equal(real.maximum, 4);
  assert.equal(real.memory.buffer.byteLength, 3 * 65536);
  assert.ok(real.memory.buffer instanceof SharedArrayBuffer);
});

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Check the Wasm-host preparation boundary independently from search behavior:
// metadata stays serializable and the returned module executes with the same status.
test("prepared serial module preserves metadata and executes the inspected bytes", async () => {
  const { prepareArtifact, inspectArtifact } = require("./wasi_artifact");
  const { WASI } = require("node:wasi");
  const bytes = command({ threaded: false, status: 23 });
  const prepared = prepareArtifact(bytes, "serial-command");
  assert.ok(prepared.module instanceof WebAssembly.Module);
  assert.deepEqual(prepared.identity, inspectArtifact(bytes, "serial-command"));
  assert.deepEqual(Object.keys(JSON.parse(JSON.stringify(prepared.identity))),
    ["kind", "sha256", "imports", "exports", "memory"]);
  const wasi = new WASI({ version: "preview1", returnOnExit: true });
  const instance = await WebAssembly.instantiate(prepared.module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  assert.equal(wasi.start(instance), 23);
  const secondWasi = new WASI({ version: "preview1", returnOnExit: true });
  const second = await WebAssembly.instantiate(prepared.module, {
    wasi_snapshot_preview1: secondWasi.wasiImport,
  });
  assert.notEqual(second.exports.memory, instance.exports.memory);
  new Uint8Array(instance.exports.memory.buffer)[0] = 99;
  assert.equal(new Uint8Array(second.exports.memory.buffer)[0], 0);
  assert.equal(secondWasi.start(second), 23);
});

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Wasm validation and ABI errors remain errors in both public entry points.
test("artifact preparation preserves validation and ABI rejection", () => {
  const { prepareArtifact, inspectArtifact } = require("./wasi_artifact");
  for (const inspect of [prepareArtifact, inspectArtifact]) {
    const api = ["losat_web_run_pair", "losat_web_alloc", "losat_web_dealloc",
      "losat_web_result_ptr", "losat_web_result_len", "losat_web_error_ptr", "losat_web_error_len"];
    for (const threaded of [false, true]) {
      const result = inspect(command({ threaded, entry: "_initialize", aliases: api }),
        `${threaded ? "threaded" : "serial"}-reactor`);
      assert.equal((result.identity || result).kind, `${threaded ? "threaded" : "serial"}-reactor`);
    }
    assert.throws(() => inspect(command({ threaded: false, entry: "_initialize", aliases: api.slice(1) })),
      /reactor is missing direct API exports/);
    assert.throws(() => inspect(command({ threaded: false, entry: null })),
      /exactly one of _start or _initialize/);
    assert.throws(() => inspect(command({ threaded: false, aliases: ["_initialize", ...api] })),
      /exactly one of _start or _initialize/);
    assert.throws(() => inspect(command({ threaded: false, exportMemory: false })),
      /must export its single memory/);
    assert.throws(() => inspect(command({ threaded: false, memoryLimits: [1, 1, 2, 1] })),
      WebAssembly.CompileError);
    assert.throws(() => inspect(command({ threaded: false, memoryLimits: [1, 3, 1, 2] })),
      /inconsistent serial artifact thread interface/);
    assert.throws(() => inspect(Buffer.from([0, 97, 115])), WebAssembly.CompileError);
    assert.throws(() => inspect(command({ threaded: false }), "serial-reactor"),
      /expected serial-reactor, got serial-command/);
    assert.throws(() => inspect(command({ threaded: true }), "serial-command"),
      /expected serial-command, got threaded-command/);
    assert.throws(() => inspect(command({ threaded: true, worker: "missing-export" })),
      /requires wasi_thread_start/);
  }
});

// Unset, the guard cache is in the user's cache directory; 0 turns it off; a default
// directory that cannot be used only means that the guard is computed again.
test("LOSAT_WASI_GUARD_CACHE is on by default in the user's cache directory, and 0 turns it off", () => {
  const { defaultGuardCacheDir } = require("./wasi_thread_host");
  const home = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-guard-home-"));
  try {
    const bytes = spawner(2, { bulk: true });
    const reference = runThreaded(bytes, [], { LOSAT_WASI_GUARD_CACHE: "0" });
    assert.equal(reference.status, 3, reference.stderr);
    const user = { HOME: home, USERPROFILE: home, XDG_CACHE_HOME: path.join(home, "xdg"), LOCALAPPDATA: path.join(home, "local") };
    const dir = defaultGuardCacheDir(user, process.platform, home);
    assert.ok(dir.startsWith(home), dir);
    for (const state of ["miss", "hit"]) {
      const result = runThreaded(bytes, [], { ...user, LOSAT_WASI_GUARD_CACHE: undefined, LOSAT_WASI_THREADS_DEBUG: "1" });
      assert.equal(result.status, 3, result.stderr);
      assert.equal(result.stdout, reference.stdout);
      assert.ok(result.stderr.includes(`guard cache ${state}: ${dir}${path.sep}guarded-`), result.stderr);
    }
    const off = runThreaded(bytes, [], { ...user, LOSAT_WASI_GUARD_CACHE: "0", LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(off.stdout, reference.stdout);
    assert.doesNotMatch(off.stderr, /guard cache/);
    // A default directory under a file cannot be created: the run computes the guard.
    const file = path.join(home, "file");
    fs.writeFileSync(file, "");
    const blocked = { ...user, XDG_CACHE_HOME: file, LOCALAPPDATA: file };
    if (process.platform === "darwin") {
      fs.rmSync(path.join(home, "Library"), { recursive: true, force: true });
      fs.writeFileSync(path.join(home, "Library"), "");
    }
    const fallback = runThreaded(bytes, [], { ...blocked, LOSAT_WASI_GUARD_CACHE: undefined, LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(fallback.status, 3, fallback.stderr);
    assert.equal(fallback.stdout, reference.stdout);
    assert.match(fallback.stderr, /guard cache unusable/);
  } finally {
    fs.rmSync(home, { recursive: true, force: true });
  }
});

test("the default guard cache directory follows each platform's convention", () => {
  const { defaultGuardCacheDir } = require("./wasi_thread_host");
  const home = os.homedir();
  assert.equal(defaultGuardCacheDir({ XDG_CACHE_HOME: "/x/cache" }, "linux"), path.join("/x/cache", "losat", "wasi-guard"));
  assert.equal(defaultGuardCacheDir({ XDG_CACHE_HOME: "relative" }, "linux"), path.join(home, ".cache", "losat", "wasi-guard"));
  assert.equal(defaultGuardCacheDir({}, "linux"), path.join(home, ".cache", "losat", "wasi-guard"));
  assert.equal(defaultGuardCacheDir({ XDG_CACHE_HOME: "/x/cache" }, "darwin"), path.join(home, "Library", "Caches", "losat", "wasi-guard"));
  assert.equal(defaultGuardCacheDir({ LOCALAPPDATA: "C:\\L" }, "win32"), path.join("C:\\L", "losat", "wasi-guard"));
  assert.equal(defaultGuardCacheDir({}, "win32"), path.join(home, "AppData", "Local", "losat", "wasi-guard"));
});

test("the guard cache keeps the entries used last", () => {
  const { guardSharedMemoryCached } = require("./wasi_thread_host");
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-guard-keep-"));
  const saved = process.env.LOSAT_WASI_GUARD_CACHE;
  try {
    process.env.LOSAT_WASI_GUARD_CACHE = directory;
    const modules = Array.from({ length: 10 }, (_, i) => spawner(i + 1, { bulk: true }));
    const sidecars = () => fs.readdirSync(directory).filter((file) => file.endsWith(".json"));
    const stamp = (bytes, seconds) => {
      const sha = require("node:crypto").createHash("sha256").update(bytes).digest("hex");
      const file = sidecars().find((name) => name.includes(sha));
      fs.utimesSync(path.join(directory, file), seconds, seconds);
    };
    modules.slice(0, 8).forEach((bytes, i) => {
      assert.equal(guardSharedMemoryCached(bytes).cache, "miss");
      stamp(bytes, 1000 + i);
    });
    // The oldest entry is used again (a hit), so the next two misses remove entries 1 and 2.
    assert.equal(guardSharedMemoryCached(modules[0]).cache, "hit");
    for (const bytes of modules.slice(8)) assert.equal(guardSharedMemoryCached(bytes).cache, "miss");
    assert.equal(sidecars().length, 8);
    assert.equal(fs.readdirSync(directory).filter((file) => file.endsWith(".wasm")).length, 8);
    assert.equal(guardSharedMemoryCached(modules[0]).cache, "hit");
    // Bytes without a sidecar and temporary files go once they are an hour old.
    const old = Date.now() / 1000 - 7200;
    const leftovers = [`guarded-${"a".repeat(64)}-${"b".repeat(64)}.wasm`, `guarded-${"c".repeat(64)}-${"d".repeat(64)}.json.1.ab.tmp`];
    for (const [i, name] of leftovers.entries()) {
      fs.writeFileSync(path.join(directory, name), "x");
      if (i === 0) fs.utimesSync(path.join(directory, name), old, old);
    }
    assert.equal(guardSharedMemoryCached(modules[1]).cache, "miss");
    const left = fs.readdirSync(directory);
    assert.equal(left.includes(leftovers[0]), false);
    assert.equal(left.includes(leftovers[1]), true);
  } finally {
    if (saved === undefined) delete process.env.LOSAT_WASI_GUARD_CACHE;
    else process.env.LOSAT_WASI_GUARD_CACHE = saved;
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:156-180
// ```c++
//     // Create the threads ...
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         thread->Reset(new CPrelimSearchThread(internal_data,
//     ...
//     // ... launch the threads ...
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     ...
//         (*thread)->Join(&result);
// ```
// A command that starts `spawns` threads (start_arg 1, 2, 4, ...; exit 90 when a spawn
// fails) and waits until every thread has added its start_arg and tid. It prints the
// thread count and 'a' + the tid sum and exits with the start_arg sum, so a lost,
// repeated or misdirected start changes stdout or the status. `bulk` adds a memory.fill
// and a memory.copy for the shared-memory guard. The host options under test change
// when workers boot, never this output.
// With `pauseAfter` k, the main thread spawns k threads, waits until they have counted
// themselves, waits 0.5 s more (memory.atomic.wait32 on a word that stays 0) and spawns the rest.
function spawner(spawns, { bulk = false, pauseAfter = 0 } = {}) {
  const types = vector([
    [0x60, 1, 0x7f, 0],                          // proc_exit(i32)
    [0x60, 0, 0],                                // _start()
    [0x60, 1, 0x7f, 1, 0x7f],                    // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],                    // wasi_thread_start(i32, i32)
    [0x60, 4, 0x7f, 0x7f, 0x7f, 0x7f, 1, 0x7f],  // fd_write(i32, i32, i32, i32) -> i32
  ]);
  const imports = [
    [...name("wasi_snapshot_preview1"), ...name("proc_exit"), 0, 0],
    [...name("env"), ...name("memory"), 2, 3, 1, ...uleb(16384)],
    [...name("wasi"), ...name("thread-spawn"), 0, 2],
    [...name("wasi_snapshot_preview1"), ...name("fd_write"), 0, 4],
  ];
  const load = (address) => [...i32(address), 0xfe, 0x10, 2, 0];                   // i32.atomic.load
  const add = (address, value) => [...i32(address), ...value, 0xfe, 0x1e, 2, 0, 0x1a]; // atomic add; drop
  const byte = (address, value) => [...i32(address), ...value, 0x3a, 0, 0];          // i32.store8
  const main = [];
  for (let i = 0; i < spawns; i++) {
    if (pauseAfter && i === pauseAfter) {
      // block { loop { if (count == pauseAfter) break; memory.atomic.wait32(0, count, -1) } }
      main.push(0x02, 0x40, 0x03, 0x40, ...load(0), ...i32(pauseAfter), 0x46, 0x0d, 1,
        ...i32(0), ...load(0), 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0c, 0, 0x0b, 0x0b);
      // memory.atomic.wait32(48, 0, 500 ms): times out
      main.push(...i32(48), ...i32(0), 0x42, ...i32(500_000_000).slice(1), 0xfe, 0x01, 2, 0, 0x1a);
    }
    // if (thread-spawn(1 << i) < 1) proc_exit(90)
    main.push(...i32(1 << i), 0x10, 1, ...i32(1), 0x48, 0x04, 0x40, ...i32(90), 0x10, 0, 0x0b);
  }
  main.push(
    // block { loop { if (count == spawns) break; memory.atomic.wait32(0, count, -1) } }
    0x02, 0x40, 0x03, 0x40, ...load(0), ...i32(spawns), 0x46, 0x0d, 1,
    ...i32(0), ...load(0), 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0c, 0, 0x0b, 0x0b,
    ...byte(32, [...load(0), ...i32(48), 0x6a]), ...byte(33, [...load(8), ...i32(97), 0x6a]),
    ...byte(34, i32(10)),
    ...i32(16), ...i32(32), 0x36, 2, 0, ...i32(20), ...i32(3), 0x36, 2, 0,  // iovec {32, 3}
    ...i32(1), ...i32(16), ...i32(1), ...i32(24), 0x10, 2, 0x1a,              // fd_write(1, ...)
    ...(bulk ? [...i32(64), ...i32(0), ...i32(0), 0xfc, 11, 0,                // memory.fill
      ...i32(64), ...i32(72), ...i32(0), 0xfc, 10, 0, 0] : []),               // memory.copy
    ...load(4), 0x10, 0,                                                       // proc_exit(sum)
  );
  const thread = [...add(4, [0x20, 1]), ...add(8, [0x20, 0]), ...add(0, i32(1)),
    ...i32(0), ...i32(1), 0xfe, 0x00, 2, 0, 0x1a];                             // notify(0, 1)
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector([[1], [3]])),
    ...section(7, vector([[...name("memory"), 2, 0], [...name("_start"), 0, 3],
      [...name("wasi_thread_start"), 0, 4]])),
    ...section(10, vector([body(main), body(thread)])),
  ]);
}

// The host's environment: no LOSAT_WASI_* of the caller, this file's guard cache, then
// `extra`, where `undefined` leaves a variable out.
function hostEnv(extra = {}) {
  const env = { ...process.env, NODE_NO_WARNINGS: "1" };
  for (const key of Object.keys(env)) {
    if (key.startsWith("LOSAT_WASI_") || key.startsWith("LOSAT_WASM_MEMORY_")) delete env[key];
  }
  const merged = { ...env, LOSAT_WASI_GUARD_CACHE: GUARD_CACHE, ...extra };
  return Object.fromEntries(Object.entries(merged).filter(([, value]) => value !== undefined));
}

// The command runner, or with `options` a driver that passes the test-only host options
// (fault, startTimeoutMs, rejectAt) through createThreadHost as check_wasi_fault.js does.
const DRIVER = `"use strict";
const { createThreadHost } = require(${JSON.stringify(path.join(__dirname, "wasi_thread_host.js"))});
(async () => {
  const [wasm, json, ...argv] = process.argv.slice(2);
  const { rejectAt, ...options } = JSON.parse(json);
  if (rejectAt) options.rejectSpawn = attempt => attempt === rejectAt;
  const host = await createThreadHost(wasm, "threaded-command", argv, options);
  try {
    process.exitCode = host.start();
    await host.waitForWorkers();
  } finally { await host.close(); }
})().catch(err => { console.error(err); process.exitCode = 1; });
`;
function runThreaded(bytes, argv, env = {}, options = null) {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-host-"));
  try {
    const fixture = path.join(directory, "command.wasm");
    assert.ok(WebAssembly.validate(bytes));
    fs.writeFileSync(fixture, bytes);
    let script = path.join(__dirname, "run_losat_wasi_threads.js");
    if (options) fs.writeFileSync(script = path.join(directory, "driver.js"), DRIVER);
    const result = spawnSync(process.execPath, [script, fixture, ...(options ? [JSON.stringify(options)] : []), ...argv], {
      env: hostEnv(env), encoding: "utf8", timeout: 20000, killSignal: "SIGKILL",
    });
    assert.ifError(result.error);
    return result;
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
}

const hostEvents = (stderr) => [...stderr.matchAll(/^\[losat-wasi-event\] (.+)$/gm)].map((m) => JSON.parse(m[1]));
const tidsOf = (events, type) => events.filter((e) => e.event === type).map((e) => e.tid).sort((a, b) => a - b);
const countOf = (events, type) => events.filter((e) => e.event === type).length;
const outcome = (result) => ({ status: result.status, signal: result.signal, stdout: result.stdout, stderr: result.stderr });

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:145-188
// TBlastThreads the_threads(GetNumberOfThreads()); (*thread)->Run(); (*thread)->Join(&result);
// The lifecycle events keep their meaning for the guest's threads; standby workers have
// their own events, and the unused ones are terminated at close.
function checkLifecycle(events, spawned, standby, label) {
  const tids = Array.from({ length: spawned }, (_, i) => i + 1);
  for (const type of ["spawn_attempt", "spawned", "ready", "exited"]) {
    assert.deepEqual(tidsOf(events, type), tids, `${label} ${type}`);
  }
  assert.ok(events.filter((e) => e.event === "exited").every((e) => e.code === 0), label);
  assert.equal(countOf(events, "spawn_rejected"), 0, label);
  assert.equal(countOf(events, "standby_spawned"), standby, label);
  assert.deepEqual(events.filter((e) => e.event === "spawned" && e.standby !== undefined).map((e) => e.tid),
    tids.slice(0, standby), label);
  assert.equal(countOf(events, "standby_unused"), Math.max(0, standby - spawned), label);
  assert.equal(countOf(events, "standby_exited"), Math.max(0, standby - spawned), label);
  if (standby && spawned) {
    assert.ok(events.findIndex((e) => e.event === "standby_spawned")
      < events.findIndex((e) => e.event === "spawn_attempt"), `${label}: standby boots before _start`);
  }
}

test("LOSAT_WASI_PRESPAWN keeps output and status with auto, fewer and more workers", () => {
  const bytes = spawner(3);
  const argv = ["-num_threads", "4"];
  const reference = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: "0" });
  assert.deepEqual(outcome(reference), { status: 7, signal: null, stdout: "3g\n", stderr: "" });
  for (const [prespawn, standby] of [["auto", 3], ["1", 1], ["5", 5], ["0", 0]]) {
    assert.deepEqual(outcome(runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: prespawn })), outcome(reference), prespawn);
    const debugged = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: prespawn, LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(debugged.status, 7, debugged.stderr);
    assert.equal(debugged.stdout, reference.stdout);
    checkLifecycle(hostEvents(debugged.stderr), 3, standby, `PRESPAWN=${prespawn}`);
  }
});

test("LOSAT_WASI_PRESPAWN=auto reads -num_threads, and unused standby workers are reaped", () => {
  const bytes = spawner(0);
  for (const [argv, standby] of [[["-num_threads", "3"], 2], [["-num_threads=5"], 4], [["--num_threads", "2"], 1],
    [[], 0], [["-num_threads", "1"], 0], [["-num_threads", "2147483647"], 0], [["-num_threads", "x"], 0]]) {
    const result = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: "auto", LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(result.signal, null, result.stderr);
    assert.equal(result.status, 0, result.stderr);
    assert.equal(result.stdout, "0a\n");
    checkLifecycle(hostEvents(result.stderr), 0, standby, argv.join(" "));
  }
});

// A spawn that finds no standby worker left boots the next set of -num_threads workers together
// with `auto`, and takes its first; a number of workers still boots one new worker per spawn after them.
test("LOSAT_WASI_PRESPAWN=auto boots the next set together when its standby workers run out", () => {
  const argv = ["-num_threads", "3"];
  for (const [spawns, prespawn, standby, refills] of [[5, "auto", 5, [3]], [4, "auto", 5, [3]], [6, "auto", 8, [3, 6]],
    [5, "2", 2, []]]) {
    const label = `${spawns} spawns PRESPAWN=${prespawn}`;
    const reference = runThreaded(spawner(spawns), argv, { LOSAT_WASI_PRESPAWN: "0" });
    // Without the reuse of returned workers (the threads of this guest return at once, so whether
    // the next spawn finds one is a race; the next test has the reuse).
    const result = runThreaded(spawner(spawns), argv,
      { LOSAT_WASI_PRESPAWN: prespawn, LOSAT_WASI_RECYCLE: "0", LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(result.status, reference.status, `${label}\n${result.stderr}`);
    assert.equal(result.stdout, reference.stdout, label);
    const events = hostEvents(result.stderr);
    checkLifecycle(events, spawns, standby, label);
    assert.deepEqual(events.filter((e) => e.event === "standby_refill").map((e) => [e.tid, e.count]),
      refills.map((tid) => [tid, 3]), label);
  }
});

// LOSAT_WASI_RECYCLE in a command (on by default): the standby workers whose threads returned run
// the next set of threads, before a refill boots new ones; `0` boots the refill as before.
test("a command's returned standby workers run its next set of threads (LOSAT_WASI_RECYCLE)", () => {
  const argv = ["-num_threads", "4"];
  const bytes = spawner(6, { pauseAfter: 3 });
  const reference = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: "0" });
  assert.equal(reference.status, 63, reference.stderr);
  for (const [recycle, reused, refills] of [[undefined, 3, 0], ["1", 3, 0], ["0", 0, 1]]) {
    const label = `LOSAT_WASI_RECYCLE=${recycle}`;
    const result = runThreaded(bytes, argv, { LOSAT_WASI_RECYCLE: recycle, LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.deepEqual(outcome({ ...result, stderr: "" }), outcome({ ...reference, stderr: "" }), label);
    const events = hostEvents(result.stderr);
    const tids = [1, 2, 3, 4, 5, 6];
    for (const type of ["spawn_attempt", "spawned", "ready", "exited"]) {
      assert.deepEqual(tidsOf(events, type), tids, `${label} ${type}`);
    }
    assert.ok(events.filter((e) => e.event === "exited").every((e) => e.code === 0), label);
    assert.deepEqual(events.filter((e) => e.event === "spawned" && e.standby !== undefined).map((e) => e.tid),
      recycle === "0" ? [1, 2, 3, 4, 5, 6] : [1, 2, 3], label);
    assert.deepEqual(events.filter((e) => e.event === "spawned" && e.recycled !== undefined).map((e) => e.tid),
      reused ? [4, 5, 6] : [], label);
    assert.equal(countOf(events, "standby_refill"), refills, label);
  }
});

test("invalid LOSAT_WASI_PRESPAWN values fail before the command starts", () => {
  for (const value of ["", "-1", "255", "two", "1.5"]) {
    const result = runThreaded(spawner(1), [], { LOSAT_WASI_PRESPAWN: value });
    assert.equal(result.status, 1, value);
    assert.equal(result.stdout, "", value);
    assert.match(result.stderr, /LOSAT_WASI_PRESPAWN must be auto or an integer from 0 to 254/, value);
  }
});

// NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
// LOG_POST(Error << "Error: " << e.what()); exit_code = BLAST_UNKNOWN_ERROR;
// Each test-only fault ends the command the same way with and without standby workers.
for (const [fault, message] of [
  ["startup-timeout", /timed out waiting for WASI thread 1 to start/],
  ["startup-error", /injected worker startup error/],
  ["trap", /injected worker trap/],
  ["abnormal-exit", /WASI worker 1 exited unexpectedly with code 23/],
]) {
  test(`fault ${fault} has the same result with LOSAT_WASI_PRESPAWN`, () => {
    for (const env of [{ LOSAT_WASI_PRESPAWN: "0" }, {}, { LOSAT_WASI_PRESPAWN: "auto" }, { LOSAT_WASI_PRESPAWN: "3" }]) {
      const result = runThreaded(spawner(1), ["-num_threads", "2"], env, { fault, startTimeoutMs: 1000 });
      const label = `${fault} ${JSON.stringify(env)}`;
      assert.equal(result.signal, "SIGTERM", `${label}\n${result.stderr}`);
      assert.equal(result.stdout, "", label);
      assert.match(result.stderr, message, label);
    }
  });
}

test("a rejected spawn has the same result with LOSAT_WASI_PRESPAWN", () => {
  const bytes = spawner(3);
  const argv = ["-num_threads", "4"];
  const reference = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: "0" }, { rejectAt: 2 });
  assert.deepEqual(outcome(reference), { status: 90, signal: null, stdout: "", stderr: "" });
  for (const [prespawn, standby] of [["auto", 3], ["1", 1]]) {
    assert.deepEqual(outcome(runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: prespawn }, { rejectAt: 2 })),
      outcome(reference), prespawn);
    const result = runThreaded(bytes, argv, { LOSAT_WASI_PRESPAWN: prespawn, LOSAT_WASI_THREADS_DEBUG: "1" },
      { rejectAt: 2 });
    assert.equal(result.status, 90, result.stderr);
    const events = hostEvents(result.stderr);
    assert.deepEqual(tidsOf(events, "spawn_attempt"), [1, 2], prespawn);
    assert.deepEqual(tidsOf(events, "spawn_rejected"), [2], prespawn);
    for (const type of ["spawned", "ready", "exited"]) assert.deepEqual(tidsOf(events, type), [1], `${prespawn} ${type}`);
    assert.equal(countOf(events, "standby_spawned"), standby, prespawn);
    assert.equal(countOf(events, "standby_unused"), standby - 1, prespawn);
    assert.equal(countOf(events, "standby_exited"), standby - 1, prespawn);
  }
});

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// The guard cache returns exactly the guard's bytes, cold or warm, and rebuilds an entry
// whose bytes or sidecar do not check out.
test("LOSAT_WASI_GUARD_CACHE returns the guard's bytes and rebuilds a bad entry", () => {
  const { guardSharedMemory } = require("./wasi_shared_memory");
  const { guardSharedMemoryCached } = require("./wasi_thread_host");
  const sha256 = (data) => require("node:crypto").createHash("sha256").update(data).digest("hex");
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-guard-cache-"));
  const saved = process.env.LOSAT_WASI_GUARD_CACHE;
  try {
    const cache = process.env.LOSAT_WASI_GUARD_CACHE = path.join(directory, "cache");
    for (const [bytes, guards] of [[spawner(1, { bulk: true }), [1, 1]], [spawner(1), [0, 0]]]) {
      const expected = guardSharedMemory(bytes);
      assert.deepEqual([expected.fill, expected.copy], guards);
      assert.equal(expected.bytes.equals(bytes), guards[0] + guards[1] === 0);
      const check = (state) => {
        const result = guardSharedMemoryCached(bytes);
        assert.equal(result.cache, state);
        assert.ok(result.bytes.equals(expected.bytes), state);
        assert.deepEqual([result.fill, result.copy], [expected.fill, expected.copy], state);
      };
      check("miss");
      check("hit");
      const files = fs.readdirSync(cache).filter((file) => file.includes(sha256(bytes))).sort();
      assert.equal(files.length, 2);
      const [sidecar, guarded] = files.map((file) => path.join(cache, file));
      assert.ok(fs.readFileSync(guarded).equals(expected.bytes));
      const record = JSON.parse(fs.readFileSync(sidecar, "utf8"));
      // A sidecar whose own SHA-256 matches its fields, so only the field checks reject it.
      const resealed = (changes) => {
        const { record_sha256: _, ...fields } = { ...record, ...changes };
        fs.writeFileSync(sidecar, JSON.stringify({ ...fields, record_sha256: sha256(JSON.stringify(fields)) }));
      };
      const corruptions = [
        () => { const data = fs.readFileSync(guarded); data[data.length - 1] ^= 1; fs.writeFileSync(guarded, data); },
        () => fs.truncateSync(guarded, 10),
        () => fs.writeFileSync(sidecar, "{"),
        () => fs.writeFileSync(sidecar, "null"),
        () => fs.writeFileSync(sidecar, JSON.stringify({ ...record, fill: record.fill + 1 })),
        () => fs.writeFileSync(sidecar, JSON.stringify({ ...record, copy: record.copy + 1 })),
        () => resealed({ fill: guards[0] + guards[1] === 0 ? 1 : 0, copy: 0 }),
        () => resealed({ fill: -1 }),
        () => resealed({ output_bytes: record.output_bytes + 1 }),
        () => resealed({ input_sha256: "0".repeat(64) }),
        () => resealed({ guard_version: "0".repeat(64) }),
        () => resealed({ format: 2 }),
      ];
      for (const corrupt of corruptions) {
        corrupt();
        check("invalid");
        check("hit");
      }
      fs.rmSync(guarded);
      check("miss");
      check("hit");
    }
    assert.deepEqual(fs.readdirSync(cache).filter((file) => file.endsWith(".tmp")), []);
  } finally {
    if (saved === undefined) delete process.env.LOSAT_WASI_GUARD_CACHE;
    else process.env.LOSAT_WASI_GUARD_CACHE = saved;
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

test("LOSAT_WASI_GUARD_CACHE runs give the default host's output, cold and warm", () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-guard-run-"));
  try {
    const bytes = spawner(2, { bulk: true });
    const reference = runThreaded(bytes, []);
    assert.deepEqual(outcome(reference), { status: 3, signal: null, stdout: "2d\n", stderr: "" });
    const env = (cache) => ({ LOSAT_WASI_GUARD_CACHE: path.join(directory, cache), LOSAT_WASI_PRESPAWN: "2" });
    assert.deepEqual(outcome(runThreaded(bytes, [], env("a"))), outcome(reference), "cold");
    assert.deepEqual(outcome(runThreaded(bytes, [], env("a"))), outcome(reference), "warm");
    for (const [cache, state] of [["b", "miss"], ["b", "hit"]]) {
      const result = runThreaded(bytes, [], { ...env(cache), LOSAT_WASI_THREADS_DEBUG: "1" });
      assert.equal(result.status, 3, result.stderr);
      assert.equal(result.stdout, reference.stdout);
      assert.match(result.stderr, new RegExp(`guard cache ${state}: `));
      assert.match(result.stderr, /shared-memory guards: 1 fills, 1 copies/);
    }
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

// NCBI reference (598d8ae6): c++/src/algo/blast/api/prelim_stage.cpp:156-180
// ```c++
//     // Create the threads ...
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         thread->Reset(new CPrelimSearchThread(internal_data,
//     ...
//     // ... launch the threads ...
//     NON_CONST_ITERATE(TBlastThreads, thread, the_threads) {
//         (*thread)->Run();
//     ...
//         (*thread)->Join(&result);
// ```
// A threaded reactor whose `run(spawns)` export clears its counters, starts `spawns` threads
// (start_arg 1, 2, 4, ...; -1 when a spawn fails), waits until every thread has added its
// start_arg and tid, releases them, and returns the start_arg sum (the tid sum stays at
// address 8). The threads wait for the release, so the threads of one call run at the same
// time. Each call is one search of `spawns` + 1 threads; the API exports that the artifact
// check asks for are aliases of `_initialize` and are never called.
function reactorSpawner() {
  const types = vector([
    [0x60, 0, 0],                                // _initialize()
    [0x60, 1, 0x7f, 1, 0x7f],                    // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],                    // wasi_thread_start(i32, i32)
  ]);
  const imports = [
    [...name("env"), ...name("memory"), 2, 3, 1, ...uleb(16384)],
    [...name("wasi"), ...name("thread-spawn"), 0, 1],
  ];
  const load = (address) => [...i32(address), 0xfe, 0x10, 2, 0];                   // i32.atomic.load
  const store = (address, value) => [...i32(address), ...value, 0xfe, 0x17, 2, 0];  // i32.atomic.store
  const add = (address, value) => [...i32(address), ...value, 0xfe, 0x1e, 2, 0, 0x1a]; // atomic add; drop
  const spawns = [0x20, 0], i = 1;                                             // local 0 = spawns, local 1 = i
  const run = [...store(0, i32(0)), ...store(4, i32(0)), ...store(8, i32(0)), ...store(12, i32(0)),
    // block { loop { if (i >= spawns) break; if (thread-spawn(1 << i) < 1) return -1; i += 1 } }
    0x02, 0x40, 0x03, 0x40, 0x20, i, ...spawns, 0x4e, 0x0d, 1,
    ...i32(1), 0x20, i, 0x74, 0x10, 0, ...i32(1), 0x48, 0x04, 0x40, ...i32(-1), 0x0f, 0x0b,
    0x20, i, ...i32(1), 0x6a, 0x21, i, 0x0c, 0, 0x0b, 0x0b,
    // block { loop { if (count == spawns) break; memory.atomic.wait32(0, count, -1) } }
    0x02, 0x40, 0x03, 0x40, ...load(0), ...spawns, 0x46, 0x0d, 1,
    ...i32(0), ...load(0), 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0c, 0, 0x0b, 0x0b,
    ...store(12, i32(1)), ...i32(12), ...spawns, 0xfe, 0x00, 2, 0, 0x1a,          // release; notify(12, spawns)
    ...load(4),
  ];
  const withLocal = (instructions) => {                                          // one i32 local
    const bytes = [1, 1, 0x7f, ...instructions, 0x0b];
    return [...uleb(bytes.length), ...bytes];
  };
  const thread = [...add(4, [0x20, 1]), ...add(8, [0x20, 0]), ...add(0, i32(1)),
    ...i32(0), ...i32(1), 0xfe, 0x00, 2, 0, 0x1a,                              // notify(0, 1)
    // block { loop { if (released) break; memory.atomic.wait32(12, 0, -1) } }
    0x02, 0x40, 0x03, 0x40, ...load(12), ...i32(1), 0x46, 0x0d, 1,
    ...i32(12), ...i32(0), 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0c, 0, 0x0b, 0x0b];
  const api = ["losat_web_run_pair", "losat_web_alloc", "losat_web_dealloc",
    "losat_web_result_ptr", "losat_web_result_len", "losat_web_error_ptr", "losat_web_error_len"];
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector([[0], [1], [2]])),
    ...section(7, vector([[...name("memory"), 2, 0], [...name("_initialize"), 0, 1],
      ...api.map((alias) => [...name(alias), 0, 1]), [...name("run"), 0, 2],
      [...name("wasi_thread_start"), 0, 3]])),
    ...section(10, vector([body([]), withLocal(run), body(thread)])),
  ]);
}

// A reactor's calls get the same threads, in the same order, with recycled workers
// (LOSAT_WASI_RECYCLE, on by default) and without them. Recycled workers run the threads of
// later calls; a call that needs more threads than there are idle workers boots new ones
// instead of waiting for a worker whose thread runs (the threads of a call wait for each other).
test("a reactor's thread workers outlive their threads (LOSAT_WASI_RECYCLE)", async () => {
  const { createThreadHost } = require("./wasi_thread_host");
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-recycle-"));
  const saved = Object.fromEntries(["LOSAT_WASI_RECYCLE", "LOSAT_WASI_PRESPAWN"].map((key) => [key, process.env[key]]));
  try {
    const fixture = path.join(directory, "reactor.wasm");
    const bytes = reactorSpawner();
    assert.ok(WebAssembly.validate(bytes));
    fs.writeFileSync(fixture, bytes);
    delete process.env.LOSAT_WASI_PRESPAWN;
    for (const recycle of [undefined, "1", "0"]) {
      if (recycle === undefined) delete process.env.LOSAT_WASI_RECYCLE;
      else process.env.LOSAT_WASI_RECYCLE = recycle;
      const label = `LOSAT_WASI_RECYCLE=${recycle}`;
      const host = await createThreadHost(fixture, "threaded-reactor");
      try {
        // Threads per call, and the recycled workers each call takes: none in the first call, the
        // three of the first call (plus two new ones) for five threads, then two of the five.
        let next = 1;
        for (const [call, spawns, reusedCount] of [[0, 3, 0], [1, 5, 3], [2, 2, 2], [3, 5, 5]]) {
          assert.equal(host.instance.exports.run(spawns), 2 ** spawns - 1, `${label} call ${call}`);
          const events = await host.waitForWorkers();
          const tids = Array.from({ length: spawns }, (_, k) => next + k);
          next += spawns;
          assert.equal(new Int32Array(host.memory.buffer, 8, 1)[0], tids.reduce((a, b) => a + b), `${label} call ${call}`);
          for (const type of ["spawn_attempt", "spawned", "ready", "exited"]) {
            assert.deepEqual(tidsOf(events, type), tids, `${label} call ${call} ${type}`);
          }
          assert.ok(events.filter((e) => e.event === "exited").every((e) => e.code === 0), label);
          const reused = events.filter((e) => e.event === "spawned" && e.recycled !== undefined).length;
          assert.equal(reused, recycle === "0" ? 0 : reusedCount, `${label} call ${call} recycled`);
        }
      } finally { await host.close(); }
    }
    process.env.LOSAT_WASI_RECYCLE = "yes";
    await assert.rejects(createThreadHost(fixture, "threaded-reactor"), /LOSAT_WASI_RECYCLE must be 0 or 1/);
  } finally {
    for (const [key, value] of Object.entries(saved)) {
      if (value === undefined) delete process.env[key];
      else process.env[key] = value;
    }
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

test("LOSAT_WASI_PRESPAWN is auto when unset and 0 boots no worker", () => {
  const bytes = spawner(3);
  const argv = ["-num_threads", "4"];
  for (const [env, standby] of [[{}, 3], [{ LOSAT_WASI_PRESPAWN: "0" }, 0]]) {
    const result = runThreaded(bytes, argv, { ...env, LOSAT_WASI_THREADS_DEBUG: "1" });
    assert.equal(result.status, 7, result.stderr);
    assert.equal(result.stdout, "3g\n");
    checkLifecycle(hostEvents(result.stderr), 3, standby, JSON.stringify(env));
  }
});

// NCBI reference (598d8ae6): c++/src/objtools/blast/seqdb_reader/seqdbatlas.cpp:114-127
// ```c++
// CMemoryFile* CSeqDBAtlas::GetMemoryFile(const string& fileName)
// ...
//     m_FileMemMap[fileName].reset(file);
// ```
// A command whose main thread opens `file` (path_open on the preopened "/") and whose one
// thread reads up to 64 bytes of it with that descriptor (fd_read); the main thread waits,
// writes what the thread read to standard output and exits with the thread's errno.
function descriptorReader(file) {
  const types = vector([
    [0x60, 1, 0x7f, 0],                          // proc_exit(i32)
    [0x60, 0, 0],                                // _start()
    [0x60, 1, 0x7f, 1, 0x7f],                    // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],                    // wasi_thread_start(i32, i32)
    [0x60, 4, 0x7f, 0x7f, 0x7f, 0x7f, 1, 0x7f],  // fd_write / fd_read(i32, i32, i32, i32) -> i32
    // path_open(fd, dirflags, path, path_len, oflags, rights_base, rights_inheriting, fdflags, fd_out) -> i32
    [0x60, 9, 0x7f, 0x7f, 0x7f, 0x7f, 0x7f, 0x7e, 0x7e, 0x7f, 0x7f, 1, 0x7f],
  ]);
  const imports = [
    [...name("wasi_snapshot_preview1"), ...name("proc_exit"), 0, 0],
    [...name("env"), ...name("memory"), 2, 3, 1, ...uleb(16384)],
    [...name("wasi"), ...name("thread-spawn"), 0, 2],
    [...name("wasi_snapshot_preview1"), ...name("fd_write"), 0, 4],
    [...name("wasi_snapshot_preview1"), ...name("path_open"), 0, 5],
    [...name("wasi_snapshot_preview1"), ...name("fd_read"), 0, 4],
  ];
  const load = (address) => [...i32(address), 0x28, 2, 0];                          // i32.load
  const store = (address, value) => [...i32(address), ...value, 0x36, 2, 0];         // i32.store
  const pathBytes = Buffer.from(file);
  // 0: thread done, 4: its errno, 8: the descriptor, 16: read iovec, 24: bytes read,
  // 32: write iovec, 40: bytes written, 64: the path, 256: the bytes.
  const main = [...pathBytes].flatMap((value, i) => [...i32(64 + i), ...i32(value), 0x3a, 0, 0]);
  main.push(
    // if (path_open(3, 0, 64, len, 0, FD_READ, 0, 0, 8) != 0) proc_exit(100)
    ...i32(3), ...i32(0), ...i32(64), ...i32(pathBytes.length), ...i32(0), 0x42, 2, 0x42, 0, ...i32(0), ...i32(8),
    0x10, 3, 0x04, 0x40, ...i32(100), 0x10, 0, 0x0b,
    // if (thread-spawn(0) < 1) proc_exit(90)
    ...i32(0), 0x10, 1, ...i32(1), 0x48, 0x04, 0x40, ...i32(90), 0x10, 0, 0x0b,
    // block { loop { if (done) break; memory.atomic.wait32(0, 0, -1) } }
    0x02, 0x40, 0x03, 0x40, ...i32(0), 0xfe, 0x10, 2, 0, 0x0d, 1,
    ...i32(0), ...i32(0), 0x42, 0x7f, 0xfe, 0x01, 2, 0, 0x1a, 0x0c, 0, 0x0b, 0x0b,
    ...store(32, i32(256)), ...store(36, load(24)),
    ...i32(1), ...i32(32), ...i32(1), ...i32(40), 0x10, 2, 0x1a,               // fd_write(1, ...)
    ...load(4), 0x10, 0,                                                         // proc_exit(errno)
  );
  const thread = [...store(16, i32(256)), ...store(20, i32(64)),
    ...i32(4), ...load(8), ...i32(16), ...i32(1), ...i32(24), 0x10, 4, 0x36, 2, 0, // errno = fd_read(fd, ...)
    ...i32(0), ...i32(1), 0xfe, 0x17, 2, 0,                                       // done (atomic store)
    ...i32(0), ...i32(1), 0xfe, 0x00, 2, 0, 0x1a];                               // notify(0, 1)
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector([[1], [3]])),
    ...section(7, vector([[...name("memory"), 2, 0], [...name("_start"), 0, 5],
      [...name("wasi_thread_start"), 0, 6]])),
    ...section(10, vector([body(main), body(thread)])),
  ]);
}

test("a descriptor that the main thread opened reads on a thread worker (-num_threads 2)", () => {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-descriptor-"));
  try {
    const file = path.join(directory, "data.txt");
    fs.writeFileSync(file, "read on another thread\n");
    const bytes = descriptorReader(path.relative("/", file));
    for (const env of [{}, { LOSAT_WASI_PRESPAWN: "0" }]) {
      const result = runThreaded(bytes, ["-num_threads", "2"], env);
      assert.deepEqual(outcome(result), { status: 0, signal: null, stdout: "read on another thread\n", stderr: "" },
        JSON.stringify(env));
    }
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
});
