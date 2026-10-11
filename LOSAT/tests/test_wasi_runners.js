"use strict";

const assert = require("node:assert/strict");
const { spawnSync } = require("node:child_process");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { test } = require("node:test");

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
const body = (instructions) => {
  const bytes = [0, ...instructions, 0x0b];
  return [...uleb(bytes.length), ...bytes];
};

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Entry/export/memory variations exercise Wasm ABI inspection only.
function command({ threaded, status = 0, worker = null, entry = "_start",
  aliases = [], exportMemory = true, memoryLimits = [1, 0, 1] }) {
  const types = vector([
    [0x60, 1, 0x7f, 0],                 // proc_exit(i32)
    [0x60, 0, 0],                       // _start()
    [0x60, 1, 0x7f, 1, 0x7f],          // thread-spawn(i32) -> i32
    [0x60, 2, 0x7f, 0x7f, 0],          // wasi_thread_start(i32, i32)
  ]);
  const imports = [[...name("wasi_snapshot_preview1"), ...name("proc_exit"), 0, 0]];
  if (threaded) {
    imports.push([...name("env"), ...name("memory"), 2, 3, 21, ...uleb(16384)]);
    imports.push([...name("wasi"), ...name("thread-spawn"), 0, 2]);
  }
  const exports = exportMemory ? [[...name("memory"), 2, 0]] : [];
  for (const exported of [entry, ...aliases].filter(Boolean)) {
    exports.push([...name(exported), 0, threaded ? 2 : 1]);
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
    if (worker !== "missing-export") exports.push([...name("wasi_thread_start"), 0, 3]);
    bodies.push(body(worker === "trap"
      ? [...i32(-1), ...i32(0), ...i32(1), 0xfc, 11, 0] // OOB memory.fill
      : [...i32(0), ...i32(1), 0xfe, 0x17, 2, 0,      // atomic store + notify
        ...i32(0), ...i32(1), 0xfe, 0, 2, 0, 0x1a]));
  }
  if (threaded && !worker) {
    functions.push(3);
    exports.push([...name("wasi_thread_start"), 0, 3]);
    bodies.push(body([]));
  }
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, types),
    ...section(2, vector(imports)),
    ...section(3, vector(functions.map((type) => [type]))),
    ...(threaded ? [] : section(5, memoryLimits)),
    ...section(7, vector(exports)),
    ...section(10, vector([body(main), ...bodies])),
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
    const env = { ...process.env, NODE_NO_WARNINGS: "1" };
    for (const key of Object.keys(env)) {
      if (key.startsWith("LOSAT_WASI_") || key.startsWith("LOSAT_WASM_MEMORY_")) delete env[key];
    }
    const result = spawnSync(process.execPath, [path.join(__dirname, runner), fixture], {
      env, encoding: "utf8", timeout: 10000, killSignal: "SIGKILL",
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
