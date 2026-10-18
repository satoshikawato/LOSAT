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

function run(options) {
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-wasi-runner-"));
  try {
    const fixture = path.join(directory, "command.wasm");
    const bytes = command(options);
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
function spawner(spawns, { bulk = false } = {}) {
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

function hostEnv(extra = {}) {
  const env = { ...process.env, NODE_NO_WARNINGS: "1" };
  for (const key of Object.keys(env)) {
    if (key.startsWith("LOSAT_WASI_") || key.startsWith("LOSAT_WASM_MEMORY_")) delete env[key];
  }
  return { ...env, ...extra };
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
