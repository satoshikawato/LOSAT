"use strict";

const assert = require("node:assert/strict");
const { execFileSync, spawnSync } = require("node:child_process");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { Worker } = require("node:worker_threads");
const { test } = require("node:test");
const { guardSharedMemory, guardRecord, verifyGuarded, GUARD_SECTION } = require("./wasi_shared_memory");

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Check the host memory boundary independently of the parallel BLAST search.
// Fixtures use the WebAssembly Core binary encoding and original bulk-memory
// operations as the oracle, including traps; no alignment behavior is simulated.
const uleb = (value) => {
  const out = [];
  do { const b = value & 127; value >>>= 7; out.push(b | (value ? 128 : 0)); } while (value);
  return out;
};
const vector = (items) => [...uleb(items.length), ...items.flat()];
const name = (s) => [...uleb(Buffer.byteLength(s)), ...Buffer.from(s)];
const section = (id, bytes) => [id, ...uleb(bytes.length), ...bytes];
const body = (ops) => { const b = [0, ...ops, 0x0b]; return [...uleb(b.length), ...b]; };
const args = [0x20, 0, 0x20, 1, 0x20, 2];

function fixture({ shared = true, extra = [], custom = true, maximum = 4 } = {}) {
  // Imported function, table and global exercise index counting. All original
  // calls, exports and element-segment function indices must remain valid.
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...(custom ? section(0, [...name("fixture"), 0xfc, 11, 0, 0xfc, 10, 0, 0]) : []),
    ...section(1, [1, 0x60, 3, 0x7f, 0x7f, 0x7f, 0]),
    ...section(2, vector([
      [...name("env"), ...name("noop"), 0, 0],
      [...name("env"), ...name("memory"), 2, shared ? 3 : 1, 1, ...uleb(maximum)],
      [...name("env"), ...name("table"), 1, 0x70, 1, 1, 1],
      [...name("env"), ...name("global"), 3, 0x7f, 0],
    ])),
    ...section(3, vector([[0], [0], [0], ...extra.map(() => [0])])),
    ...section(7, vector([
      [...name("fill"), 0, 1], [...name("copy"), 0, 2], [...name("indirect"), 0, 3],
      ...extra.map((_, i) => [...name(`extra${i}`), 0, 4 + i]),
    ])),
    ...section(9, [1, 0, 0x41, 0, 0x0b, 1, 1]), // table[0] = original fill
    ...section(10, vector([
      body([...args, 0x10, 0, ...args, 0xfc, 11, 0]),
      body([...args, 0xfc, 10, 0, 0]),
      body([...args, 0x41, 0, 0x11, 0, 0]),
      ...extra.map(body),
    ])),
    ...section(11, [1, 0, 0x41, 0, 0x0b, 7, 0xfc, 11, 0, 0xfc, 10, 0, 0]),
  ]);
}

function instance(bytes, memory = new WebAssembly.Memory({ initial: 1, maximum: 4, shared: true })) {
  let calls = 0;
  const module = new WebAssembly.Module(bytes);
  const exports = new WebAssembly.Instance(module, { env: {
    memory, noop() { calls += 1; }, global: 0,
    table: new WebAssembly.Table({ element: "anyfunc", initial: 1, maximum: 1 }),
  } }).exports;
  return { exports, memory, calls: () => calls };
}

test("bulk guards preserve import, export, call and element indices", () => {
  const bytes = fixture();
  const original = Buffer.from(bytes);
  const guarded = guardSharedMemory(bytes);
  assert.deepEqual(bytes, original, "input artifact is immutable");
  assert.equal(guarded.fill, 1);
  assert.equal(guarded.copy, 1);
  const i = instance(guarded.bytes);
  i.exports.indirect(100, 257, 7);
  assert.equal(i.calls(), 1);
  assert.deepEqual([...new Uint8Array(i.memory.buffer, 100, 7)], Array(7).fill(1));
  assert.deepEqual([...new Uint8Array(i.memory.buffer, 0, 7)], [0xfc, 11, 0, 0xfc, 10, 0, 0]);
});

for (const operation of ["fill", "copy"]) {
  test(`${operation}: unsigned ranges, overlap, zero length and traps match original Wasm`, () => {
    const original = instance(fixture());
    const guarded = instance(guardSharedMemory(fixture()).bytes);
    const cases = operation === "fill"
      ? [[10, -1, 50], [65530, 258, 6], [65536, 1, 0], [65537, 1, 0],
        [65535, 1, 2], [-1, 1, 0], [0, 1, -1], [-8, 1, 16]]
      : [[11, 10, 30], [10, 11, 30], [65530, 20, 6], [20, 65530, 6],
        [65536, 65536, 0], [65537, 0, 0], [0, 65537, 0], [65535, 0, 2],
        [0, 65535, 2], [-1, 0, 0], [0, -1, 0], [0, 0, -1], [-8, 0, 16]];
    for (const parameters of cases) {
      for (const i of [original, guarded]) {
        const view = new Uint8Array(i.memory.buffer);
        for (let p = 0; p < view.length; p += 1) view[p] = p & 255;
      }
      let trapped = false;
      try { original.exports[operation](...parameters); }
      catch (err) { assert.ok(err instanceof WebAssembly.RuntimeError); trapped = true; }
      if (trapped) assert.throws(() => guarded.exports[operation](...parameters), WebAssembly.RuntimeError);
      else guarded.exports[operation](...parameters);
      assert.deepEqual(new Uint8Array(guarded.memory.buffer), new Uint8Array(original.memory.buffer),
        `${operation}(${parameters})`);
      assert.equal(guarded.memory.buffer.byteLength, 65536, "grow(0) allocates no pages");
    }
  });
}

test("opcode-looking constants, SIMD immediates, locals and memargs are not rewritten", () => {
  const ops = [
    0x41, 0xfc, 0x0b, 0x1a, // signed LEB constant containing a fill prefix
    0x44, 0xfc, 11, 0, 0xfc, 10, 0, 0, 0, 0x1a, // f64 constant
    0xfd, 12, 0xfc, 11, 0, 0xfc, 10, 0, 0, ...Array(9).fill(0),
    0xfd, 0x15, 0, 0x1a, // i8x16.extract_lane_s
    0x41, 0, 0x28, 2, 0xfc, 0x0b, 0x1a, // load offset LEB containing fill prefix
    0x41, 0, 0xfe, 0x10, 2, 0, 0x1a, // i32.atomic.load
    0xfe, 3, 0, // atomic.fence
    0x02, 0x40, 0x41, 0, 0x0e, 1, 0, 0, 0x0b, // br_table
  ];
  const bytes = fixture({ extra: [ops] });
  assert.ok(WebAssembly.validate(bytes));
  const guarded = guardSharedMemory(bytes);
  assert.equal(guarded.fill, 1);
  assert.equal(guarded.copy, 1);
  instance(guarded.bytes).exports.extra0(0, 0, 0);
});

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host memory boundary only. LOSAT's threaded modules declare 65,535 pages
// (LOSAT/.cargo/config.toml). The guard's ranges are unsigned 64-bit, so fills and
// copies across 2 GiB and up to the last byte behave as in the original module.
test("fill and copy across 2 GiB and at the top of a 65,535-page memory match original Wasm", (t) => {
  const pages = 65535, end = pages * 65536, half = 2 ** 31;
  const windows = [[half - 64, 128], [end - 64, 64]];
  const grown = (bytes) => {
    const memory = new WebAssembly.Memory({ initial: 1, maximum: pages, shared: true });
    memory.grow(pages - 1); // address space only: untouched pages are not resident
    return instance(bytes, memory);
  };
  let original, guarded;
  try {
    original = grown(fixture({ maximum: pages }));
    guarded = grown(guardSharedMemory(fixture({ maximum: pages })).bytes);
  } catch (err) {
    if (!(err instanceof RangeError)) throw err;
    t.skip(`this machine cannot reserve two 4 GiB memories: ${err.message}`);
    return;
  }
  assert.equal(guarded.memory.buffer.byteLength, end);
  const cases = [
    ["fill", [end - 16, 0x5a, 16]], ["fill", [half - 8, 0x33, 16]], ["fill", [end, 1, 0]],
    ["fill", [end - 16, 1, 17]], ["fill", [end + 1, 1, 0]], ["fill", [2 ** 32 - 16, 1, 32]],
    ["copy", [end - 16, half - 8, 16]], ["copy", [half - 8, end - 16, 16]], ["copy", [half, half - 4, 32]],
    ["copy", [end, end, 0]], ["copy", [end - 8, half, 16]], ["copy", [half, end - 8, 16]],
  ];
  for (const [operation, parameters] of cases) {
    for (const i of [original, guarded]) {
      for (const [at, length] of windows) {
        const view = new Uint8Array(i.memory.buffer, at, length);
        for (let p = 0; p < length; p += 1) view[p] = (at + p) & 255;
      }
    }
    let trapped = false;
    try { original.exports[operation](...parameters); }
    catch (err) { assert.ok(err instanceof WebAssembly.RuntimeError); trapped = true; }
    if (trapped) assert.throws(() => guarded.exports[operation](...parameters), WebAssembly.RuntimeError);
    else guarded.exports[operation](...parameters);
    for (const [at, length] of windows) {
      assert.deepEqual(new Uint8Array(guarded.memory.buffer, at, length),
        new Uint8Array(original.memory.buffer, at, length), `${operation}(${parameters}) at ${at}`);
    }
  }
});

test("invalid modules and unsupported memory layout fail explicitly", () => {
  assert.throws(() => guardSharedMemory(fixture().subarray(0, 30)), /invalid input Wasm/);
  assert.throws(() => guardSharedMemory(fixture({ shared: false })), /shared wasm32 memory/);
});

for (const [operation, parameters] of [
  ["fill", [65536, 83, 8]],
  ["copy destination", [65536, 83, 8]],
  ["copy source", [64, 65536, 8]],
]) {
  test(`another isolate grows memory before ${operation} in a running Wasm consumer`, { timeout: 15000 }, async () => {
    // Signal the grower, then remain inside Wasm while waiting for publication.
    // Copy from and fill within the new page without returning through JS first.
    const consumer = [
      0x41, 0, 0x41, 1, 0xfe, 0x17, 2, 0,
      0x03, 0x40,
      0x41, 0, 0xfe, 0x10, 2, 0, 0x41, 2, 0x47, 0x0d, 0,
      0x0b,
      ...args, 0xfc, ...(operation === "fill" ? [11, 0] : [10, 0, 0]),
    ];
    const memory = new WebAssembly.Memory({ initial: 1, maximum: 4, shared: true });
    const i = instance(guardSharedMemory(fixture({ extra: [consumer] })).bytes, memory);
    new Uint8Array(memory.buffer, 83, 8).fill(83);
    const worker = new Worker(`
      const { workerData, parentPort } = require("node:worker_threads");
      const view = new Int32Array(workerData.buffer);
      while (Atomics.load(view, 0) !== 1) Atomics.wait(view, 0, 0, 10);
      workerData.grow(1);
      new Uint8Array(workerData.buffer, 65536, 8).fill(83);
      Atomics.store(view, 0, 2);
      parentPort.postMessage("grown");
    `, { eval: true, workerData: memory });
    // A separate worker deadline protects the synchronous Wasm call too.
    const watchdog = new Worker(`
      const { workerData, parentPort } = require("node:worker_threads");
      const timer = setTimeout(() => { Atomics.store(new Int32Array(workerData), 0, 2); }, 8000);
      parentPort.on("message", () => { clearTimeout(timer); parentPort.close(); });
    `, { eval: true, workerData: memory.buffer });
    const done = new Promise((resolve, reject) => {
      worker.once("message", resolve);
      worker.once("error", reject);
      worker.once("exit", (code) => { if (code) reject(new Error(`grower exit ${code}`)); });
    });
    try {
      i.exports.extra0(...parameters);
      assert.equal(await done, "grown");
      assert.equal(memory.buffer.byteLength, 131072);
      assert.deepEqual([...new Uint8Array(memory.buffer, parameters[0], 8)], Array(8).fill(83));
    } finally {
      watchdog.postMessage("done");
      await Promise.all([worker.terminate(), watchdog.terminate()]);
    }
  });
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Build step (`wasi_shared_memory.js guard`, release v0.2.0 DECISIONS 52): the guarded module
// carries its record, is never guarded twice, verifies, and a module that cannot be guarded
// fails the build without an output.
const marker = (record) => Buffer.from(section(0, [...name(GUARD_SECTION), ...Buffer.from(JSON.stringify(record))]));
// A module with shared memory and no memory.fill or memory.copy.
const plain = () => Buffer.from([
  0, 97, 115, 109, 1, 0, 0, 0,
  ...section(1, [1, 0x60, 0, 0]),
  ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, 4]])),
  ...section(3, [1, 0]),
  ...section(10, vector([body([])])),
]);
// The code section of `bytes` with function `index` (of the defined ones) replaced by `replacement`.
function replaceBody(bytes, index, replacement) {
  const leb = (at) => { let value = 0, shift = 0, b; do { b = bytes[at.pos++]; value |= (b & 127) << shift; shift += 7; } while (b & 128); return value; };
  const at = { pos: 8 };
  const parts = [bytes.subarray(0, 8)];
  while (at.pos < bytes.length) {
    const start = at.pos, id = bytes[at.pos++], size = leb(at), end = at.pos + size;
    if (id !== 10) { parts.push(bytes.subarray(start, end)); at.pos = end; continue; }
    const count = leb(at), bodies = [];
    for (let i = 0; i < count; i += 1) { const length = leb(at); bodies.push(bytes.subarray(at.pos, at.pos + length)); at.pos += length; }
    bodies[index] = Buffer.from(replacement);
    const code = [...uleb(count), ...bodies.flatMap((b) => [...uleb(b.length), ...b])];
    parts.push(Buffer.from(section(10, code)));
  }
  return Buffer.concat(parts);
}

test("the guarded module carries its record and is returned as it is when guarded again", () => {
  const guarded = guardSharedMemory(fixture());
  assert.deepEqual(guarded.record, { version: 3, fill: 1, copy: 1, grow: 0, atomic: 0, atomic_helpers: 0 });
  assert.deepEqual(guardRecord(guarded.bytes), guarded.record);
  assert.deepEqual(guarded.bytes.subarray(-marker(guarded.record).length), marker(guarded.record));
  assert.equal(guardRecord(fixture()), null);
  const again = guardSharedMemory(guarded.bytes);
  assert.equal(again.bytes, guarded.bytes);
  assert.deepEqual([again.fill, again.copy, again.record], [0, 0, guarded.record]);
  assert.deepEqual(verifyGuarded(guarded.bytes), guarded.record);
  // Nothing to guard: the input comes back without a record.
  const unchanged = guardSharedMemory(plain());
  assert.deepEqual([unchanged.bytes, unchanged.record], [plain(), null]);
});

test("verifyGuarded rejects unguarded, forged, partly guarded and miscounted modules", () => {
  assert.throws(() => verifyGuarded(fixture()), /no losat\.shared_memory_guard section/);
  const functions = (...types) => Buffer.from([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...section(1, vector([[0x60, 0, 0]])),
    ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, 4]])),
    ...section(3, vector(types.map((type) => [type]))),
    ...section(10, vector(types.map(() => body([])))),
    ...marker({ version: 3, fill: 1, copy: 0, grow: 0, atomic: 0, atomic_helpers: 0 }),
  ]);
  assert.ok(WebAssembly.validate(functions(0)));
  assert.throws(() => verifyGuarded(functions(0)), /the guard's helpers are missing/);
  assert.throws(() => verifyGuarded(functions(0, 0)), /the guard's helpers are missing/);
  assert.throws(() => verifyGuarded(functions(0, 0, 0)), /do not have the guard's type/);
  assert.throws(() => verifyGuarded(Buffer.concat([fixture(), marker({ version: 3, fill: 1, copy: 1, grow: 0, atomic: 0, atomic_helpers: 0 })])),
    /not the guard's bodies/);
  const original = fixture({ extra: [[...args, 0xfc, 10, 0, 0]] });
  const guarded = guardSharedMemory(original);
  assert.deepEqual(guarded.record, { version: 3, fill: 1, copy: 2, grow: 0, atomic: 0, atomic_helpers: 0 });
  // extra0 (defined function 3) with its original memory.copy instead of the helper call.
  const partly = replaceBody(guarded.bytes, 3, body([...args, 0xfc, 10, 0, 0]).slice(1));
  assert.ok(WebAssembly.validate(partly));
  assert.throws(() => verifyGuarded(partly), /function 4 has an unguarded memory\.copy/);
  const tail = marker(guarded.record).length;
  const recount = (record) => Buffer.concat([guarded.bytes.subarray(0, -tail), marker(record)]);
  assert.throws(() => verifyGuarded(recount({ version: 3, fill: 2, copy: 2, grow: 0, atomic: 0, atomic_helpers: 0 })), /record counts 2 fills, 2 copies, 0 grows and 0 atomics/);
  assert.throws(() => guardSharedMemory(recount({ version: 2, fill: 1, copy: 2, atomic: 0, atomic_helpers: 0 })), /unsupported losat\.shared_memory_guard record/);
  assert.throws(() => guardSharedMemory(recount({ version: 3, fill: 0, copy: 0, grow: 0, atomic: 0, atomic_helpers: 0 })), /unsupported losat\.shared_memory_guard record/);
});

test("the guard command writes a verified module; a decline or nothing to guard exits 1 without output", () => {
  const command = require.resolve("./wasi_shared_memory");
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "losat-guard-command-"));
  const file = (name) => path.join(dir, name);
  const json = (...argv) => JSON.parse(execFileSync(process.execPath, [command, ...argv], { encoding: "utf8" }));
  const fails = (argv, pattern) => {
    const result = spawnSync(process.execPath, [command, ...argv], { encoding: "utf8" });
    assert.equal(result.status, 1, result.stderr);
    assert.match(result.stderr, pattern);
    assert.equal(result.stdout, "");
  };
  try {
    fs.writeFileSync(file("in.wasm"), fixture());
    const written = json("guard", file("in.wasm"), file("out.wasm"));
    assert.deepEqual([written.version, written.fill, written.copy, written.guarded_before], [3, 1, 1, false]);
    assert.deepEqual(fs.readFileSync(file("out.wasm")), guardSharedMemory(fixture()).bytes);
    assert.equal(json("verify", file("out.wasm")).sha256, written.output_sha256);
    // A guarded input is written as it is.
    const again = json("guard", file("out.wasm"), file("again.wasm"));
    assert.deepEqual([again.guarded_before, again.output_sha256, again.fill], [true, written.output_sha256, 1]);
    // An instruction that the reader declines (GC ref.i31), and a module with nothing to guard.
    fs.writeFileSync(file("gc.wasm"), fixture({ extra: [[0x41, 0, 0xfb, 0x1c, 0x1a]] }));
    assert.ok(WebAssembly.validate(fs.readFileSync(file("gc.wasm"))));
    fails(["guard", file("gc.wasm"), file("gc-out.wasm")], /unsupported Wasm opcode 0xfb/);
    fs.writeFileSync(file("plain.wasm"), plain());
    fails(["guard", file("plain.wasm"), file("plain-out.wasm")], /no memory\.fill, memory\.copy, memory\.grow or atomic instruction to guard/);
    fails(["verify", file("in.wasm")], /no losat\.shared_memory_guard section/);
    fails(["guard", file("in.wasm")], /usage/);
    assert.deepEqual(fs.readdirSync(dir).sort(), ["again.wasm", "gc.wasm", "in.wasm", "out.wasm", "plain.wasm"]);
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
});

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Atomic memory instructions: the guarded helper is the original instruction, and refreshes a
// stale bound first (the header of wasi_shared_memory.js).
const I32 = 0x7f, I64 = 0x7e;
const WIDTHS = [[I32, 4], [I64, 8], [I32, 1], [I32, 2], [I64, 1], [I64, 2], [I64, 4]];
const ATOMICS = [
  { sub: 0, values: [I32], result: [I32], size: 4, name: "memory.atomic.notify" },
  { sub: 1, values: [I32, I64], result: [I32], size: 4, name: "memory.atomic.wait32" },
  { sub: 2, values: [I64, I64], result: [I32], size: 8, name: "memory.atomic.wait64" },
  ...WIDTHS.flatMap(([type, size], k) => [
    { sub: 0x10 + k, values: [], result: [type], size, name: `load ${k}` },
    { sub: 0x17 + k, values: [type], result: [], size, name: `store ${k}` },
    ...[0x1e, 0x25, 0x2c, 0x33, 0x3a, 0x41].map((base) => ({ sub: base + k, values: [type], result: [type], size, name: `rmw 0x${base.toString(16)} ${k}` })),
    { sub: 0x48 + k, values: [type, type], result: [type], size, name: `cmpxchg ${k}` },
  ]),
];
const alignOf = (size) => Math.log2(size);

// One exported function per entry of `ops`: its parameters are the instruction's operands.
function atomicFixture(ops) {
  const types = ops.map(({ values, result }) => [0x60, values.length + 1, I32, ...values, result.length, ...result]);
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, vector(types)),
    ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, 4]])),
    ...section(3, vector(ops.map((_, i) => uleb(i)))),
    ...section(7, vector(ops.map((_, i) => [...name(`op${i}`), 0, ...uleb(i)]))),
    ...section(10, vector(ops.map(({ sub, values, size, offset }) => body([
      0x20, 0, ...values.flatMap((_, i) => [0x20, i + 1]),
      0xfe, ...uleb(sub), alignOf(size), ...uleb(offset),
    ])))),
  ]);
}

function runAtomic(fn, memory, args) {
  try { return { result: fn(...args), bytes: Buffer.from(new Uint8Array(memory.buffer)).toString("hex") }; }
  catch (err) { return { trap: err.message, bytes: Buffer.from(new Uint8Array(memory.buffer)).toString("hex") }; }
}

test("every atomic instruction: guarded results, memory and traps match original Wasm", () => {
  for (const offset of [0, 8, 65528, 2 ** 32 - 8]) {
    const ops = ATOMICS.map((op) => ({ ...op, offset }));
    const original = new WebAssembly.Module(atomicFixture(ops));
    const guardedResult = guardSharedMemory(atomicFixture(ops));
    assert.equal(guardedResult.atomic, ops.length);
    assert.equal(guardedResult.record.atomic_helpers, ops.length);
    assert.ok(verifyGuarded(guardedResult.bytes));
    const guarded = new WebAssembly.Module(guardedResult.bytes);
    ops.forEach((op, i) => {
      const value = (type, v) => (type === I64 ? BigInt(v) : v);
      // A timed wait returns at once: not-equal when the operand differs, timed-out at 0 ns.
      const operands = op.sub === 1 || op.sub === 2
        ? [[value(op.values[0], 7), 0n], [value(op.values[0], 0), 0n]]
        : [op.values.map((type, k) => value(type, 3 + k))];
      const top = 65536 - op.size - offset;
      // In bounds, the last address in bounds, one page and one step past it, misaligned, wrapped.
      const addresses = [offset >= 65528 ? 0 : 16, top, top + op.size, 2 ** 31 - op.size, -op.size,
        ...(op.size > 1 ? [1] : [])];
      for (const args of operands) {
        for (const address of addresses) {
          const make = (module) => {
            const memory = new WebAssembly.Memory({ initial: 1, maximum: 4, shared: true });
            new Uint8Array(memory.buffer).fill(0);
            const view = new DataView(memory.buffer);
            for (let at = 0; at < 65536; at += 4096) view.setUint32(at, at + 1, true);
            return { memory, fn: new WebAssembly.Instance(module, { env: { memory } }).exports[`op${i}`] };
          };
          const a = make(original);
          const b = make(guarded);
          assert.deepEqual(runAtomic(b.fn, b.memory, [address, ...args]), runAtomic(a.fn, a.memory, [address, ...args]),
            `${op.name} offset ${offset} address ${address}`);
        }
      }
    });
  }
});

test("helpers are shared by opcode and memarg, and the guard adds no import or export", () => {
  const op = { sub: 0x1e, values: [I32], result: [I32], size: 4 };
  const ops = [{ ...op, offset: 0 }, { ...op, offset: 0 }, { ...op, offset: 8 }, { ...ATOMICS[3], offset: 0 }];
  const { record, bytes } = guardSharedMemory(atomicFixture(ops));
  assert.deepEqual([record.atomic, record.atomic_helpers], [4, 3]);
  const before = new WebAssembly.Module(atomicFixture(ops));
  const after = new WebAssembly.Module(bytes);
  assert.deepEqual(WebAssembly.Module.imports(after), WebAssembly.Module.imports(before));
  assert.deepEqual(WebAssembly.Module.exports(after), WebAssembly.Module.exports(before));
});

test("atomic.fence is not rewritten", () => {
  const fence = [0xfe, 3, 0];
  const original = fixture({ extra: [[...fence]] });
  const guarded = guardSharedMemory(original);
  assert.equal(guarded.atomic, 0);
  assert.ok(verifyGuarded(guarded.bytes));
});

test("verifyGuarded rejects an unguarded atomic, a forged helper and a wrong atomic count", () => {
  const add = { sub: 0x1e, values: [I32], result: [I32], size: 4, offset: 8 };
  const ops = [add, { ...add, offset: 0 }];
  const guarded = guardSharedMemory(atomicFixture(ops));
  assert.deepEqual(verifyGuarded(guarded.bytes), guarded.record);
  const original = [0x20, 0, 0x20, 1, 0xfe, 0x1e, 2, 8];
  const partly = replaceBody(guarded.bytes, 0, body(original).slice(1));
  assert.ok(WebAssembly.validate(partly));
  assert.throws(() => verifyGuarded(partly), /function 0 has an unguarded atomic instruction 0xfe 30/);
  // The first atomic helper without its bound check is the same operation, and is not the guard's.
  const forged = replaceBody(guarded.bytes, 5, body(original).slice(1));
  assert.ok(WebAssembly.validate(forged));
  assert.throws(() => verifyGuarded(forged), /not the guard's bodies/);
  const tail = marker(guarded.record).length;
  const recount = (record) => Buffer.concat([guarded.bytes.subarray(0, -tail), marker(record)]);
  assert.throws(() => verifyGuarded(recount({ ...guarded.record, atomic: 3 })), /0 fills, 0 copies, 0 grows and 3 atomics/);
  assert.throws(() => guardSharedMemory(recount({ ...guarded.record, atomic_helpers: 3 })), /unsupported losat\.shared_memory_guard record/);
});

// The stale bound itself, on the engine, for every atomic instruction. On V8 a helper's own entry
// stack check usually handles the pending growth first, and a tier-up check in a loop may too, so
// this runs in a child process that stays in Liftoff with stack checks off (--liftoff-only
// --no-wasm-stack-checks): only the helper's compares and memory.grow(0) can refresh the bound.
// Each probe tells the main thread that it entered, waits in a plain-load loop until the main
// thread has grown the memory through the module and published an address in the new page, and
// returns memory.size, the same bound. That the original Wasm is stale here (or traps) is
// reported, not asserted: an engine may refresh the bound itself. A failure ends in the child's
// timeout, never in a hung test process.
const STALE_CHILD = `
  const { Worker } = require("node:worker_threads");
  const module = new WebAssembly.Module(require("node:fs").readFileSync(process.argv[1]));
  const memory = new WebAssembly.Memory({ initial: 1, maximum: 256, shared: true });
  const view = new Int32Array(memory.buffer);
  const { grow, ...probes } = new WebAssembly.Instance(module, { env: { memory } }).exports;
  const worker = new Worker(\`
    const { workerData, parentPort } = require("node:worker_threads");
    const { module, memory } = workerData;
    const view = new Int32Array(memory.buffer);
    const { grow, ...probes } = new WebAssembly.Instance(module, { env: { memory } }).exports;
    const outcomes = {};
    for (const [label, probe] of Object.entries(probes)) {
      try { outcomes[label] = probe(); } catch (err) { outcomes[label] = err.message; }
      Atomics.store(view, 1, 0);
    }
    parentPort.postMessage(outcomes);
  \`, { eval: true, workerData: { module, memory } });
  worker.once("message", (outcomes) => { console.log(JSON.stringify(outcomes)); process.exit(0); });
  worker.once("error", (err) => { console.log("worker error: " + err.message); process.exit(1); });
  for (const _ of Object.keys(probes)) {
    while (Atomics.load(view, 0) === 0) Atomics.wait(view, 0, 0, 10);
    Atomics.store(view, 0, 0);
    const old = grow(1);
    Atomics.store(view, 1, old * 65536 + 64);
  }
`;

function runStaleChild(bytes, directory, label) {
  const file = path.join(directory, `${label}.wasm`);
  fs.writeFileSync(file, bytes);
  const result = spawnSync(process.execPath, ["--liftoff-only", "--no-wasm-stack-checks", "-e", STALE_CHILD, file],
    { encoding: "utf8", timeout: 30000 });
  assert.equal(result.status, 0, `${label}: status ${result.status} signal ${result.signal}: ${result.stderr}`);
  return JSON.parse(result.stdout);
}

// The atomic instructions that read memory: all but the stores and notify.
const readsMemory = ({ sub }) => sub !== 0 && !(sub >= 0x17 && sub <= 0x1d);

// probes: { label: [ops of the function] }, called in this order, one page grown for each.
function staleProbes(probes) {
  const entered = 0, pointer = 4;
  const labels = Object.keys(probes);
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, vector([[0x60, 0, 1, I32], [0x60, 1, I32, 1, I32]])),
    ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, ...uleb(256)]])),
    ...section(3, vector([...labels.map(() => [0]), [1]])),
    ...section(7, vector([...labels.map((label, i) => [...name(label), 0, i]), [...name("grow"), 0, labels.length]])),
    ...section(10, vector([...labels.map((label) => body([
      0x41, entered, 0x41, 1, 0xfe, 0x17, 2, 0,
      0x41, entered, 0x41, 1, 0xfe, 0, 2, 0, 0x1a,
      0x03, 0x40, 0x41, pointer, 0x28, 2, 0, 0x45, 0x0d, 0, 0x0b,
      ...probes[label],
      0x3f, 0,
    ])), body([0x20, 0, 0x40, 0])])),
  ]);
}

const published = [0x41, 4, 0x28, 2, 0];
// An instruction whose expected values differ from the memory: wait returns at once.
const atomicOn = (address, { sub, values, result, size }, offset = 0) => [
  ...address, ...values.flatMap((type) => (type === I64 ? [0x42, 1] : [0x41, 1])),
  0xfe, ...uleb(sub), alignOf(size), ...uleb(offset), ...(result.length ? [0x1a] : []),
];

test("another isolate grows memory through the module: every atomic instruction refreshes a stale bound", { timeout: 90000 }, (t) => {
  const probes = {};
  // Before the instruction: the address itself is in the new page, in the operand or in the memarg.
  probes.offset = atomicOn([0x41, 0], ATOMICS[6 + 3], 65536 + 64);
  ATOMICS.forEach((op) => { probes[`before ${op.name}`] = atomicOn(published, op); });
  // After it: the instruction reads memory elsewhere, then a plain store uses the published address.
  // V8 checks that store against the stale bound on ARM64 outside macOS.
  const store = [...published, 0x41, 7, 0x36, 2, 0];
  ATOMICS.filter(readsMemory).forEach((op) => { probes[`after ${op.name}`] = [...atomicOn([0x41, 32], op), ...store]; });
  const original = staleProbes(probes);
  const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-guard-stale-"));
  try {
    const outcomes = runStaleChild(guardSharedMemory(original).bytes, directory, "guarded");
    const expected = Object.fromEntries(Object.keys(probes).map((label, i) => [label, i + 2]));
    assert.deepEqual(outcomes, expected);
    const stale = Object.entries(runStaleChild(original, directory, "original")).filter(([label, value]) => value !== expected[label]);
    t.diagnostic(`original Wasm: ${stale.length} of ${Object.keys(probes).length} probes stale or trapped`);
  } finally {
    fs.rmSync(directory, { recursive: true, force: true });
  }
});

// One exported memory.grow, the record's only count.
const growFixture = () => Buffer.from([
  0, 97, 115, 109, 1, 0, 0, 0,
  ...section(1, vector([[0x60, 1, I32, 1, I32]])),
  ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, 4]])),
  ...section(3, vector([[0]])),
  ...section(7, vector([[...name("grow"), 0, 0]])),
  ...section(10, vector([body([0x20, 0, 0x40, 0])])),
]);

test("memory.grow: guarded results and size match original Wasm, and the largest grown size is recorded", () => {
  const CELL = 16 / 4;
  const run = (bytes, start = 0) => {
    const memory = new WebAssembly.Memory({ initial: 1, maximum: 4, shared: true });
    new Int32Array(memory.buffer)[CELL] = start;
    const { grow } = new WebAssembly.Instance(new WebAssembly.Module(bytes), { env: { memory } }).exports;
    return [99, 1, 0, 2, 1, 0x7fffffff].map((pages) => {
      const result = grow(pages);
      return [result, memory.buffer.byteLength, new Int32Array(memory.buffer)[CELL]];
    });
  };
  const guarded = guardSharedMemory(growFixture());
  assert.deepEqual(guarded.record, { version: 3, fill: 0, copy: 0, grow: 1, atomic: 0, atomic_helpers: 0 });
  assert.deepEqual(verifyGuarded(guarded.bytes), guarded.record);
  const page = 65536;
  // A failed grow writes nothing; a grow by 0 records the size too; the original leaves the cell alone.
  assert.deepEqual(run(guarded.bytes), [[-1, page, 0], [1, 2 * page, 2], [2, 2 * page, 2], [2, 4 * page, 4], [-1, 4 * page, 4], [-1, 4 * page, 4]]);
  assert.deepEqual(run(growFixture()).map(([result, bytes]) => [result, bytes]), run(guarded.bytes).map(([result, bytes]) => [result, bytes]));
  // A larger size that another thread recorded stays.
  assert.deepEqual(run(guarded.bytes, 100).map((row) => row[2]), [100, 100, 100, 100, 100, 100]);
});

test("verifyGuarded rejects an unguarded memory.grow, a forged helper and a wrong grow count", () => {
  const guarded = guardSharedMemory(growFixture());
  const original = body([0x20, 0, 0x40, 0]).slice(1);
  const partly = replaceBody(guarded.bytes, 0, original);
  assert.ok(WebAssembly.validate(partly));
  assert.throws(() => verifyGuarded(partly), /function 0 has an unguarded memory\.grow/);
  // The helper without the record of the size is the same operation, and is not the guard's.
  const forged = replaceBody(guarded.bytes, 3, original);
  assert.ok(WebAssembly.validate(forged));
  assert.throws(() => verifyGuarded(forged), /not the guard's bodies/);
  const tail = marker(guarded.record).length;
  const recount = (record) => Buffer.concat([guarded.bytes.subarray(0, -tail), marker(record)]);
  assert.throws(() => verifyGuarded(recount({ ...guarded.record, grow: 2 })), /0 fills, 0 copies, 2 grows and 0 atomics/);
  assert.throws(() => guardSharedMemory(recount({ version: 3, fill: 0, copy: 0, atomic: 0, atomic_helpers: 0 })), /unsupported losat\.shared_memory_guard record/);
});
