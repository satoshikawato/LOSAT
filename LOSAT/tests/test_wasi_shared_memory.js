"use strict";

const assert = require("node:assert/strict");
const { Worker } = require("node:worker_threads");
const { test } = require("node:test");
const { guardSharedMemory } = require("./wasi_shared_memory");

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

function fixture({ shared = true, extra = [], custom = true } = {}) {
  // Imported function, table and global exercise index counting. All original
  // calls, exports and element-segment function indices must remain valid.
  return Buffer.from([
    0, 97, 115, 109, 1, 0, 0, 0,
    ...(custom ? section(0, [...name("fixture"), 0xfc, 11, 0, 0xfc, 10, 0, 0]) : []),
    ...section(1, [1, 0x60, 3, 0x7f, 0x7f, 0x7f, 0]),
    ...section(2, vector([
      [...name("env"), ...name("noop"), 0, 0],
      [...name("env"), ...name("memory"), 2, shared ? 3 : 1, 1, 4],
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
