"use strict";
// Line WO, V8 alone (no LOSAT code). Does an atomic in a Liftoff function trap on a page that
// another thread grew, while a plain access to the same address passes?
//
// The worker runs `probe(kind)`, one function with no call and no loop after its entry, so its
// prologue is its only stack check: it tells the main thread it has entered, runs a long straight
// line of atomic.fence, reads the address the main thread published, then touches that address.
// The main thread grows the memory while the fences run, then publishes the first new page.
// Prediction if the cause is the stale per-instance bound: plain = ok, atomic = trap,
// refresh (memory.grow(0) first) = ok. Each kind runs ROUNDS times. On ARM64 outside macOS,
// V8 also bounds checks plain multi-byte stores (kPartialOOBWritesAreNoops): i32.store = trap,
// i32.load and i32.store8 = ok.
//
// usage: node stale_atomic.js [rounds] [fences]
const { Worker, isMainThread, workerData, parentPort } = require("worker_threads");

const uleb = v => { const o = []; do { let b = v & 127; v >>>= 7; if (v) b |= 128; o.push(b); } while (v); return o; };
const sleb = v => { const o = []; for (;;) { const b = v & 127; v >>= 7; if ((v === 0 && !(b & 64)) || (v === -1 && (b & 64))) { o.push(b); return o; } o.push(b | 128); } };
const vec = items => [...uleb(items.length), ...items.flat()];
const sec = (id, body) => [id, ...uleb(body.length), ...body];
const str = s => [...uleb(s.length), ...Buffer.from(s)];

const ENTERED = 0; // i32 cell: worker -> main
const POINTER = 4; // i32 cell: main -> worker

function moduleBytes(fences) {
  const fence = [];
  for (let i = 0; i < fences; i += 1) fence.push(0xfe, 0x03, 0x00);
  const c = v => [0x41, ...sleb(v)];
  // (func $probe (param $kind i32) (result i32) (local $p i32) ...)
  const body = [
    1, 1, 0x7f, // one i32 local ($p = local 1)
    ...c(ENTERED), ...c(1), 0xfe, 0x17, 2, 0, // i32.atomic.store ENTERED = 1
    ...c(ENTERED), ...c(1), 0xfe, 0x00, 2, 0, 0x1a, // memory.atomic.notify, drop
    ...fence,
    ...c(POINTER), 0x28, 2, 0, 0x21, 1, // $p = i32.load POINTER (plain)
    // kind 2: refresh this instance's view first
    0x20, 0, ...c(2), 0x46, 0x04, 0x40, ...c(0), 0x40, 0, 0x1a, 0x0b,
    // kind 0: plain i32.store, 3: plain i32.load, 4: plain i32.store8, each returning;
    // kinds 1, 2: i32.atomic.rmw.add
    0x20, 0, 0x45, 0x04, 0x40, 0x20, 1, ...c(7), 0x36, 2, 0, ...c(0), 0x0f, 0x0b,
    0x20, 0, ...c(3), 0x46, 0x04, 0x40, 0x20, 1, 0x28, 2, 0, 0x0f, 0x0b,
    0x20, 0, ...c(4), 0x46, 0x04, 0x40, 0x20, 1, ...c(7), 0x3a, 0, 0, ...c(0), 0x0f, 0x0b,
    0x20, 1, ...c(7), 0xfe, 0x1e, 2, 0, 0x0b,
  ];
  return Buffer.from([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...sec(1, vec([[0x60, 1, 0x7f, 1, 0x7f]])),
    ...sec(2, vec([[...str("env"), ...str("memory"), 2, 3, 1, ...uleb(1024)]])),
    ...sec(3, vec([[0]])),
    ...sec(7, vec([[...str("probe"), 0, 0]])),
    ...sec(10, vec([[...uleb(body.length), ...body]])),
  ]);
}

if (!isMainThread) {
  const { module, memory } = workerData;
  const instance = new WebAssembly.Instance(module, { env: { memory } });
  parentPort.on("message", ({ kind }) => {
    let outcome;
    try { instance.exports.probe(kind); outcome = "ok"; }
    catch (err) { outcome = err.message; }
    parentPort.postMessage({ outcome });
  });
  parentPort.postMessage({ ready: true });
} else {
  (async () => {
    const rounds = Number(process.argv[2] || 50);
    const fences = Number(process.argv[3] || 100000);
    const module = new WebAssembly.Module(moduleBytes(fences));
    const names = ["plain i32.store", "atomic rmw.add", "atomic rmw.add after memory.grow(0)", "plain i32.load", "plain i32.store8"];
    const counts = names.map(() => new Map());
    const order = [];
    for (let r = 0; r < rounds; r += 1) order.push(0, 1, 2, 3, 4);
    for (const kind of order) {
      // A new worker and a new memory each time: every instance is fresh, so `probe` is in Liftoff.
      const memory = new WebAssembly.Memory({ initial: 1, maximum: 1024, shared: true });
      const view = new Int32Array(memory.buffer);
      const worker = new Worker(__filename, { workerData: { module, memory } });
      await new Promise(resolve => worker.once("message", resolve));
      const done = new Promise(resolve => worker.once("message", resolve));
      worker.postMessage({ kind });
      while (Atomics.load(view, ENTERED / 4) === 0) Atomics.wait(view, ENTERED / 4, 0, 100);
      const old = memory.grow(1);
      Atomics.store(view, POINTER / 4, old * 65536 + 64);
      const { outcome } = await done;
      await worker.terminate();
      counts[kind].set(outcome, (counts[kind].get(outcome) || 0) + 1);
    }
    console.log(`node ${process.version} v8 ${process.versions.v8} ${process.platform} ${process.arch}, ${rounds} rounds per kind, ${fences} fences`);
    names.forEach((n, k) => console.log(`${n}: ${[...counts[k]].map(([o, c]) => `${o} x${c}`).join(", ")}`));
  })();
}
