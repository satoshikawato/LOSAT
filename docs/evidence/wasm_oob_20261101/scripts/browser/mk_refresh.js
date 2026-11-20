"use strict";
// Writes refresh_original.wasm and refresh_guarded.wasm: a plain store to an address in memory that
// another thread grew, after an i32.atomic.load of another address; probe returns memory.size.
// The same module as the tests in LOSAT/tests/test_wasi_shared_memory.js.
// usage: node mk_refresh.js <wasi_shared_memory.js>
const fs = require("fs");
const uleb = (v) => { const o = []; do { const b = v & 127; v >>>= 7; o.push(b | (v ? 128 : 0)); } while (v); return o; };
const vec = (items) => [...uleb(items.length), ...items.flat()];
const sec = (id, b) => [id, ...uleb(b.length), ...b];
const str = (s) => [...uleb(s.length), ...Buffer.from(s)];
const body = (ops) => { const b = [0, ...ops, 0x0b]; return [...uleb(b.length), ...b]; };
const fence = Array.from({ length: 60000 }, () => [0xfe, 3, 0]).flat();
const bytes = Buffer.from([0, 97, 115, 109, 1, 0, 0, 0,
  ...sec(1, vec([[0x60, 1, 0x7f, 1, 0x7f]])),
  ...sec(2, vec([[...str("env"), ...str("memory"), 2, 3, 1, 4]])),
  ...sec(3, vec([[0], [0]])),
  ...sec(7, vec([[...str("probe"), 0, 0], [...str("grow"), 0, 1]])),
  ...sec(10, vec([body([
    0x41, 0, 0x41, 1, 0xfe, 0x17, 2, 0,
    0x41, 0, 0x41, 1, 0xfe, 0, 2, 0, 0x1a,
    ...fence,
    0x03, 0x40, 0x41, 4, 0x28, 2, 0, 0x45, 0x0d, 0, 0x0b,
    0x41, 8, 0xfe, 0x10, 2, 0, 0x1a,
    0x41, 4, 0x28, 2, 0, 0x41, 7, 0x36, 2, 0,
    0x3f, 0,
  ]), body([0x20, 0, 0x40, 0])])),
]);
fs.writeFileSync(`${__dirname}/refresh_original.wasm`, bytes);
const guarded = require(process.argv[2]).guardSharedMemory(bytes);
fs.writeFileSync(`${__dirname}/refresh_guarded.wasm`, guarded.bytes);
console.log(guarded.record);
