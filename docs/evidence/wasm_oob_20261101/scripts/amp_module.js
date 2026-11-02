"use strict";
// Line WO: widens one window of the real module, for the trap-rate experiments only.
// Finds the crossbeam `Collector` clone in a thread-local's get_or_init_slow:
//   i32.const 0; i32.atomic.load A (the Once); i32.eqz; br_if 0; call X; end;
//   i32.const 0; i32.load A+4 (the Arc); local.tee L; i32.const 1; i32.atomic.rmw.add 0
// and inserts FENCES atomic.fence (straight-line code: no stack check) before the Once load,
// optionally followed by memory.grow(0) (arm "refresh").
// usage: node amp_module.js <in.wasm> <out.wasm|-> <fences> [refresh]
const fs = require("fs");
const [inPath, outPath, fencesArg, refreshArg] = process.argv.slice(2);
const fences = Number(fencesArg);
const bytes = fs.readFileSync(inPath);

const uleb = v => { const o = []; do { let b = v & 127; v = Math.floor(v / 128); if (v) b |= 128; o.push(b); } while (v); return Buffer.from(o); };
function uint(buf, pos) {
  let value = 0; let shift = 0; let at = pos;
  for (;;) { const b = buf[at++]; value += (b & 127) * 2 ** shift; shift += 7; if (!(b & 128)) return [value, at]; }
}

const sections = [];
for (let pos = 8; pos < bytes.length;) {
  const id = bytes[pos];
  const [size, start] = uint(bytes, pos + 1);
  sections.push({ id, bytes: bytes.subarray(start, start + size) });
  pos = start + size;
}

let imports = 0;
{
  const s = sections.find(x => x.id === 2).bytes;
  let [count, at] = uint(s, 0);
  for (let i = 0; i < count; i += 1) {
    let n; [n, at] = uint(s, at); at += n; [n, at] = uint(s, at); at += n;
    const kind = s[at++];
    if (kind === 0) { imports += 1; [, at] = uint(s, at); }
    else if (kind === 1) { at += 1; const f = s[at++]; [, at] = uint(s, at); if (f & 1) [, at] = uint(s, at); }
    else if (kind === 2) { const f = s[at++]; [, at] = uint(s, at); if (f & 1) [, at] = uint(s, at); }
    else if (kind === 3) { at += 2; }
  }
}

const names = new Map();
for (const s of sections) {
  if (s.id !== 0) continue;
  let [n, at] = uint(s.bytes, 0);
  if (s.bytes.subarray(at, at + n).toString() !== "name") continue;
  at += n;
  while (at < s.bytes.length) {
    const sub = s.bytes[at++];
    let size; [size, at] = uint(s.bytes, at);
    if (sub === 1) {
      let [count, p] = uint(s.bytes, at);
      for (let i = 0; i < count; i += 1) {
        let idx; let len; [idx, p] = uint(s.bytes, p); [len, p] = uint(s.bytes, p);
        names.set(idx, s.bytes.subarray(p, p + len).toString()); p += len;
      }
    }
    at += size;
  }
}

// The site, or null, at `p` of a body.
function site(b, p) {
  if (b[p] !== 0x41 || b[p + 1] !== 0 || b[p + 2] !== 0xfe || b[p + 3] !== 0x10 || b[p + 4] !== 2) return null;
  let [a, q] = uint(b, p + 5);
  if (b[q] !== 0x45 || b[q + 1] !== 0x0d || b[q + 2] !== 0 || b[q + 3] !== 0x10) return null;
  let x; [x, q] = uint(b, q + 4);
  if (b[q] !== 0x0b || b[q + 1] !== 0x41 || b[q + 2] !== 0 || b[q + 3] !== 0x28 || b[q + 4] !== 2) return null;
  let a2; [a2, q] = uint(b, q + 5);
  if (a2 !== a + 4 || b[q] !== 0x22) return null;
  let l; [l, q] = uint(b, q + 1);
  if (b[q] !== 0x41 || b[q + 1] !== 1 || b[q + 2] !== 0xfe || b[q + 3] !== 0x1e || b[q + 4] !== 2 || b[q + 5] !== 0) return null;
  return { p, a, x };
}

const code = sections.find(x => x.id === 10);
let [count, at] = uint(code.bytes, 0);
const out = [uleb(count)];
const found = [];
for (let i = 0; i < count; i += 1) {
  let len; let body; [len, body] = uint(code.bytes, at);
  const b = code.bytes.subarray(body, body + len);
  at = body + len;
  const name = names.get(i + imports) || "";
  let hit = null;
  if (name.includes("get_or_init_slow")) {
    for (let p = 0; p < b.length - 12 && !hit; p += 1) hit = site(b, p);
  }
  if (hit) {
    found.push({ index: i + imports, name, ptr: hit.a + 4, once: hit.a, call: hit.x });
    const add = Buffer.alloc(fences * 3);
    for (let k = 0; k < fences; k += 1) { add[k * 3] = 0xfe; add[k * 3 + 1] = 0x03; add[k * 3 + 2] = 0; }
    const refresh = refreshArg === "refresh" ? Buffer.from([0x41, 0, 0x40, 0, 0x1a]) : Buffer.alloc(0);
    const nb = Buffer.concat([b.subarray(0, hit.p), add, refresh, b.subarray(hit.p)]);
    out.push(uleb(nb.length), nb);
  } else {
    out.push(uleb(len), b);
  }
}
for (const f of found) console.log(`site: wasm-function[${f.index}] once=${f.once} ptr=${f.ptr} calls ${f.call} ${f.name}`);
if (found.length === 0) { console.error("no site found"); process.exit(1); }
if (outPath !== "-") {
  code.bytes = Buffer.concat(out);
  const parts = [bytes.subarray(0, 8)];
  for (const s of sections) parts.push(Buffer.from([s.id]), uleb(s.bytes.length), s.bytes);
  const result = Buffer.concat(parts);
  if (!WebAssembly.validate(result)) throw new Error("the rewritten module does not validate");
  fs.writeFileSync(outPath, result);
  console.log(`wrote ${outPath}: ${result.length} bytes, ${fences} fences${refreshArg === "refresh" ? " + memory.grow(0)" : ""}`);
}
