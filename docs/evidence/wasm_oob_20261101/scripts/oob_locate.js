"use strict";
// Prints the instructions around a byte offset of the host-guarded module (V8's stack offsets
// refer to the guarded module, not the artifact).
// usage: node oob_locate.js <module.wasm> <wasi_shared_memory.js> <hex offset> [more offsets]
const fs = require("fs");
const [wasmPath, guardPath, ...offsets] = process.argv.slice(2);
const { guardSharedMemory } = require(guardPath);
const bytes = guardSharedMemory(fs.readFileSync(wasmPath)).bytes;

function uint(buf, pos) {
  let value = 0; let shift = 0; let at = pos;
  for (;;) {
    const b = buf[at++];
    value += (b & 127) * 2 ** shift;
    shift += 7;
    if (!(b & 128)) return [value, at];
  }
}
function sint(buf, pos) {
  let result = 0n; let shift = 0n; let at = pos; let b;
  do {
    b = buf[at++];
    result |= BigInt(b & 127) << shift;
    shift += 7n;
  } while (b & 128);
  if (b & 64) result |= -1n << shift;
  return [result, at];
}

const PLAIN = {
  0x00: "unreachable", 0x01: "nop", 0x05: "else", 0x0b: "end", 0x0f: "return", 0x1a: "drop", 0x1b: "select",
  0x45: "i32.eqz", 0x46: "i32.eq", 0x47: "i32.ne", 0x48: "i32.lt_s", 0x49: "i32.lt_u", 0x4a: "i32.gt_s", 0x4b: "i32.gt_u",
  0x4c: "i32.le_s", 0x4d: "i32.le_u", 0x4e: "i32.ge_s", 0x4f: "i32.ge_u", 0x50: "i64.eqz", 0x51: "i64.eq", 0x52: "i64.ne",
  0x6a: "i32.add", 0x6b: "i32.sub", 0x6c: "i32.mul", 0x71: "i32.and", 0x72: "i32.or", 0x73: "i32.xor", 0x74: "i32.shl",
  0x75: "i32.shr_s", 0x76: "i32.shr_u", 0x7c: "i64.add", 0x7d: "i64.sub", 0x7e: "i64.mul", 0x83: "i64.and", 0x84: "i64.or",
  0x86: "i64.shl", 0x88: "i64.shr_u", 0xa7: "i32.wrap_i64", 0xac: "i64.extend_i32_s", 0xad: "i64.extend_i32_u",
};
const MEM = {
  0x28: "i32.load", 0x29: "i64.load", 0x2a: "f32.load", 0x2b: "f64.load", 0x2c: "i32.load8_s", 0x2d: "i32.load8_u", 0x2e: "i32.load16_s",
  0x2f: "i32.load16_u", 0x30: "i64.load8_s", 0x31: "i64.load8_u", 0x32: "i64.load16_s", 0x33: "i64.load16_u", 0x34: "i64.load32_s",
  0x35: "i64.load32_u", 0x36: "i32.store", 0x37: "i64.store", 0x38: "f32.store", 0x39: "f64.store", 0x3a: "i32.store8", 0x3b: "i32.store16",
  0x3c: "i64.store8", 0x3d: "i64.store16", 0x3e: "i64.store32",
};
const ATOMIC = {
  0x00: "memory.atomic.notify", 0x01: "memory.atomic.wait32", 0x02: "memory.atomic.wait64",
  0x10: "i32.atomic.load", 0x11: "i64.atomic.load", 0x12: "i32.atomic.load8_u", 0x13: "i32.atomic.load16_u",
  0x14: "i64.atomic.load8_u", 0x15: "i64.atomic.load16_u", 0x16: "i64.atomic.load32_u",
  0x17: "i32.atomic.store", 0x18: "i64.atomic.store", 0x19: "i32.atomic.store8", 0x1a: "i32.atomic.store16",
  0x1b: "i64.atomic.store8", 0x1c: "i64.atomic.store16", 0x1d: "i64.atomic.store32",
  0x1e: "i32.atomic.rmw.add", 0x1f: "i64.atomic.rmw.add", 0x25: "i32.atomic.rmw.sub", 0x26: "i64.atomic.rmw.sub",
  0x2c: "i32.atomic.rmw.and", 0x2d: "i64.atomic.rmw.and", 0x33: "i32.atomic.rmw.or", 0x34: "i64.atomic.rmw.or",
  0x3a: "i32.atomic.rmw.xor", 0x3b: "i64.atomic.rmw.xor", 0x41: "i32.atomic.rmw.xchg", 0x42: "i64.atomic.rmw.xchg",
  0x48: "i32.atomic.rmw.cmpxchg", 0x49: "i64.atomic.rmw.cmpxchg",
};

// One instruction at `pos`: [text, next position]. Throws on an encoding it does not know.
function decode(buf, pos) {
  const op = buf[pos];
  let at = pos + 1;
  if ([0x02, 0x03, 0x04].includes(op)) { const [v, n] = sint(buf, at); return [["block", "loop", "if"][op - 2] + ` ${v}`, n]; }
  if ([0x0c, 0x0d, 0x10, 0x20, 0x21, 0x22, 0x23, 0x24].includes(op)) {
    const [v, n] = uint(buf, at);
    return [{ 0x0c: "br", 0x0d: "br_if", 0x10: "call", 0x20: "local.get", 0x21: "local.set", 0x22: "local.tee", 0x23: "global.get", 0x24: "global.set" }[op] + ` ${v}`, n];
  }
  if (op === 0x0e) {
    let [count, n] = uint(buf, at);
    for (let i = 0; i <= count; i += 1) [, n] = uint(buf, n);
    return ["br_table", n];
  }
  if (op === 0x11) { let [t, n] = uint(buf, at); let tb; [tb, n] = uint(buf, n); return [`call_indirect type ${t} table ${tb}`, n]; }
  if (op === 0x41 || op === 0x42) { const [v, n] = sint(buf, at); return [`${op === 0x41 ? "i32" : "i64"}.const ${v}`, n]; }
  if (op === 0x43 || op === 0x44) return [op === 0x43 ? "f32.const" : "f64.const", at + (op === 0x43 ? 4 : 8)];
  if (op === 0x3f || op === 0x40) return [op === 0x3f ? "memory.size" : "memory.grow", at + 1];
  if (MEM[op]) { let [al, n] = uint(buf, at); let off; [off, n] = uint(buf, n); return [`${MEM[op]} align=${al} offset=${off}`, n]; }
  if (PLAIN[op]) return [PLAIN[op], at];
  if (op >= 0x45 && op <= 0xc4) return [`numeric 0x${op.toString(16)}`, at];
  if (op === 0xfe) {
    let [sub, n] = uint(buf, at);
    if (sub === 3) return ["atomic.fence", n + 1];
    let al; let off; [al, n] = uint(buf, n); [off, n] = uint(buf, n);
    return [`${ATOMIC[sub] || `atomic 0x${sub.toString(16)}`} align=${al} offset=${off}`, n];
  }
  if (op === 0xfc) {
    const [sub, n] = uint(buf, at);
    if (sub === 10) return ["memory.copy", n + 2];
    if (sub === 11) return ["memory.fill", n + 1];
    if (sub <= 7) return [`trunc_sat ${sub}`, n];
    if (sub === 8) { const [, m] = uint(buf, n); return ["memory.init", m + 1]; }
    if (sub === 9) { const [, m] = uint(buf, n); return ["data.drop", m]; }
    throw new Error(`0xfc ${sub}`);
  }
  if (op === 0xfd) {
    const [sub, n] = uint(buf, at);
    if (sub <= 11 || sub === 92 || sub === 93) { let [al, m] = uint(buf, n); let off; [off, m] = uint(buf, m); return [`v128 mem ${sub} align=${al} offset=${off}`, m]; }
    if (sub === 12 || sub === 13) return [`v128 imm ${sub}`, n + 16];
    if (sub >= 21 && sub <= 34) return [`v128 lane ${sub}`, n + 1];
    if (sub >= 84 && sub <= 91) { let [al, m] = uint(buf, n); let off; [off, m] = uint(buf, m); return [`v128 memlane ${sub} align=${al} offset=${off}`, m + 1]; }
    return [`v128 op ${sub}`, n];
  }
  throw new Error(`opcode 0x${op.toString(16)} at ${pos}`);
}

// Walk the sections to the function body containing the offset, then decode that body.
function bodyAround(target) {
  let pos = 8;
  while (pos < bytes.length) {
    const id = bytes[pos];
    const [size, start] = uint(bytes, pos + 1);
    if (id === 10) {
      let [count, at] = uint(bytes, start);
      for (let i = 0; i < count; i += 1) {
        const [len, body] = uint(bytes, at);
        if (target >= body && target < body + len) return { index: i, body, end: body + len };
        at = body + len;
      }
    }
    pos = start + size;
  }
  throw new Error("offset is outside the code section");
}

for (const hex of offsets) {
  const target = parseInt(hex, 16);
  const { index, body, end } = bodyAround(target);
  let [groups, at] = uint(bytes, body);
  for (let i = 0; i < groups; i += 1) { [, at] = uint(bytes, at); at += 1; }
  const list = [];
  while (at < end) {
    const [text, next] = decode(bytes, at);
    list.push({ at, text });
    at = next;
  }
  const hit = list.findIndex(x => x.at === target);
  console.log(`offset 0x${target.toString(16)}: defined function ${index}, body 0x${body.toString(16)}-0x${end.toString(16)}, ${list.length} instructions, exact=${hit >= 0}`);
  const centre = hit >= 0 ? hit : list.findIndex(x => x.at > target);
  for (let i = Math.max(0, centre - 40); i < Math.min(list.length, centre + 8); i += 1) {
    console.log(`${i === hit ? "=>" : "  "} 0x${list[i].at.toString(16)} ${list[i].text}`);
  }
}

// usage addition: env OOB_STATS=1 prints the module's atomic sites and memory.grow sites;
// env OOB_FUNC=<defined index> prints that whole function.
function* bodies() {
  let pos = 8;
  while (pos < bytes.length) {
    const id = bytes[pos];
    const [size, start] = uint(bytes, pos + 1);
    if (id === 10) {
      let [count, at] = uint(bytes, start);
      for (let i = 0; i < count; i += 1) {
        const [len, body] = uint(bytes, at);
        yield { index: i, body, end: body + len };
        at = body + len;
      }
    }
    pos = start + size;
  }
}
function instructions(f) {
  let [groups, at] = uint(bytes, f.body);
  for (let i = 0; i < groups; i += 1) { [, at] = uint(bytes, at); at += 1; }
  const list = [];
  while (at < f.end) { const [text, next] = decode(bytes, at); list.push({ at, text }); at = next; }
  return list;
}
if (process.env.OOB_FUNC) {
  for (const f of bodies()) {
    if (f.index === Number(process.env.OOB_FUNC)) for (const x of instructions(f)) console.log(`0x${x.at.toString(16)} ${x.text}`);
  }
}
if (process.env.OOB_STATS) {
  const perOp = new Map(); const pairs = new Set(); const grow = []; let sites = 0; let constLoad = 0; let functions = 0; let withAtomic = 0;
  for (const f of bodies()) {
    functions += 1;
    const list = instructions(f);
    let any = false;
    list.forEach((x, i) => {
      if (x.text === "memory.grow") grow.push(f.index);
      if (/atomic\.|^memory\.atomic/.test(x.text) && !/fence/.test(x.text)) {
        sites += 1; any = true;
        const m = x.text.match(/^(\S+) align=\d+ offset=(\d+)/);
        perOp.set(m[1], (perOp.get(m[1]) || 0) + 1);
        pairs.add(`${m[1]}/${m[2]}`);
        if (/load/.test(m[1]) && i > 0 && /^i32\.const/.test(list[i - 1].text)) constLoad += 1;
      }
    });
    if (any) withAtomic += 1;
  }
  console.log(`functions ${functions}, with atomics ${withAtomic}, atomic sites ${sites}, distinct (op, offset) ${pairs.size}, loads of a constant address ${constLoad}`);
  console.log([...perOp].sort((a, b) => b[1] - a[1]).map(([k, v]) => `${k}:${v}`).join(" "));
  console.log(`memory.grow sites in defined functions: ${grow.join(",")}`);
}
