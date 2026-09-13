"use strict";

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// int num_threads = args[kArgNumThreads].AsInteger();
// if (num_threads > kMaxValue) { m_NumThreads = kMaxValue; }
// else { m_NumThreads = num_threads; }
// Host compatibility only: preserve the requested parallel search and every
// original memory operation. No search code or on-disk artifact is changed.
//
// Node v24.21.0, deps/v8/src/wasm/wasm-external-refs.cc:778-810:
//   size_t mem_size = trusted_data->memory_size(mem_index);
//   if (!base::IsInBounds(dst, size, mem_size)) return 0;
// The bulk-operation helpers can see stale per-instance bounds after another
// isolate grows shared memory. wasm-objects.cc:1221-1240 broadcasts growth and
// refreshes the current isolate synchronously, including memory.grow(0).
// Guard using unsigned i64 ranges; refresh only when cached bounds are too small,
// then execute the original operation once. Genuine invalid ranges still trap.
// Binary encoding: WebAssembly Core spec, binary/modules and binary/instructions.
// Accept the wasm32 instruction set emitted by this command target; reject
// unsupported encodings explicitly instead of guessing instruction boundaries.

function uleb(value) {
  const bytes = [];
  do {
    const byte = value & 127;
    value >>>= 7;
    bytes.push(byte | (value ? 128 : 0));
  } while (value);
  return Buffer.from(bytes);
}

class Reader {
  constructor(bytes) { this.bytes = bytes; this.pos = 0; }
  byte() {
    if (this.pos >= this.bytes.length) throw new Error("truncated Wasm encoding");
    return this.bytes[this.pos++];
  }
  take(size) {
    if (size > this.bytes.length - this.pos) throw new Error("truncated Wasm section");
    const result = this.bytes.subarray(this.pos, this.pos + size);
    this.pos += size;
    return result;
  }
  uint() {
    let value = 0;
    for (let shift = 0; shift < 35; shift += 7) {
      const byte = this.byte();
      value += (byte & 127) * 2 ** shift;
      if (!(byte & 128) && value <= 0xffffffff) return value;
    }
    throw new Error("invalid Wasm u32 LEB");
  }
  signed() {
    for (let i = 0; i < 10; i += 1) if (!(this.byte() & 128)) return;
    throw new Error("invalid Wasm signed LEB");
  }
  type() {
    if (![0x7f, 0x7e, 0x7d, 0x7c, 0x7b, 0x70, 0x6f].includes(this.byte())) {
      throw new Error("unsupported Wasm value type");
    }
  }
  memory() {
    if (this.uint() !== 0) throw new Error("only memory index 0 is supported");
  }
  memarg() {
    if (this.uint() & ~0x3f) throw new Error("unsupported Wasm memory addressing");
    this.uint();
  }
  end() {
    if (this.pos !== this.bytes.length) throw new Error("unexpected Wasm section suffix");
  }
}

function sized(bytes) { return Buffer.concat([uleb(bytes.length), bytes]); }
function section(id, bytes) { return Buffer.concat([Buffer.from([id]), sized(bytes)]); }

function instruction(r) {
  const op = r.byte();
  if ([0x02, 0x03, 0x04, 0x41, 0x42, 0xd0].includes(op)) r.signed();
  else if ([0x0c, 0x0d, 0x10, 0x12, 0x14, 0x15, 0xd2, 0xd5, 0xd6].includes(op)
      || (op >= 0x20 && op <= 0x26)) r.uint();
  else if (op === 0x0e) {
    const count = r.uint();
    for (let i = 0; i <= count; i += 1) r.uint();
  } else if (op === 0x11 || op === 0x13) { r.uint(); r.uint(); }
  else if (op === 0x1c) {
    const count = r.uint();
    for (let i = 0; i < count; i += 1) r.type();
  } else if (op >= 0x28 && op <= 0x3e) r.memarg();
  else if (op === 0x3f || op === 0x40) r.memory();
  else if (op === 0x43 || op === 0x44) r.take(op === 0x43 ? 4 : 8);
  else if (op === 0xfc) {
    const sub = r.uint();
    if (sub <= 7) return null;
    if (sub === 8) { r.uint(); r.memory(); }
    else if (sub === 9 || sub === 13 || (sub >= 15 && sub <= 17)) r.uint();
    else if (sub === 10) { r.memory(); r.memory(); return "copy"; }
    else if (sub === 11) { r.memory(); return "fill"; }
    else if (sub === 12 || sub === 14) { r.uint(); r.uint(); }
    else throw new Error(`unsupported Wasm 0xfc opcode ${sub}`);
  } else if (op === 0xfd) {
    const sub = r.uint();
    // Node v24.21.0 deps/v8/src/wasm/wasm-opcodes.h:300-554 defines the
    // supported SIMD encodings. Reject reserved holes even if a future engine
    // assigns an instruction with new immediates to one of them.
    if ([0x9a, 0xa2, 0xa5, 0xa6, 0xaf, 0xb0, 0xb2, 0xb3, 0xb4, 0xbb,
      0xc2, 0xc5, 0xc6, 0xcf, 0xd0, 0xd2, 0xd3, 0xd4, 0xe2, 0xee].includes(sub)) {
      throw new Error(`unsupported Wasm SIMD opcode ${sub}`);
    }
    if (sub <= 11 || sub === 92 || sub === 93) r.memarg();
    else if (sub === 12 || sub === 13) r.take(16);
    else if (sub >= 21 && sub <= 34) r.take(1);
    else if (sub >= 84 && sub <= 91) { r.memarg(); r.take(1); }
    else if (sub > 0x113) throw new Error(`unsupported Wasm SIMD opcode ${sub}`);
  } else if (op === 0xfe) {
    const sub = r.uint();
    if (sub <= 2 || (sub >= 0x10 && sub <= 0x4e)) r.memarg();
    else if (sub === 3) {
      if (r.byte() !== 0) throw new Error("unsupported Wasm atomic fence");
    } else throw new Error(`unsupported Wasm atomic opcode ${sub}`);
  } else if (![0x00, 0x01, 0x05, 0x0b, 0x0f, 0x1a, 0x1b, 0xd1, 0xd3, 0xd4].includes(op)
      && !(op >= 0x45 && op <= 0xc4)) {
    throw new Error(`unsupported Wasm opcode 0x${op.toString(16)}`);
  }
  return null;
}

function guardedBody(copy) {
  // (param dst i32) (param src-or-value i32) (param len i32) (local bytes i64)
  const outside = (index) => [0x20, index, 0xad, 0x20, 2, 0xad, 0x7c, 0x20, 3, 0x56];
  return Buffer.from([
    1, 1, 0x7e,
    0x3f, 0, 0xad, 0x42, 16, 0x86, 0x21, 3, // bytes = u64(memory.size) << 16
    ...outside(0), ...(copy ? [...outside(1), 0x72] : []),
    0x04, 0x40, 0x41, 0, 0x40, 0, 0x1a, 0x0b, // if outside: grow(0); drop
    0x20, 0, 0x20, 1, 0x20, 2,
    0xfc, ...(copy ? [10, 0, 0] : [11, 0]), 0x0b,
  ]);
}

function guardSharedMemory(bytes) {
  if (!WebAssembly.validate(bytes)) throw new Error("invalid input Wasm module");
  const reader = new Reader(bytes);
  const header = reader.take(8);
  const sections = [];
  while (reader.pos < bytes.length) {
    const id = reader.byte();
    sections.push({ id, bytes: reader.take(reader.uint()) });
  }
  const get = (id) => {
    const found = sections.find((s) => s.id === id);
    if (!found) throw new Error(`missing Wasm section ${id}`);
    return found;
  };
  let imports = 0;
  let memories = 0;
  const ir = new Reader(get(2).bytes);
  const importCount = ir.uint();
  for (let i = 0; i < importCount; i += 1) {
    ir.take(ir.uint()); ir.take(ir.uint());
    const kind = ir.byte();
    if (kind === 0) { imports += 1; ir.uint(); }
    else if (kind === 1 || kind === 2) {
      if (kind === 1) ir.type();
      const flags = ir.uint();
      if (kind === 2 ? flags !== 3 : flags > 1) {
        throw new Error("expected a shared wasm32 memory and ordinary tables");
      }
      ir.uint();
      if (flags & 1) ir.uint();
      if (kind === 2) memories += 1;
    } else if (kind === 3) { ir.type(); ir.byte(); }
    else throw new Error(`unsupported Wasm import kind ${kind}`);
  }
  ir.end();
  if (memories !== 1 || sections.some((s) => s.id === 5)) {
    throw new Error("expected exactly one imported shared wasm32 memory");
  }
  const types = get(1);
  const tr = new Reader(types.bytes);
  const typeCount = tr.uint();
  // One type index per entry (recursive/GC type groups are unsupported).
  for (let i = 0; i < typeCount; i += 1) {
    if (tr.byte() !== 0x60) throw new Error("unsupported Wasm function type");
    for (let j = 0; j < 2; j += 1) {
      const count = tr.uint();
      for (let k = 0; k < count; k += 1) tr.type();
    }
  }
  tr.end();
  const functions = get(3);
  const fr = new Reader(functions.bytes);
  const functionCount = fr.uint();
  const functionEntries = functions.bytes.subarray(fr.pos);
  const code = get(10);
  const cr = new Reader(code.bytes);
  if (cr.uint() !== functionCount) throw new Error("Wasm function/code count mismatch");
  const counts = { fill: 0, copy: 0 };
  const bodies = [];
  for (let i = 0; i < functionCount; i += 1) {
    const body = cr.take(cr.uint());
    const r = new Reader(body);
    const locals = r.uint();
    for (let j = 0; j < locals; j += 1) { r.uint(); r.type(); }
    const chunks = [];
    let start = 0;
    while (r.pos < body.length) {
      const before = r.pos;
      const bulk = instruction(r);
      if (bulk) {
        chunks.push(body.subarray(start, before), Buffer.from([0x10]),
          uleb(imports + functionCount + (bulk === "copy" ? 1 : 0)));
        counts[bulk] += 1;
        start = r.pos;
      }
    }
    chunks.push(body.subarray(start));
    bodies.push(sized(Buffer.concat(chunks)));
  }
  cr.end();
  if (counts.fill + counts.copy === 0) return { bytes, ...counts };
  const originalTypes = new Reader(types.bytes);
  originalTypes.uint();
  types.bytes = Buffer.concat([uleb(typeCount + 1), types.bytes.subarray(originalTypes.pos),
    Buffer.from([0x60, 3, 0x7f, 0x7f, 0x7f, 0])]);
  functions.bytes = Buffer.concat([uleb(functionCount + 2), functionEntries,
    uleb(typeCount), uleb(typeCount)]);
  code.bytes = Buffer.concat([uleb(functionCount + 2), ...bodies,
    sized(guardedBody(false)), sized(guardedBody(true))]);
  const guarded = Buffer.concat([header, ...sections.map((s) => section(s.id, s.bytes))]);
  if (!WebAssembly.validate(guarded)) throw new Error("invalid guarded Wasm module");
  return { bytes: guarded, ...counts };
}

module.exports = { guardSharedMemory };
