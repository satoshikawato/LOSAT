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
//
// Atomic memory instructions read the same stale bound. Node v26.8.2 (V8 14.6),
// deps/v8/src/wasm/baseline/liftoff-compiler.cc: every atomic access calls
// BoundsCheckMem(..., kDoForceCheck, kCheckAlignment), a compare with the instance's
// memory size, where a plain load or store relies on the trap handler; Turboshaft does the
// same for wait and notify (turboshaft-graph-interface.cc, kNeedsBoundsCheck). Another thread's
// growth reaches the instance's size only when its isolate handles the growth interrupt, at a
// stack check or a tier-up check, and Liftoff::Loop has no stack check under dynamic tiering.
// So a thread that reads a pointer into memory grown by another thread can trap on it although
// the memory is there: rayon's first Stealer::steal at Arc::clone of crossbeam's Collector, in
// dc-megablast at 16 threads, 1 of the 483 threaded runs on record (2026-11-01,
// docs/evidence/wasm_oob_20261101/).
// Each atomic instruction becomes a call of a helper with the same bound check as the bulk
// helpers, one helper per opcode and memarg, whose only effect is that refresh (on V8 the
// helper's entry stack check usually handles the interrupt first; the compare does not depend
// on that). The instruction itself runs unchanged, once.
//
// V8 compares plain multi-byte stores too on ARM64 outside macOS, in both tiers
// (wasm/compilation-environment.h, kPartialOOBWritesAreNoops): real blastp at 16 threads trapped
// in dlmalloc in 64 and 100 of 100 runs on ubuntu-24.04-arm (Node 26.8.2, 24.21) with only the bulk
// instructions guarded, 11 and 29 with the atomics guarded too. A module has
// about 200,000 such stores, too many to wrap. But a thread learns an address in pages that
// another thread grew only from what that thread wrote, and an atomic read (a lock, or the atomic
// that published the address) is what orders the two. So memory.grow also raises GROWN_PAGES,
// before the grower can publish anything, and each atomic instruction that reads memory compares
// memory.size with GROWN_PAGES afterwards: smaller means that this instance's bound is stale, and
// memory.grow(0) refreshes it before the thread goes on.
//
// A wait (memory.atomic.wait32, wait64) without a timeout never returns in WebKit once the page
// ends the waiting thread's worker. JavaScriptCore wakes the waiter for the termination, but the wait
// loop re-checks a flag that only the handling of the termination sets, which a waiting thread cannot
// do (WebKit d8fa5ad85d, Playwright's WebKit 2311: WaiterListManager::waitSyncImpl loops while
// !vm.hasTerminationRequest(); VM::notifyNeedTermination only fires the trap). The worker, its
// instance and the shared memory then live until the page closes. A page can reserve a memory of
// 4 GiB about 15 times, and each cancel of a threaded search whose threads were waiting kept one, until
// no threaded search could start (web/app/tests/e2e/worker-end.spec.ts). A wait that has a timeout
// ends, and a thread that calls a host function is ended there; but a finite timeout costs each wait
// a kernel timer (about 1 us in all three browsers, ct3 ping-pong), so waits keep their own timeout
// while the memory's workers live. The helper of a wait records the address it waits at, plus 1 (so
// that 0 means a free cell), in WAIT_TABLE (one of WAIT_SLOTS cells, chosen from the address) and
// clears it when the wait returns.
// When the host ends the memory's workers (the memory keeper, ThreadHost.terminate: endWaits in
// web/app/src/infra/reactor/shared-memory.ts) it sets WAITS_ENDING and notifies every recorded address,
// again each millisecond, until the table is empty (or nothing changed in it and no notify woke a
// thread for a second: a cell that stays set without a waiter, after an overflow of the main thread's
// stack into these cells). From then on, and while the table is full, the helper waits in slices of
// at most WAIT_SLICE_NS and calls the module's import wasi_snapshot_preview1.sched_yield between them;
// the browser ends the worker in that call. The helper re-reads WAITS_ENDING after it recorded its
// address, so a notify that comes before its wait is repeated by the host (it notifies until the table
// is empty). Before it records the address it loads from it, which traps as the wait would (bounds,
// alignment), so no address stays recorded after a trap. A cell taken after a trap cannot happen; the
// cells, like GROWN_PAGES, are reached by the main thread's stack only at the end of an overflow. The helper returns what the instruction
// returns (0 woken, 1 not equal, 2 timed out; the host's notify is a wake that the guest's futexes
// treat as spurious), and waits at least a finite timeout in full. A module that does not import
// sched_yield (the hosts' test fixtures) has no host function to call, and its wait helpers are the
// single wait, as before; the web build refuses a threaded module without it (web/app/build/reactors.ts).
//
// memory.init stays as it is: only __wasm_init_memory (the module's start) and __wasm_init_tls
// use it, the latter in the spawning thread, called from pthread_create right after the malloc
// of the block, with function entries (stack checks) in between.
// Binary encoding: WebAssembly Core spec, binary/modules and binary/instructions.
// Accept the wasm32 instruction set emitted by this command target; reject
// unsupported encodings explicitly instead of guessing instruction boundaries.
//
// The guarded module carries the custom section GUARD_SECTION (the guard's version and
// counts). A module that has it is returned as it is, so a host may guard every module it
// compiles and a module guarded at build time is not guarded twice. The build guards the
// threaded artifacts with the command line at the end of this file (`guard`, `verify`), so
// that every host, browsers included, gets guarded bytes (release v0.2.0 DECISIONS 52).

const crypto = require("node:crypto");
const fs = require("node:fs");

const GUARD_SECTION = "losat.shared_memory_guard";
const GUARD_VERSION = 4;
// The longest slice of a wait once the host ends the workers, in nanoseconds: how long such a
// worker may still hold the shared memory.
const WAIT_SLICE_NS = 100_000_000;
// Nonzero once the host ends the memory's workers (endWaits): waits wait in slices (the header).
const WAITS_ENDING = 20;
// WAIT_SLOTS i32 cells from WAIT_TABLE: the address of a wait in progress plus 1, 0 for none (the header).
const WAIT_TABLE = 64;
const WAIT_SLOTS = 128;
// Address of the largest memory size, in pages, that a thread of the guest grew the memory to
// (0 = none yet). WAITS_ENDING and WAIT_TABLE lie below 1024 too. Addresses below 1024 are unused (wasm-ld's default global base) or, with
// --stack-first as rustc links this target, the last bytes of the main thread's stack, reached
// only at the end of an overflow (docs/evidence/wasm_oob_20261101/). A wrong value here costs
// one more memory.grow(0), or leaves the bound check of the helpers as the only refresh.
const GROWN_PAGES = 16;
const WASM_HEADER = Buffer.from([0, 0x61, 0x73, 0x6d, 1, 0, 0, 0]);

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
    const align = this.uint();
    if (align & ~0x3f) throw new Error("unsupported Wasm memory addressing");
    return { align, offset: this.uint() };
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
  else if (op === 0x3f) r.memory();
  else if (op === 0x40) { r.memory(); return "grow"; }
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
    if (sub <= 2 || (sub >= 0x10 && sub <= 0x4e)) return { sub, ...r.memarg() };
    if (sub === 3) {
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

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host compatibility only. An atomic instruction's operands after the address, its result and
// the bytes it accesses (WebAssembly threads proposal, binary encoding: 0x10 loads, 0x17 stores,
// 0x1e-0x41 read-modify-writes, 0x48 compare-exchanges, seven widths each).
const ATOMIC_WIDTHS = [[0x7f, 4], [0x7e, 8], [0x7f, 1], [0x7f, 2], [0x7e, 1], [0x7e, 2], [0x7e, 4]];
function atomicShape(sub) {
  if (sub === 0) return { values: [0x7f], result: [0x7f], size: 4 }; // memory.atomic.notify
  if (sub === 1) return { values: [0x7f, 0x7e], result: [0x7f], size: 4 }; // memory.atomic.wait32
  if (sub === 2) return { values: [0x7e, 0x7e], result: [0x7f], size: 8 }; // memory.atomic.wait64
  const family = Math.floor((sub - 0x10) / 7);
  const [type, size] = ATOMIC_WIDTHS[(sub - 0x10) % 7];
  if (family === 0) return { values: [], result: [type], size };
  if (family === 1) return { values: [type], result: [], size };
  return { values: family === 8 ? [type, type] : [type], result: [type], size };
}

function atomicType({ sub }) {
  const { values, result } = atomicShape(sub);
  return Buffer.from([0x60, values.length + 1, 0x7f, ...values, result.length, ...result]);
}

function signedLeb(value) {
  const bytes = [];
  for (;;) {
    const byte = value % 128;
    value = Math.floor(value / 128);
    if (value === 0 && !(byte & 64)) return [...bytes, byte];
    bytes.push(byte | 128);
  }
}

// The atomic instructions that read memory: all but the stores and notify.
const atomicAcquires = (sub) => sub !== 0 && !(sub >= 0x17 && sub <= 0x1d);

const isWait = (sub) => sub === 1 || sub === 2;

// The bound check of an atomic helper: bytes = u64(memory.size) << 16, then memory.grow(0) when
// address + offset + size is past them.
function boundCheck(bytes, offset, size) {
  return [
    0x3f, 0, 0xad, 0x42, 16, 0x86, 0x21, bytes,
    0x20, 0, 0xad, 0x42, ...signedLeb(offset + size), 0x7c, 0x20, bytes, 0x56,
    0x04, 0x40, 0x41, 0, 0x40, 0, 0x1a, 0x0b,
  ];
}

// After an atomic instruction that reads memory: memory.grow(0) when memory.size < GROWN_PAGES.
const ACQUIRE_CHECK = [0x3f, 0, 0x41, GROWN_PAGES, 0xfe, 0x10, 2, 0, 0x49, 0x04, 0x40, 0x41, 0, 0x40, 0, 0x1a, 0x0b];

// i32.const of a 32-bit value (signed LEB128 of its two's complement).
function i32Const(value) {
  let v = value | 0;
  const bytes = [];
  for (;;) {
    const byte = v & 127;
    v >>= 7;
    if ((v === 0 && !(byte & 64)) || (v === -1 && (byte & 64))) return [0x41, ...bytes, byte];
    bytes.push(byte | 128);
  }
}

// A wait (the header): the bound check, a load at the address (it traps as the wait would), then
//   loop: if WAITS_ENDING or no free cell (sliced):
//           slice = timeout < 0 || timeout > SLICE ? SLICE : timeout
//           result = wait(address, expected, slice); result != 2 (woken, not equal): done
//           a finite timeout: timeout -= slice; timeout <= 0: done (result is 2)
//           sched_yield(); again
//         record address + 1 in a free cell (from cell (address >> 2) & (SLOTS - 1) on); none free:
//           sliced = 1; again
//         WAITS_ENDING now: clear the cell; again
//         result = wait(address, expected, timeout); clear the cell; done
// then the acquire check, and the result.
function waitBody({ sub, align, offset }, yieldIndex) {
  // (param address i32) (param expected) (param timeout i64) (local bytes i64) (local result i32)
  // (local slice i64) (local cell i32) (local probe i32) (local at = address + offset + 1, i32) (local sliced i32)
  const { size } = atomicShape(sub);
  const slice = signedLeb(WAIT_SLICE_NS);
  const wait = [0xfe, ...uleb(sub), ...uleb(align), ...uleb(offset)];
  const ending = [0x41, 0, 0xfe, 0x10, 2, ...uleb(WAITS_ENDING)];
  const clearCell = [0x20, 6, 0x41, 2, 0x74, 0x41, 0, 0xfe, 0x17, 2, ...uleb(WAIT_TABLE)];
  return Buffer.from([
    4, 1, 0x7e, 1, 0x7f, 1, 0x7e, 4, 0x7f,
    ...boundCheck(3, offset, size),
    0x20, 0, 0xfe, sub === 1 ? 0x10 : 0x11, ...uleb(align), ...uleb(offset), 0x1a,
    0x20, 0, ...i32Const(offset), 0x6a, 0x41, 1, 0x6a, 0x21, 8,
    0x02, 0x40, 0x03, 0x40,
    0x20, 9, ...ending, 0x72, 0x04, 0x40,
    0x20, 2, 0x42, 0, 0x53, 0x20, 2, 0x42, ...slice, 0x55, 0x72,
    0x04, 0x7e, 0x42, ...slice, 0x05, 0x20, 2, 0x0b, 0x21, 5,
    0x20, 0, 0x20, 1, 0x20, 5, ...wait,
    0x22, 4, 0x41, 2, 0x47, 0x0d, 2,
    0x20, 2, 0x42, 0, 0x59, 0x04, 0x40,
    0x20, 2, 0x20, 5, 0x7d, 0x22, 2, 0x42, 0, 0x57, 0x0d, 3,
    0x0b,
    0x10, ...uleb(yieldIndex), 0x1a, 0x0c, 1,
    0x0b,
    0x41, 0, 0x21, 7,
    0x02, 0x40, 0x03, 0x40,
    0x20, 7, ...i32Const(WAIT_SLOTS), 0x4f, 0x0d, 1,
    0x20, 8, 0x41, 2, 0x76, 0x20, 7, 0x6a, ...i32Const(WAIT_SLOTS - 1), 0x71, 0x21, 6,
    0x20, 6, 0x41, 2, 0x74, 0x41, 0, 0x20, 8, 0xfe, 0x48, 2, ...uleb(WAIT_TABLE),
    0x45, 0x0d, 1,
    0x20, 7, 0x41, 1, 0x6a, 0x21, 7, 0x0c, 0,
    0x0b, 0x0b,
    0x20, 7, ...i32Const(WAIT_SLOTS), 0x4f, 0x04, 0x40, 0x41, 1, 0x21, 9, 0x0c, 1, 0x0b,
    ...ending, 0x04, 0x40, ...clearCell, 0x0c, 1, 0x0b,
    0x20, 0, 0x20, 1, 0x20, 2, ...wait, 0x21, 4,
    ...clearCell,
    0x0c, 1,
    0x0b, 0x0b,
    ...ACQUIRE_CHECK,
    0x20, 4, 0x0b,
  ]);
}

function atomicBody({ sub, align, offset }, yieldIndex) {
  if (isWait(sub) && yieldIndex !== undefined) return waitBody({ sub, align, offset }, yieldIndex);
  // (param address i32) (param values...) (local bytes i64) (local result)
  const { values, result, size } = atomicShape(sub);
  const bytes = values.length + 1;
  const kept = bytes + 1;
  const after = atomicAcquires(sub) ? ACQUIRE_CHECK : [];
  return Buffer.from([
    ...(result.length ? [2, 1, 0x7e, 1, result[0]] : [1, 1, 0x7e]),
    ...boundCheck(bytes, offset, size),
    0x20, 0, ...values.flatMap((_, i) => [0x20, i + 1]),
    0xfe, ...uleb(sub), ...uleb(align), ...uleb(offset),
    ...(result.length ? [0x21, kept] : []),
    ...after,
    ...(result.length ? [0x20, kept] : []),
    0x0b,
  ]);
}

// The original memory.grow, then GROWN_PAGES = max(GROWN_PAGES, memory.size) when it succeeded.
const GROW_TYPE = Buffer.from([0x60, 1, 0x7f, 1, 0x7f]);
function growBody() {
  // (param pages i32) (local old i32) (local seen i32)
  return Buffer.from([
    1, 2, 0x7f,
    0x20, 0, 0x40, 0, 0x22, 1, 0x41, 0x7f, 0x47, // (old = memory.grow(pages)) != -1
    0x04, 0x40,
    0x03, 0x40,
    0x41, GROWN_PAGES, 0xfe, 0x10, 2, 0, 0x22, 2, // seen = GROWN_PAGES
    0x3f, 0, 0x49, // seen < memory.size
    0x04, 0x40,
    0x41, GROWN_PAGES, 0x20, 2, 0x3f, 0, 0xfe, 0x48, 2, 0, 0x20, 2, 0x47, 0x0d, 1, // cmpxchg failed: again
    0x0b,
    0x0b,
    0x0b,
    0x20, 1, 0x0b,
  ]);
}

function sectionsOf(bytes) {
  const reader = new Reader(bytes);
  const header = reader.take(8);
  if (!WASM_HEADER.equals(Buffer.from(header))) throw new Error("not a Wasm 1.0 module");
  const sections = [];
  while (reader.pos < bytes.length) {
    const id = reader.byte();
    sections.push({ id, bytes: reader.take(reader.uint()) });
  }
  return { header, sections };
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host compatibility only. The record of a module that guardSharedMemory wrote:
// { version, fill, copy, grow, atomic, atomic_helpers }, or null. A module whose encoding cannot be read has none
// (guardSharedMemory then rejects it as invalid).
function guardRecord(bytes) {
  let payload = null;
  try {
    for (const s of sectionsOf(bytes).sections) {
      if (s.id !== 0) continue;
      const r = new Reader(s.bytes);
      if (Buffer.from(r.take(r.uint())).toString("utf8") === GUARD_SECTION) payload = s.bytes.subarray(r.pos);
    }
  } catch {
    return null;
  }
  if (payload === null) return null;
  let record;
  try { record = JSON.parse(Buffer.from(payload).toString("utf8")); } catch { record = null; }
  const count = (value) => Number.isSafeInteger(value) && value >= 0;
  if (record === null || typeof record !== "object" || record.version !== GUARD_VERSION
      || !count(record.fill) || !count(record.copy) || !count(record.grow) || !count(record.atomic) || !count(record.atomic_helpers)
      || record.fill + record.copy + record.grow + record.atomic === 0 || record.atomic_helpers > record.atomic) {
    throw new Error(`unsupported ${GUARD_SECTION} record: ${Buffer.from(payload).toString("utf8").slice(0, 80)}`);
  }
  return { version: record.version, fill: record.fill, copy: record.copy, grow: record.grow, atomic: record.atomic,
    atomic_helpers: record.atomic_helpers };
}

// The sections that the guard rewrites, read once for guarding and for verifying.
function layout(bytes) {
  const { header, sections } = sectionsOf(bytes);
  const get = (id) => {
    const found = sections.find((s) => s.id === id);
    if (!found) throw new Error(`missing Wasm section ${id}`);
    return found;
  };
  let imports = 0;
  let memories = 0;
  // "module.name" of each imported function, by function index, and its type index.
  const importedFunctions = [];
  const ir = new Reader(get(2).bytes);
  const importCount = ir.uint();
  for (let i = 0; i < importCount; i += 1) {
    const module = Buffer.from(ir.take(ir.uint())).toString("utf8");
    const field = Buffer.from(ir.take(ir.uint())).toString("utf8");
    const kind = ir.byte();
    if (kind === 0) { imports += 1; importedFunctions.push({ name: `${module}.${field}`, type: ir.uint() }); }
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
  const typeEntries = [];
  for (let i = 0; i < typeCount; i += 1) {
    const at = tr.pos;
    if (tr.byte() !== 0x60) throw new Error("unsupported Wasm function type");
    for (let j = 0; j < 2; j += 1) {
      const count = tr.uint();
      for (let k = 0; k < count; k += 1) tr.type();
    }
    typeEntries.push(types.bytes.subarray(at, tr.pos));
  }
  tr.end();
  const functions = get(3);
  const fr = new Reader(functions.bytes);
  const functionCount = fr.uint();
  const functionEntries = functions.bytes.subarray(fr.pos);
  const functionTypes = [];
  for (let i = 0; i < functionCount; i += 1) functionTypes.push(fr.uint());
  fr.end();
  const code = get(10);
  const cr = new Reader(code.bytes);
  if (cr.uint() !== functionCount) throw new Error("Wasm function/code count mismatch");
  const bodies = [];
  for (let i = 0; i < functionCount; i += 1) bodies.push(cr.take(cr.uint()));
  cr.end();
  return { header, sections, types, typeCount, typeEntries, functions, functionCount, functionEntries,
    functionTypes, code, bodies, imports, importedFunctions };
}

// The function index of the module's import wasi_snapshot_preview1.sched_yield, which the helper
// of a wait calls between its slices (the header), or undefined if the module does not import it.
const SCHED_YIELD = "wasi_snapshot_preview1.sched_yield";
function schedYield(m) {
  const index = m.importedFunctions.findIndex((f) => f.name === SCHED_YIELD);
  if (index < 0) return undefined;
  if (!Buffer.from(m.typeEntries[m.importedFunctions[index].type]).equals(Buffer.from([0x60, 0, 1, 0x7f]))) {
    throw new Error(`${SCHED_YIELD} is not a function () -> i32`);
  }
  return index;
}

// Calls `visit(guarded, before, after)` for every instruction of a function body, after its
// locals; `guarded` is "fill", "copy", "grow", an atomic instruction's { sub, align, offset }, or null.
function walk(body, visit) {
  const r = new Reader(body);
  const locals = r.uint();
  for (let j = 0; j < locals; j += 1) { r.uint(); r.type(); }
  while (r.pos < body.length) {
    const before = r.pos;
    visit(instruction(r), before, r.pos);
  }
}

const GUARD_TYPE = Buffer.from([0x60, 3, 0x7f, 0x7f, 0x7f, 0]);

function guardSharedMemory(bytes) {
  // A module that this guard wrote is returned as it is (nothing guarded by this call).
  const record = guardRecord(bytes);
  if (record !== null) return { bytes, fill: 0, copy: 0, grow: 0, atomic: 0, record };
  if (!WebAssembly.validate(bytes)) throw new Error("invalid input Wasm module");
  const m = layout(bytes);
  const { types, typeCount, functions, functionCount, functionEntries, code, imports } = m;
  const counts = { fill: 0, copy: 0, grow: 0, atomic: 0 };
  // The defined functions keep their indices; then fill, copy, grow, and one helper per atomic
  // opcode and memarg, in the order of their first use.
  const helpers = new Map();
  const bodies = m.bodies.map((body) => {
    const chunks = [];
    let start = 0;
    walk(body, (guarded, before, after) => {
      if (!guarded) return;
      let target;
      if (typeof guarded === "string") {
        target = imports + functionCount + ["fill", "copy", "grow"].indexOf(guarded);
        counts[guarded] += 1;
      } else {
        const key = `${guarded.sub}/${guarded.align}/${guarded.offset}`;
        if (!helpers.has(key)) helpers.set(key, { ...guarded, index: imports + functionCount + 3 + helpers.size });
        target = helpers.get(key).index;
        counts.atomic += 1;
      }
      chunks.push(body.subarray(start, before), Buffer.from([0x10]), uleb(target));
      start = after;
    });
    chunks.push(body.subarray(start));
    return sized(Buffer.concat(chunks));
  });
  if (counts.fill + counts.copy + counts.grow + counts.atomic === 0) return { bytes, ...counts, record: null };
  const originalTypes = new Reader(types.bytes);
  originalTypes.uint();
  const atomicHelpers = [...helpers.values()];
  const yieldIndex = atomicHelpers.some((helper) => isWait(helper.sub)) ? schedYield(m) : undefined;
  const newTypes = new Map();
  const helperTypes = atomicHelpers.map((helper) => {
    const type = atomicType(helper);
    const key = type.toString("hex");
    if (!newTypes.has(key)) newTypes.set(key, { type, index: typeCount + 2 + newTypes.size });
    return newTypes.get(key).index;
  });
  types.bytes = Buffer.concat([uleb(typeCount + 2 + newTypes.size), types.bytes.subarray(originalTypes.pos), GUARD_TYPE, GROW_TYPE,
    ...[...newTypes.values()].map((entry) => entry.type)]);
  const total = functionCount + 3 + atomicHelpers.length;
  functions.bytes = Buffer.concat([uleb(total), functionEntries,
    uleb(typeCount), uleb(typeCount), uleb(typeCount + 1), ...helperTypes.map(uleb)]);
  code.bytes = Buffer.concat([uleb(total), ...bodies,
    sized(guardedBody(false)), sized(guardedBody(true)), sized(growBody()), ...atomicHelpers.map((helper) => sized(atomicBody(helper, yieldIndex)))]);
  const written = { version: GUARD_VERSION, ...counts, atomic_helpers: atomicHelpers.length };
  const marker = section(0, Buffer.concat([sized(Buffer.from(GUARD_SECTION)), Buffer.from(JSON.stringify(written))]));
  const guarded = Buffer.concat([m.header, ...m.sections.map((s) => section(s.id, s.bytes)), marker]);
  if (!WebAssembly.validate(guarded)) throw new Error("invalid guarded Wasm module");
  return { bytes: guarded, ...counts, record: written };
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host compatibility only. Throws unless `bytes` is a module that guardSharedMemory wrote:
// its record, the helpers last with the guard's types and bodies (fill, copy, then the atomic
// ones), no memory.fill, memory.copy or atomic memory instruction anywhere else, and as many
// direct calls of the helpers as the record counts. So no unguarded instruction remains; which
// helper replaced which instruction is not checked (a build that needs that guards the input
// again and compares: prepare_release_candidate_v020.py). A host trusts the record; the build
// checks the artifact it ships with this.
function verifyGuarded(bytes) {
  const record = guardRecord(bytes);
  if (record === null) throw new Error(`the module has no ${GUARD_SECTION} section`);
  if (!WebAssembly.validate(bytes)) throw new Error("invalid Wasm module");
  const m = layout(bytes);
  const n = m.functionCount;
  const fill = n - 3 - record.atomic_helpers;
  if (fill < 0) throw new Error("the guard's helpers are missing");
  const typeOf = (i) => Buffer.from(m.typeEntries[m.functionTypes[i]]);
  if (!GUARD_TYPE.equals(typeOf(fill)) || !GUARD_TYPE.equals(typeOf(fill + 1))) {
    throw new Error("the guard's helpers do not have the guard's type");
  }
  if (!Buffer.from(m.bodies[fill]).equals(guardedBody(false)) || !Buffer.from(m.bodies[fill + 1]).equals(guardedBody(true))) {
    throw new Error("the guard's helpers are not the guard's bodies");
  }
  if (!GROW_TYPE.equals(typeOf(fill + 2))) throw new Error("the guard's helpers do not have the guard's type");
  if (!Buffer.from(m.bodies[fill + 2]).equals(growBody())) throw new Error("the guard's helpers are not the guard's bodies");
  for (let i = fill + 3; i < n; i += 1) {
    const found = [];
    walk(m.bodies[i], (guarded) => { if (guarded && typeof guarded !== "string") found.push(guarded); });
    // The helper's operation: its first atomic instruction, or the wait of a wait helper (whose first
    // atomic instructions are the load at the wait's address and the load of WAITS_ENDING).
    const op = found.find((g) => isWait(g.sub)) ?? found[0];
    const yieldIndex = op !== undefined && isWait(op.sub) ? schedYield(m) : undefined;
    if (op === undefined || !Buffer.from(m.bodies[i]).equals(atomicBody(op, yieldIndex))) {
      throw new Error("the guard's helpers are not the guard's bodies");
    }
    if (!atomicType(op).equals(typeOf(i))) throw new Error("the guard's helpers do not have the guard's type");
  }
  const calls = [0, 0, 0, 0];
  for (let i = 0; i < fill; i += 1) {
    const body = m.bodies[i];
    walk(body, (guarded, before) => {
      if (typeof guarded === "string") throw new Error(`function ${m.imports + i} has an unguarded memory.${guarded}`);
      if (guarded) throw new Error(`function ${m.imports + i} has an unguarded atomic instruction 0xfe ${guarded.sub}`);
      if (body[before] !== 0x10) return;
      const r = new Reader(body);
      r.pos = before + 1;
      const target = r.uint() - m.imports - fill;
      if (target >= 0) calls[Math.min(target, 3)] += 1;
    });
  }
  if (calls[0] !== record.fill || calls[1] !== record.copy || calls[2] !== record.grow || calls[3] !== record.atomic) {
    throw new Error(`the record counts ${record.fill} fills, ${record.copy} copies, ${record.grow} grows and ${record.atomic} atomics, the code calls the helpers ${calls.join(", ")} times`);
  }
  return record;
}

// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host compatibility only. Whether the guard's waits in `bytes` can end with their worker (the
// header): the module imports wasi_snapshot_preview1.sched_yield, () -> i32.
function slicesWaits(bytes) {
  return schedYield(layout(bytes)) !== undefined;
}

const sha256 = (bytes) => crypto.createHash("sha256").update(bytes).digest("hex");

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); return status;
// Build step: `guard <input.wasm> <output.wasm>` writes the guarded module (through a
// temporary file and a rename) and prints its record as JSON; `verify <file.wasm>` checks a
// guarded module. Either exits 1 with the reason, and `guard` writes nothing, when the
// module cannot be guarded (an encoding the reader declines), has nothing to guard, or does
// not verify: a build never ships an unguarded threaded module without noticing.
function main(argv) {
  const [command, input, output] = argv;
  if (command === "guard" && input && output && argv.length === 3) {
    const bytes = fs.readFileSync(input);
    const result = guardSharedMemory(bytes);
    if (result.record === null) {
      throw new Error(`${input}: no memory.fill, memory.copy, memory.grow or atomic instruction to guard (is it a threaded module?)`);
    }
    verifyGuarded(result.bytes);
    const temp = `${output}.${process.pid}.${crypto.randomBytes(6).toString("hex")}.tmp`;
    try {
      fs.writeFileSync(temp, result.bytes, { flag: "wx" });
      fs.renameSync(temp, output);
    } catch (err) {
      fs.rmSync(temp, { force: true });
      throw err;
    }
    return { input_sha256: sha256(bytes), output_sha256: sha256(result.bytes), output_bytes: result.bytes.length,
      ...result.record, guarded_before: result.fill + result.copy + result.grow + result.atomic === 0, guard_sha256: sha256(fs.readFileSync(__filename)) };
  }
  if (command === "verify" && input && argv.length === 2) {
    const bytes = fs.readFileSync(input);
    return { sha256: sha256(bytes), output_bytes: bytes.length, ...verifyGuarded(bytes) };
  }
  throw new Error("usage: wasi_shared_memory.js guard <input.wasm> <output.wasm> | verify <file.wasm>");
}

if (require.main === module) {
  try {
    process.stdout.write(`${JSON.stringify(main(process.argv.slice(2)))}\n`);
  } catch (err) {
    process.stderr.write(`wasi_shared_memory.js: ${err.message}\n`);
    process.exitCode = 1;
  }
}

module.exports = { guardSharedMemory, guardRecord, verifyGuarded, slicesWaits, GUARD_SECTION, WAITS_ENDING, WAIT_TABLE, WAIT_SLOTS };
