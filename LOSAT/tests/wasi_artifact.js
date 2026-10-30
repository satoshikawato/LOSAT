"use strict";
const crypto = require("node:crypto");

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); return status;
// WASI ABI boundary: command uses _start; reactor uses _initialize once.
// Inspect actual imports/exports and binary memory limits, never the filename.
// Module validation and ownership are Wasm-host concerns. Retain the validated
// module for execution; each invocation still creates its own WASI instance.
// A caller that has compiled a module with the same imports and exports from these bytes (the
// threaded host's guarded module) passes it, and the bytes' SHA-256 if it has it, to skip a
// second compile.
function prepareArtifact(bytes, expected, module = new WebAssembly.Module(bytes), sha256 = null) {
  const imports = WebAssembly.Module.imports(module);
  const exports = WebAssembly.Module.exports(module);
  const names = new Set(exports.filter(x => x.kind === "function").map(x => x.name));
  const memories = [];
  let pos = 8;
  const uint = () => {
    let v = 0, shift = 0, b;
    do { b = bytes[pos++]; v += (b & 127) * 2 ** shift; shift += 7; } while (b & 128);
    return v;
  };
  const string = () => { const n = uint(); const s = bytes.subarray(pos, pos+n).toString(); pos += n; return s; };
  const limits = () => {
    const flags = uint(), initial = uint(), maximum = flags & 1 ? uint() : null;
    if (flags & ~3) throw new Error("unsupported WASI memory limits");
    return { initial, maximum, shared: Boolean(flags & 2) };
  };
  while (pos < bytes.length) {
    const id = bytes[pos++], size = uint(), end = pos + size;
    if (id === 2) {
      const count = uint();
      for (let i = 0; i < count; i++) {
        const mod = string(), name = string(), kind = bytes[pos++];
        if (kind === 0) uint();
        else if (kind === 1) { pos++; limits(); }
        else if (kind === 2) memories.push({ module: mod, name, imported: true, ...limits() });
        else if (kind === 3) pos += 2;
        else throw new Error("unsupported WASI import kind");
      }
    } else if (id === 5) {
      const count = uint();
      for (let i = 0; i < count; i++) memories.push({ imported: false, ...limits() });
    }
    pos = end;
  }
  const threaded = imports.some(x => x.module === "wasi" && x.name === "thread-spawn" && x.kind === "function");
  const command = names.has("_start"), reactor = names.has("_initialize");
  if (command === reactor) throw new Error("artifact must have exactly one of _start or _initialize");
  if (reactor && !["losat_web_run_pair", "losat_web_alloc", "losat_web_dealloc", "losat_web_result_ptr", "losat_web_result_len", "losat_web_error_ptr", "losat_web_error_len"].every(x => names.has(x))) {
    throw new Error("reactor is missing direct API exports");
  }
  if (memories.length !== 1 || !exports.some(x => x.name === "memory" && x.kind === "memory")) throw new Error("artifact must export its single memory");
  if (threaded && (!names.has("wasi_thread_start") || !memories[0].shared || !memories[0].imported || memories[0].module !== "env" || memories[0].name !== "memory")) {
    throw new Error("threaded artifact requires wasi_thread_start and imported shared memory");
  }
  if (!threaded && (memories[0].shared || names.has("wasi_thread_start"))) throw new Error("inconsistent serial artifact thread interface");
  const kind = `${threaded ? "threaded" : "serial"}-${command ? "command" : "reactor"}`;
  if (expected && kind !== expected) throw new Error(`expected ${expected}, got ${kind}`);
  return { module, identity: { kind, sha256: sha256 ?? crypto.createHash("sha256").update(bytes).digest("hex"), imports, exports, memory: memories[0] } };
}
// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); ... return status;
// Preserve the metadata-only API and CLI JSON; validation errors still propagate.
function inspectArtifact(bytes, expected) {
  return prepareArtifact(bytes, expected).identity;
}
// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); return status;
// Host compatibility only. A Wasm i32 reaches JavaScript as a signed number, so a
// pointer or a length of 2 GiB or more arrives negative. Node v24.21.0 and v26.8.2,
// src/node_wasi.cc:300-301,393-394:
//   bool CheckType<uint32_t>(Local<Value> value) { return value->IsUint32(); }
//   if (!CheckTypes<Args...>(args)) { args.GetReturnValue().Set(UVWASI_EINVAL); return; }
// so node:wasi answers EINVAL to such a call. Every i32 parameter of WASI preview1 is
// unsigned (pointer, size, descriptor, flags); pass each as its u32 value. The i64
// parameters are BigInt and pass unchanged.
function unsignedWasiImports(wasiImport) {
  const imports = {};
  for (const [name, fn] of Object.entries(wasiImport)) {
    imports[name] = (...args) => {
      for (let i = 0; i < args.length; i += 1) if (typeof args[i] === "number") args[i] >>>= 0;
      return fn(...args);
    };
  }
  return imports;
}
// NCBI reference: c++/src/algo/blast/blastinput/blast_args.cpp:3205-3222
// m_NumThreads = num_threads;
// Host memory only. The module declares the largest shared memory it can use (65,535
// pages, LOSAT/.cargo/config.toml), and a provided memory may have a smaller maximum
// (WebAssembly import matching). An engine reserves a shared memory's maximum as address
// space when it creates it. JavaScriptCore throws RangeError when it cannot. V8 without
// guard regions (Android, or an address-space limit) keeps a smaller reservation instead,
// Node v26.8.2 deps/v8/src/objects/backing-store.cc:454-466:
//   // Retry with smaller maximum pages at each retry.
//   size_t sizes[] = {maximum_pages - delta, maximum_pages - 2 * delta,
//                     maximum_pages - 3 * delta, initial_pages};
// and with the last of these the memory cannot grow at all. So create each candidate one
// page short of the initial size and grow that page; a maximum that fails either step is
// replaced by the next smaller one, down to the 1 GiB that the modules declared before.
const MEMORY_STEP_DOWN_PAGES = [49152, 32768, 24576, 16384];

function createSharedMemory(initial, maximum, stepDown, create = descriptor => new WebAssembly.Memory(descriptor)) {
  const smaller = stepDown ? MEMORY_STEP_DOWN_PAGES.filter(pages => pages < maximum && pages >= initial) : [];
  let failure;
  for (const pages of [maximum, ...smaller]) {
    try {
      const memory = create({ initial: initial - 1, maximum: pages, shared: true });
      memory.grow(1);
      return { memory, maximum: pages };
    } catch (error) {
      if (!(error instanceof RangeError)) throw error;
      failure = error;
    }
  }
  throw failure;
}

module.exports = { inspectArtifact, prepareArtifact, unsignedWasiImports, createSharedMemory };
if (require.main === module) {
  try { console.log(JSON.stringify(inspectArtifact(require("node:fs").readFileSync(process.argv[2]), process.argv[3]), null, 2)); }
  catch (error) { console.error(error); process.exitCode = 1; }
}
