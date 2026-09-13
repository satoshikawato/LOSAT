"use strict";
const crypto = require("node:crypto");

// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); return status;
// WASI ABI boundary: command uses _start; reactor uses _initialize once.
// Inspect actual imports/exports and binary memory limits, never the filename.
function inspectArtifact(bytes, expected) {
  const module = new WebAssembly.Module(bytes);
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
  return { kind, sha256: crypto.createHash("sha256").update(bytes).digest("hex"), imports, exports, memory: memories[0] };
}
module.exports = { inspectArtifact };
if (require.main === module) {
  try { console.log(JSON.stringify(inspectArtifact(require("node:fs").readFileSync(process.argv[2]), process.argv[3]), null, 2)); }
  catch (error) { console.error(error); process.exitCode = 1; }
}
