"use strict";
// Line WO: real searches on the threaded reactor, then the low bytes of its memory: the cell of the
// shared-memory guard (address 16), and how deep the stack below 1 MiB went.
// usage: node reactor_low_memory.js <repo root> <losat-threaded-reactor.wasm>
const fs = require("node:fs");
const path = require("node:path");
const [root, artifact] = process.argv.slice(2);
const { createThreadHost } = require(path.join(path.resolve(root), "LOSAT/tests/wasi_thread_host"));
const read = (name) => fs.readFileSync(path.join(root, "LOSAT/tests/fasta", name));

function runPair(host, values) {
  const api = host.instance.exports;
  const allocations = values.map((value) => {
    const bytes = Buffer.from(value), ptr = api.losat_web_alloc(bytes.length) >>> 0;
    new Uint8Array(host.memory.buffer, ptr, bytes.length).set(bytes);
    return [ptr, bytes.length];
  });
  try {
    const status = api.losat_web_run_pair(...allocations.flat());
    const text = (ptr, len) => Buffer.from(new Uint8Array(host.memory.buffer, ptr >>> 0, len >>> 0)).toString();
    return { status, rows: text(api.losat_web_result_ptr(), api.losat_web_result_len()).split("\n").length - 1, error: text(api.losat_web_error_ptr(), api.losat_web_error_len()).slice(0, 80) };
  }
  finally { for (const [ptr, length] of allocations) api.losat_web_dealloc(ptr, length); }
}

(async () => {
  const host = await createThreadHost(artifact, "threaded-reactor");
  const cases = [
    ["blastn", "AP027131.fasta", "NZ_CP006932.fasta", ["-task", "blastn"]],
    ["blastn", "AP027131.fasta", "AP027132.fasta", ["-task", "megablast"]],
    ["blastn", "AP027131.fasta", "AP027132.fasta", ["-task", "dc-megablast"]],
    ["blastp", "AP027131.faa", "AP027132.faa", []],
    ["tblastn", "AP027131.faa", "AP027132.fasta", []],
    ["blastx", "AP027131.fasta", "AP027132.faa", []],
  ];
  try {
    for (const [program, query, subject, args] of cases) {
      const { status, rows, error } = runPair(host, [program, read(query), read(subject), "6", [...args, "-num_threads", "4"].join("\0")]);
      await host.waitForWorkers();
      const low = new Uint8Array(host.memory.buffer, 0, 1 << 20);
      let first = 20;
      while (first < low.length && low[first] === 0) first += 1;
      console.log(program, args.join(" "), "status", status, "rows", rows, "cell", new Int32Array(host.memory.buffer)[4], "pages", host.memory.grow(0),
        "deepest_stack", (1 << 20) - first, "below16", Buffer.from(low.subarray(0, 16)).toString("hex"), error);
    }
  } finally { await host.close(); }
})();
