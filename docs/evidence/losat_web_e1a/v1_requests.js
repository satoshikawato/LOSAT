"use strict";
// Sends the same BLASTP requests to a web ABI v1 threaded reactor and prints one JSON
// line per request: the status, the result length and SHA-256, and the error text.
// Run it with the reactor built before a change and with the one built after it; the
// two outputs must be identical (the v1 behaviour, including which error is reported
// first when a request has two errors, must not change).
//
// Usage: node v1_requests.js REACTOR.wasm FASTA
const crypto = require("node:crypto");
const fs = require("node:fs");
const path = require("node:path");
const { createThreadHost } = require(path.join(__dirname, "../../../LOSAT/tests/wasi_thread_host"));

function runPair(host, program, query, subject, format, args) {
  const api = host.instance.exports;
  const values = [program, query, subject, format, args.join("\0")];
  const allocations = values.map((value) => {
    const bytes = Buffer.from(value);
    const ptr = api.losat_web_alloc(bytes.length);
    new Uint8Array(host.memory.buffer, ptr, bytes.length).set(bytes);
    return [ptr, bytes.length];
  });
  try {
    const status = api.losat_web_run_pair(...allocations.flat());
    const view = (ptr, len) => Buffer.from(new Uint8Array(host.memory.buffer, ptr, len));
    const result = view(api.losat_web_result_ptr(), api.losat_web_result_len());
    const error = view(api.losat_web_error_ptr(), api.losat_web_error_len()).toString();
    return {
      status,
      length: result.length,
      sha256: crypto.createHash("sha256").update(result).digest("hex"),
      error,
    };
  } finally {
    for (const [ptr, length] of allocations) api.losat_web_dealloc(ptr, length);
  }
}

// [name, outfmt, other options]. Requests with two errors show which one is reported.
const REQUESTS = [
  ["bad fields + bad matrix", "6 nosuch", ["-matrix", "FOO"]],
  ["custom outfmt 0 + unsupported matrix", "0 qseqid", ["-matrix", "BLOSUM45"]],
  ["unsupported outfmt + ungapped", "9", ["-ungapped"]],
  ["unsupported outfmt + zero threads", "9", ["-num_threads", "0"]],
  ["unsupported outfmt + bad evalue", "9", ["-evalue", "abc"]],
  ["bad fields only", "6 nosuch", []],
  ["custom outfmt 0 only", "0 qseqid", []],
  ["unsupported outfmt only", "9", []],
  ["unsupported matrix only", "6", ["-matrix", "BLOSUM45"]],
  ["outfmt 0", "0", []],
  ["outfmt 0, max_target_seqs 1", "0", ["-max_target_seqs", "1"]],
  ["outfmt 6", "6", []],
  ["outfmt 7, 2 threads", "7", ["-num_threads", "2"]],
  ["custom fields", "6 qseqid sseqid qseq sseq btop", []],
];

(async () => {
  const [wasm, fasta] = process.argv.slice(2);
  const host = await createThreadHost(wasm, "threaded-reactor");
  const sequences = fs.readFileSync(fasta, "utf8");
  for (const [name, outfmt, args] of REQUESTS) {
    const result = runPair(host, "blastp", sequences, sequences, outfmt, args);
    await host.waitForWorkers();
    console.log(JSON.stringify({ name, outfmt, args, ...result }));
  }
  process.exit(0);
})().catch((error) => {
  console.error(error);
  process.exit(1);
});
