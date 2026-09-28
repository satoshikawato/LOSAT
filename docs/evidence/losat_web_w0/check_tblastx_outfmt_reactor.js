"use strict";
// Evidence check for S01: the rebuilt serial v1 reactor rejects TBLASTX outfmt 0 and 7
// and still returns outfmt 6. Usage: node check_tblastx_outfmt_reactor.js <repo> <reactor.wasm>
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");

const [repo, artifact] = process.argv.slice(2);
const { WASI } = require("node:wasi");

function fasta() {
  // Same source sequence as fixtures() in LOSAT/tests/check_wasm_threading.py (nuc1).
  const text = fs.readFileSync(path.join(repo, "LOSAT/tests/fasta/small_test.fasta"), "utf8");
  const seq = text.split(/\r?\n/).filter((l) => !l.startsWith(">")).map((l) => l.trim()).join("");
  return `>seq0\n${seq.slice(0, 900)}\n`;
}

function runPair(host, program, query, subject, format, args) {
  const api = host.instance.exports;
  const values = [program, query, subject, format, args.join("\0")];
  const allocations = values.map((value) => {
    const bytes = Buffer.from(value), ptr = api.losat_web_alloc(bytes.length);
    new Uint8Array(host.memory.buffer, ptr, bytes.length).set(bytes);
    return [ptr, bytes.length];
  });
  try {
    const status = api.losat_web_run_pair(...allocations.flat());
    const result = Buffer.from(new Uint8Array(host.memory.buffer, api.losat_web_result_ptr(), api.losat_web_result_len()));
    const error = Buffer.from(new Uint8Array(host.memory.buffer, api.losat_web_error_ptr(), api.losat_web_error_len())).toString();
    return { status, result, error };
  } finally {
    for (const [ptr, length] of allocations) api.losat_web_dealloc(ptr, length);
  }
}

(async () => {
  // Same serial-reactor hosting as LOSAT/tests/check_wasi_api_limits.js.
  const module = await WebAssembly.compile(fs.readFileSync(artifact));
  const wasi = new WASI({ version: "preview1", env: {}, preopens: {}, returnOnExit: true });
  const instance = await WebAssembly.instantiate(module, wasi.getImportObject());
  wasi.initialize(instance);
  const host = { instance, memory: instance.exports.memory };
  const seq = fasta();
  const rows = [];
  for (const [format, extra] of [["6", []], ["0", []], ["7", []], ["", ["-outfmt", "0"]], ["6 qseqid", []]]) {
    const r = runPair(host, "tblastx", seq, seq, format, extra);
    rows.push({
      outfmt_argument: format, extra_args: extra, status: r.status,
      result_bytes: r.result.length, result_sha256: crypto.createHash("sha256").update(r.result).digest("hex"),
      error: r.error,
    });
  }
  console.log(JSON.stringify({ artifact: path.basename(artifact), rows }, null, 2));
  const ok = rows[0].status === 0 && rows[0].result_bytes > 0 && rows.slice(1).every((r) => r.status === -1 && /unsupported TBLASTX outfmt/.test(r.error));
  process.exit(ok ? 0 : 1);
})().catch((e) => { console.error(e); process.exit(2); });
