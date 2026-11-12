"use strict";
// V-ABI (plan §6.2): runs the two adapter reactors under Node through ABI v2
// (docs/web/abi_v2.md) and compares every output format with the native CLI of the same
// commit run for that format alone, by SHA-256. The standard error of the CLI is
// compared with stream 3, and the HSP records (stream 1) are checked against the
// outfmt 6 rows and outfmt 0 sections they point to.
//
// One subject is registered once per program and reused by consecutive runs with other
// queries and options (plan §4.6). The serial reactor runs with 1 thread, the threaded
// reactor with 1, 2 and 4 (the host appends -num_threads, docs/web/abi_v2.md §7).
//
// The searches come from web/adapter/tools/v_abi_cases.py; where a frozen SHA-256
// exists (Gate A, the Stage G matrix, NCBI's fixtures) the match is recorded as well, for a
// search with `"frozen_prefix"` after the frozen path prefix is put back (`toFrozen`). A
// search marked `"expect": "failure"` must fail in the CLI and in the ABI with the same
// error. A search with `"outfmt"` (a format besides 0, 6 and 7, or 6 or 7 with fields) runs
// with that -outfmt: its output, stream 4, is compared with the CLI's run of the same argv,
// whose standard error is then stream 3 (abi_v2.md §5, §7); SAM's (`"sam": true`) without the
// CL tag of its @PG line, where each writes its own command line. The output of outfmt 13 and 14
// (`"files": true`) is a ZIP archive, whose entries are compared with the files of the CLI's run
// with `-out results`: the side files `results_<n>.json|xml`, and the master file as
// `results.json|xml`.
//
// A BLASTN search with `-window_masker_db` names a counts file: the host registers that file under
// the argv's value (`file_begin`, `file_chunk`, `file_end`; abi_v2.md §4) before the run, and keeps
// it for the next searches that name it, as the application keeps one (abi_v2.md §6); a value that
// is not a regular file is not registered, which the run reports as the CLI reports a file it
// cannot open.
//
// Usage: node v_abi.js --native losat --serial losat-web-serial.wasm
//          --threads losat-web-threads.wasm --cases CASES.json [--engines serial,threads]
//          [--out DIR]
// Several processes may run parts of the case list (for example one engine or one
// program each); each writes its own results.
const assert = require("node:assert/strict");
const { execFileSync } = require("node:child_process");
const crypto = require("node:crypto");
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const { WASI } = require("node:wasi");

const ROOT = path.resolve(__dirname, "../../..");
const { createThreadHost } = require(path.join(ROOT, "LOSAT/tests/wasi_thread_host"));
const { unsignedWasiImports } = require(path.join(ROOT, "LOSAT/tests/wasi_artifact"));

const FORMATS = { blastp: [0, 6, 7], tblastn: [0, 6, 7], blastn: [0, 6, 7], tblastx: [0, 6, 7], blastx: [0, 6, 7] };
// The stream of the output of the argv's -outfmt.
const CUSTOM_STREAM = 4;
const WITH_HITS = new Set(["blastp", "tblastn", "blastn", "tblastx", "blastx"]);

function sha256(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

// A SAM document without the CL tag of its @PG line (LOSAT/tests/db_search_fixtures.py
// `strip_sam_cl`: from the `\tCL:` after `\n@PG\tID:0\tVN:2.17.0+` up to the next `\tPN:blastn\n`).
function withoutCl(document) {
  const head = Buffer.from("\n@PG\tID:0\tVN:2.17.0+\tCL:");
  const at = document.indexOf(head);
  if (at < 0) return document;
  const start = at + head.length - "\tCL:".length;
  const end = document.indexOf(Buffer.from("\tPN:blastn\n"), start);
  if (end < 0) return document;
  return Buffer.concat([document.subarray(0, start), document.subarray(end)]);
}

// The bytes of a run as the frozen run wrote them, before their frozen hash is compared: every
// `prefix.run` that does not follow a character of `prefix.not_after` is `prefix.frozen` (the
// search's `"frozen_prefix"`, LOSAT/tests/frozen_paths.py `to_frozen`). Without one, the bytes.
function toFrozen(bytes, prefix) {
  if (!prefix) return bytes;
  const [run, frozen, notAfter] = [Buffer.from(prefix.run), Buffer.from(prefix.frozen), Buffer.from(prefix.not_after)];
  const parts = [];
  let from = 0, at = bytes.indexOf(run);
  while (at >= 0) {
    if (at === 0 || !notAfter.includes(bytes[at - 1])) {
      parts.push(bytes.subarray(from, at), frozen);
      from = at + run.length;
      at = bytes.indexOf(run, from);
    } else {
      at = bytes.indexOf(run, at + 1);
    }
  }
  parts.push(bytes.subarray(from));
  return Buffer.concat(parts);
}

// The CRC-32 of ZIP (ISO-HDLC).
function crc32(bytes) {
  let crc = ~0;
  for (const byte of bytes) {
    crc ^= byte;
    for (let k = 0; k < 8; k++) crc = crc & 1 ? 0xedb88320 ^ (crc >>> 1) : crc >>> 1;
  }
  return ~crc >>> 0;
}

// The entries of a ZIP archive of stored files (the adapter's archive of outfmt 13 and 14), as
// [name, bytes] in the order of its central directory, read as a ZIP reader reads them: the end
// record (or the ZIP64 end record), each central header, its local header with the same CRC-32 and
// sizes, and the bytes, whose CRC-32 is checked.
function storedZipEntries(archive) {
  const end = archive.length - 22;
  assert.equal(archive.readUInt32LE(end), 0x06054b50, "the end of central directory record");
  let count = archive.readUInt16LE(end + 10);
  let at = archive.readUInt32LE(end + 16);
  if (count === 0xffff) {
    const record = Number(archive.readBigUInt64LE(end - 20 + 8));
    assert.equal(archive.readUInt32LE(record), 0x06064b50, "the ZIP64 end of central directory record");
    count = Number(archive.readBigUInt64LE(record + 32));
    at = Number(archive.readBigUInt64LE(record + 48));
  }
  const entries = [];
  for (let i = 0; i < count; i++) {
    assert.equal(archive.readUInt32LE(at), 0x02014b50, "a central header");
    const [crc, size, nameLength, extraLength] = [archive.readUInt32LE(at + 16), archive.readUInt32LE(at + 24),
                                                  archive.readUInt16LE(at + 28), archive.readUInt16LE(at + 30)];
    assert.equal(archive.readUInt16LE(at + 10), 0, "a stored entry");
    assert.equal(archive.readUInt32LE(at + 20), size, "stored: the compressed size is the size");
    const name = archive.subarray(at + 46, at + 46 + nameLength).toString("utf8");
    let local = archive.readUInt32LE(at + 42);
    if (local === 0xffffffff) local = Number(archive.readBigUInt64LE(at + 46 + nameLength + 4));
    assert.equal(archive.readUInt32LE(local), 0x04034b50, `${name}: a local header`);
    assert.equal(archive.readUInt32LE(local + 14), crc, `${name}: the local CRC-32`);
    assert.equal(archive.readUInt32LE(local + 22), size, `${name}: the local size`);
    const data = local + 30 + archive.readUInt16LE(local + 26) + archive.readUInt16LE(local + 28);
    const bytes = archive.subarray(data, data + size);
    assert.equal(crc32(bytes), crc, `${name}: the CRC-32 of its bytes`);
    entries.push([name, bytes]);
    at += 46 + nameLength + extraLength;
  }
  return entries;
}

// The files of the CLI's run of outfmt 13 or 14 with `-out <dir>/results`, as the adapter's archive
// names them: the side files in order, then the master file as results.json or results.xml.
function nativeFiles(directory, outfmt) {
  const extension = outfmt.trim().startsWith("13") ? "json" : "xml";
  const files = [];
  for (let n = 1; fs.existsSync(path.join(directory, `results_${n}.${extension}`)); n++) {
    files.push([`results_${n}.${extension}`, fs.readFileSync(path.join(directory, `results_${n}.${extension}`))]);
  }
  files.push([`results.${extension}`, fs.readFileSync(path.join(directory, "results"))]);
  return files;
}

// The -outfmt numbers whose warnings are those of outfmt 6 (the tabular writer's): the structured
// formats do not write the queries' warnings (`ReportFormat::writes_query_warnings`).
const TABULAR = new Set([6, 7, 10, 20]);

function parseArgs(argv) {
  const options = {};
  for (let i = 0; i < argv.length; i += 2) options[argv[i].replace(/^--/, "")] = argv[i + 1];
  for (const required of ["native", "serial", "threads", "cases"]) {
    if (!options[required]) throw new Error(`--${required} is required`);
  }
  // Searches run from their own working directory (the repository or LOSAT/), so a relative
  // CLI path is resolved against the directory v_abi.js was started from.
  options.native = path.resolve(options.native);
  return options;
}

// ---------------------------------------------------------------------------
// The ABI, over one reactor instance.

class Reactor {
  constructor(exports, memory) {
    this.exports = exports;
    this.memory = memory;
    this.streams = new Map();
  }

  // Pointers and lengths are unsigned; JavaScript receives an i32 as a signed number
  // (docs/web/abi_v2.md §6), so a pointer of 2 GiB or more arrives negative.
  emit(stream, ptr, len) {
    const chunk = Buffer.from(new Uint8Array(this.memory().buffer, ptr >>> 0, len >>> 0));
    if (!this.streams.has(stream)) this.streams.set(stream, []);
    this.streams.get(stream).push(chunk);
  }

  take() {
    const result = new Map();
    for (const [stream, chunks] of this.streams) result.set(stream, Buffer.concat(chunks));
    this.streams.clear();
    return result;
  }

  error() {
    const e = this.exports;
    return Buffer.from(new Uint8Array(this.memory().buffer, e.losat_web2_last_error_ptr() >>> 0, e.losat_web2_last_error_len() >>> 0)).toString();
  }

  // Calls an export; Buffer and string arguments become (ptr, len) pairs.
  call(name, ...args) {
    const e = this.exports, allocations = [], flat = [];
    for (const arg of args) {
      if (typeof arg === "number") { flat.push(arg); continue; }
      const bytes = Buffer.isBuffer(arg) ? arg : Buffer.from(arg);
      const ptr = e.losat_web2_alloc(bytes.length) >>> 0;
      new Uint8Array(this.memory().buffer, ptr, bytes.length).set(bytes);
      allocations.push([ptr, bytes.length]);
      flat.push(ptr, bytes.length);
    }
    try {
      const status = e[name](...flat);
      return { status, error: status < 0 ? this.error() : "", streams: this.take() };
    } finally {
      for (const [ptr, len] of allocations) e.losat_web2_dealloc(ptr, len);
    }
  }
}

// Registers the file at `file` under `name` in chunks of `chunkSize` bytes through one buffer of
// the reactor (abi_v2.md §4 `file_begin`); returns the handle and the *file* JSON.
function registerFile(reactor, name, file, chunkSize = 1 << 20) {
  const e = reactor.exports, bytes = fs.readFileSync(file);
  const begun = reactor.call("losat_web2_file_begin", 1, name, bytes.length);
  assert.ok(begun.status > 0, `file_begin ${name}: ${begun.error}`);
  const size = Math.max(1, Math.min(chunkSize, bytes.length));
  const ptr = e.losat_web2_alloc(size) >>> 0;
  try {
    for (let at = 0; at < bytes.length; at += size) {
      const part = bytes.subarray(at, at + size);
      new Uint8Array(reactor.memory().buffer, ptr, part.length).set(part);
      const status = e.losat_web2_file_chunk(begun.status, ptr, part.length);
      assert.equal(status, 0, `file_chunk ${name}: ${status < 0 ? reactor.error() : ""}`);
    }
  } finally {
    e.losat_web2_dealloc(ptr, size);
  }
  const ended = reactor.call("losat_web2_file_end", begun.status);
  assert.equal(ended.status, 0, `file_end ${name}: ${ended.error}`);
  const response = JSON.parse(ended.streams.get(2).toString());
  assert.equal(response.handle, begun.status);
  assert.equal(response.name, name);
  assert.equal(response.kind, "windowmasker_counts");
  assert.equal(response.size, bytes.length);
  assert.equal(response.sha256, sha256(bytes), `${name}: the SHA-256 of the bytes given`);
  return { handle: begun.status, response };
}

// The counts files of the searches: the argv's `-window_masker_db` value -> its handle and path.
// A search that names a registered value with another file releases it first.
function registerCounts(reactor, counts, search) {
  const at = search.argv.indexOf("-window_masker_db");
  if (at < 0) return null;
  const name = search.argv[at + 1];
  const file = path.resolve(search.cwd, name);
  const known = counts.get(name);
  if (known && known.file === file) return { name, registered: false };
  if (known) {
    assert.equal(reactor.call("losat_web2_file_release", known.handle).status, 0);
    counts.delete(name);
  }
  if (name === "" || !fs.existsSync(file) || !fs.statSync(file).isFile()) return { name, registered: false };
  const { handle } = registerFile(reactor, name, file);
  counts.set(name, { handle, file });
  return { name, registered: true };
}

async function openSerial(wasmPath) {
  const wasi = new WASI({ version: "preview1", returnOnExit: true, args: [], env: {} });
  let reactor;
  const imports = {
    wasi_snapshot_preview1: unsignedWasiImports(wasi.wasiImport),
    losat_host: { emit: (stream, ptr, len) => reactor.emit(stream, ptr, len) },
  };
  const { instance } = await WebAssembly.instantiate(fs.readFileSync(wasmPath), imports);
  wasi.initialize(instance);
  reactor = new Reactor(instance.exports, () => instance.exports.memory);
  return { reactor, close: async () => {} };
}

async function openThreads(wasmPath) {
  let reactor;
  const host = await createThreadHost(wasmPath, "threaded-reactor", [], {
    imports: { losat_host: { emit: (stream, ptr, len) => reactor.emit(stream, ptr, len) } },
  });
  reactor = new Reactor(host.instance.exports, () => host.memory);
  return { reactor, close: () => host.close(), waitForWorkers: () => host.waitForWorkers() };
}

// ---------------------------------------------------------------------------

function inputPath(search, option) {
  return path.resolve(search.cwd, search.argv[search.argv.indexOf(option) + 1]);
}

// The CLI's outputs of a search, one run per format with -out. A search with
// `"expect": "failure"` fails in the CLI (exit 3 or 1), whose files keep what was written
// before the failure (the outfmt 0 prolog before a first query batch fails).
function nativeOutputs(native, search, scratch) {
  const outputs = new Map();
  let stderr = null, status = null;
  const stderrs = new Map();
  for (const format of FORMATS[search.program]) {
    const out = path.join(scratch, `native.${format}.out`);
    const argv = [...search.argv, "-outfmt", String(format), "-num_threads", "1", "-out", out];
    const result = require("node:child_process").spawnSync(native, argv, { cwd: search.cwd, maxBuffer: 1 << 30 });
    const failed = result.status !== 0;
    if (failed !== (search.expect === "failure")) {
      throw new Error(`native ${argv.join(" ")} ${failed ? "failed" : "succeeded"} (exit ${result.status}): ${result.error ?? result.stderr}`);
    }
    outputs.set(format, fs.existsSync(out) ? fs.readFileSync(out) : Buffer.alloc(0));
    fs.rmSync(out, { force: true });
    stderrs.set(format, result.stderr);
    if (stderr === null) [stderr, status] = [result.stderr, result.status];
  }
  // The diagnostics are the CLI's standard error for the tabular formats (abi_v2.md §5):
  // outfmt 6 and 7 write the same lines, and outfmt 0 the same without the formatting
  // handler's lines, which -num_descriptions and -num_alignments can make differ.
  const handler = /^Warning: \[\w+\] (The parameter -num_descriptions is ignored|Examining 5 or more matches)/;
  const lines = (format) => stderrs.get(format).toString().split("\n");
  const withoutHandler = (format) => lines(format).filter((line) => !handler.test(line));
  if (stderrs.has(6)) {
    stderr = stderrs.get(6);
    if (stderrs.has(7)) assert.equal(sha256(stderrs.get(7)), sha256(stderr), "outfmt 7 writes the warnings of outfmt 6");
    for (const format of FORMATS[search.program].filter((f) => f !== 6)) {
      assert.deepEqual(withoutHandler(format), withoutHandler(6), `outfmt ${format}: the same warnings but the handler's`);
      for (const line of lines(format)) assert.ok(lines(6).includes(line), `outfmt ${format}: ${line} is in the warnings of outfmt 6`);
    }
  }
  if (search.outfmt !== undefined) {
    const directory = fs.mkdtempSync(path.join(scratch, "native.custom."));
    const out = path.join(directory, search.files ? "results" : "native.custom.out");
    const argv = [...search.argv, "-outfmt", search.outfmt, "-num_threads", "1", "-out", out];
    const result = require("node:child_process").spawnSync(native, argv, { cwd: search.cwd, maxBuffer: 1 << 30 });
    if (result.status !== 0) throw new Error(`native ${argv.join(" ")} failed (exit ${result.status}): ${result.error ?? result.stderr}`);
    outputs.set(CUSTOM_STREAM, search.files ? nativeFiles(directory, search.outfmt) : fs.readFileSync(out));
    fs.rmSync(directory, { recursive: true, force: true });
    // The warnings of outfmt 6, and those of the fields' writer (NCBI's taxdb warning).
    if (TABULAR.has(Number(search.outfmt.trim().split(" ")[0]))) {
      for (const line of lines(6)) assert.ok(result.stderr.toString().split("\n").includes(line), `-outfmt ${search.outfmt}: ${line} is in its warnings`);
    }
    stderr = result.stderr;
  }
  return { outputs, stderr, status };
}

// A run that fails as the CLI fails: the CLI's error is the last message of its standard
// error (the warnings before it go to stream 3, which a failed run does not send), and
// whatever the run sent of an output before it failed is the start of the CLI's file. A search with
// `"prolog": false` fails before any output (an error of the options or of their files), so the
// CLI's files are empty; the others fail after the CLI wrote the outfmt 0 prolog.
function checkFailure(run, expected, search, label) {
  const program = search.program;
  assert.equal(run.status, -1, `${label}: the run fails as the CLI (exit ${expected.status})`);
  assert.ok(run.error.length > 0 && expected.stderr.toString().endsWith(run.error), `${label}: the CLI's error (${run.error})`);
  for (const format of FORMATS[program]) {
    const sent = run.streams.get(format) || Buffer.alloc(0);
    assert.ok(expected.outputs.get(format).subarray(0, sent.length).equals(sent), `${label}: outfmt ${format} before the failure`);
  }
  if (search.prolog === false) {
    assert.equal(expected.outputs.get(0).length, 0, `${label}: the CLI writes no outfmt 0 before the failure`);
  } else {
    assert.ok(expected.outputs.get(0).length > 0, `${label}: the CLI's outfmt 0 file has the prolog`);
  }
}

function checkHits(streams, program, label, queryIds) {
  const records = (streams.get(1) || Buffer.alloc(0)).toString().split("\n").filter(Boolean).map((line) => JSON.parse(line));
  // A stream without bytes (a search without hits) is never emitted.
  const out6 = streams.get(6) || Buffer.alloc(0), out0 = streams.get(0) || Buffer.alloc(0);
  const rows = out6.toString().split("\n").filter(Boolean);
  assert.equal(records.length, rows.length, `${label}: one record per outfmt 6 row`);
  const ranks = new Map();
  records.forEach((record, index) => {
    assert.equal(record.index, index, `${label}: index`);
    const rank = ranks.get(record.q_idx) || 0;
    assert.equal(record.rank, rank, `${label}: rank`);
    ranks.set(record.q_idx, rank + 1);
    const row = out6.subarray(record.out6[0], record.out6[1]).toString();
    assert.equal(row.trimEnd(), rows[index], `${label}: out6 of HSP ${index}`);
    const fields = row.trimEnd().split("\t");
    assert.deepEqual(fields.slice(6, 10).map(Number), [record.q_start, record.q_end, record.s_start, record.s_end], `${label}: coordinates of HSP ${index}`);
    // The records come in query order, each with the index of its query among the registered
    // records (a block search numbers a block's queries from the block's first query); where the
    // query's id (`register`'s response) is an accession-like word, it is the row's qseqid.
    assert.ok(record.q_idx < queryIds.length, `${label}: q_idx of HSP ${index} within the registered queries`);
    if (index > 0) assert.ok(record.q_idx >= records[index - 1].q_idx, `${label}: q_idx order at HSP ${index}`);
    const id = queryIds[record.q_idx];
    if (/^[A-Za-z_]+[0-9]+(\.[0-9]+)?$/.test(id)) assert.equal(fields[0], id, `${label}: q_idx of HSP ${index}`);
    // A subject beyond the alignments that outfmt 0 shows (for example BLASTN's 250) has
    // no section and no heading.
    if (record.out0 === null) {
      assert.equal(record.out0_subject, null, `${label}: out0_subject of HSP ${index} without a section`);
      return;
    }
    const section = out0.subarray(record.out0[0], record.out0[1]).toString();
    assert.ok(section.startsWith(" Score ="), `${label}: out0 of HSP ${index}`);
    assert.ok(section.includes(`bits (${record.raw_score}),`), `${label}: raw score in out0 of HSP ${index}`);
    const heading = out0.subarray(record.out0_subject[0], record.out0_subject[1]).toString();
    assert.ok(heading.startsWith("> ") && record.out0_subject[1] <= record.out0[0], `${label}: out0_subject of HSP ${index}`);
  });
  return records.length;
}

// describe, validate, register and scan: the responses and errors the application sees.
function checkSurface(reactor, native) {
  for (const program of Object.keys(FORMATS)) {
    const described = reactor.call("losat_web2_describe", program);
    assert.equal(described.status, 0, described.error);
    const description = JSON.parse(described.streams.get(2).toString());
    assert.deepEqual(description.formats, FORMATS[program], `${program}: formats`);
    const flags = description.parameters.map((parameter) => parameter.flag);
    assert.ok(flags.includes("-evalue") && flags.includes("-query"), `${program}: parameters`);
    assert.ok(!flags.includes("-outfmt") && !flags.includes("-out") && !flags.includes("-num_threads"), `${program}: adapter-owned options`);
    assert.ok(flags.includes("-num_descriptions") && flags.includes("-num_alignments"), `${program}: hit list options`);
    if (program === "tblastn") assert.ok(description.subject_gencodes.includes(32), "tblastn: genetic code 32");
    if (program === "blastx") assert.ok(description.query_gencodes.includes(11) && !("subject_gencodes" in description), "blastx: the query's genetic codes");
    if (program === "blastn") {
      assert.ok(flags.includes("-template_type") && flags.includes("-template_length"), "blastn: dc-megablast templates");
      assert.ok(flags.includes("-window_masker_db"), "blastn: -window_masker_db");
      const task = description.parameters.find((parameter) => parameter.flag === "-task");
      assert.ok(task.help.includes("dc-megablast") && task.help.includes("blastn-short"), "blastn: tasks in -task's help");
    }
  }
  // validate: accepted argv, and CLI errors with the CLI's words.
  assert.equal(reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s"].join("\0")).status, 0);
  for (const argv of [["blastp", "-query", "q", "-subject", "s", "-evalue", "abc"],
                      ["tblastn", "-query", "q", "-subject", "s", "-db_gencode", "7"],
                      ["blastn", "-query", "q", "-subject", "s", "-nosuchoption", "1"],
                      ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast", "-template_type", "coding"],
                      ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast", "-template_type", "coding",
                       "-template_length", "17"],
                      ["blastn", "-query", "q", "-subject", "s", "-task", "rmblastn"]]) {
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, argv.join(" "));
    const cli = require("node:child_process").spawnSync(native, argv);
    assert.equal(validated.error, cli.stderr.toString(), `${argv.join(" ")}: the CLI's message`);
  }
  // BLASTX's ranges (line SX, DW-39): validate reads their syntax, the subject range first.
  const blastx = ["blastx", "-query", "q", "-subject", "s"];
  assert.equal(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "1-9", "-subject_loc", "2-5"].join("\0")).status, 0);
  assert.match(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "9-1", "-subject_loc", "5-2"].join("\0")).error,
               /^BLAST engine error: Invalid specification of subject location/);
  assert.match(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "9-1"].join("\0")).error,
               /^BLAST engine error: Invalid specification of query location/);
  // The host validates the argv that it runs, with its -num_threads, so a -num_threads in
  // the user's words is a repeated option (abi_v2.md §7). The adapter parses the argv as the
  // CLI does, so the whole message (its usage line too) and an unknown program's are the CLI's.
  for (const argv of [["blastn", "-query", "q", "-subject", "s", "-num_threads", "2", "-num_threads", "1"],
                      ["nosuch", "-query", "q", "-subject", "s"]]) {
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, argv.join(" "));
    const cli = require("node:child_process").spawnSync(native, argv);
    assert.equal(validated.error, cli.stderr.toString(), `${argv.join(" ")}: the CLI's error`);
  }
  // -outfmt (WA2-OP-03, W11): the formats that the CLI writes and 6 and 7 with fields are valid (13 and 14 without
  // the CLI's -out), and the CLI's errors have the CLI's words.
  for (const outfmt of ["6", "0", "6 qseqid sseqid pident", "7 std qcovs", "5", "9", "10 qseqid", "11", "13", "14", "16", "17 SQ", "20"]) {
    assert.equal(reactor.call("losat_web2_validate", ["blastn", "-query", "q", "-subject", "s", "-outfmt", outfmt].join("\0")).status, 0, outfmt);
  }
  for (const [program, outfmt] of [["blastp", "13"], ["blastx", "11"], ["blastx", "14"]]) {
    assert.equal(reactor.call("losat_web2_validate", [program, "-query", "q", "-subject", "s", "-outfmt", outfmt].join("\0")).status, 0, `${program} ${outfmt}`);
  }
  // The CLI reads the files before its formatting handler checks the delimiter: real files.
  const files = ["-query", path.join(ROOT, "LOSAT/tests/fasta/outfmt0/mask_query.fasta"),
                 "-subject", path.join(ROOT, "LOSAT/tests/fasta/outfmt0/mask_subject.fasta")];
  for (const outfmt of ["abc", "99", "6 delim=; qseqid staxids"]) {
    const argv = ["blastn", ...files, "-outfmt", outfmt];
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, outfmt);
    const cli = require("node:child_process").spawnSync(native, argv);
    assert.equal(validated.error, cli.stderr.toString(), `-outfmt ${outfmt}: the CLI's message`);
  }
  assert.match(reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-out", "x"].join("\0")).error, /not accepted/);
  // The tasks of Session SD and NCBI's option check of a template (blast_options.c:1247-1261).
  for (const task of ["dc-megablast", "blastn-short"]) {
    assert.equal(reactor.call("losat_web2_validate", ["blastn", "-query", "q", "-subject", "s", "-task", task].join("\0")).status, 0, task);
  }
  assert.match(reactor.call("losat_web2_validate", ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast",
                                                    "-word_size", "13"].join("\0")).error,
               /^BLAST query\/options error: Invalid discontiguous template parameters: word size must be either 11 or 12/);
  // BLASTN scores without a Karlin-Altschul table row (WA-FM-01). The CLI meets the tables only
  // when some query context has an ungapped Karlin-Altschul block, which depends on the query's
  // letters (a query of N only has none: a warning and rc 0). `validate` has no query, so it
  // accepts these scores and `run` rejects them; whatever `validate` rejects, the CLI rejects
  // with the same message.
  {
    const dir = fs.mkdtempSync(path.join(os.tmpdir(), "losat-v-abi-scores-"));
    const acgt = "ACGTACGTAGCTAGCTAGGATCGATCGATCGTAGCTAGCTAGCTAGGCTAGCATCGATCGA";
    fs.writeFileSync(path.join(dir, "acgt.fa"), `>mq1 query\n${acgt}\n`);
    fs.writeFileSync(path.join(dir, "n.fa"), `>mq1 query\n${"N".repeat(60)}\n`);
    fs.writeFileSync(path.join(dir, "b.fa"), `>mq1 query\n${acgt}B\n`);
    fs.writeFileSync(path.join(dir, "s.fa"), `>ms1\n${acgt}\n`);
    const env = Object.fromEntries(Object.entries(process.env).filter(([name]) => !name.startsWith("LOSAT_")));
    const cases = [
      // [query, extra options, the CLI's exit status]
      ["acgt.fa", ["-reward", "7"], 0], ["acgt.fa", ["-reward", "3", "-penalty", "-1"], 0],
      ["acgt.fa", ["-task", "blastn", "-reward", "5", "-penalty", "-1"], 0],
      ["acgt.fa", ["-task", "blastn", "-reward", "7"], 3], ["n.fa", ["-task", "blastn", "-reward", "7"], 0],
      ["b.fa", ["-reward", "3", "-penalty", "-1"], null], ["acgt.fa", ["-reward", "2", "-penalty", "-1"], 3],
    ];
    for (const [query, extra, status] of cases) {
      const label = `${query} ${extra.join(" ")}`;
      const argv = ["blastn", "-query", path.join(dir, query), "-subject", path.join(dir, "s.fa"), ...extra];
      const cli = require("node:child_process").spawnSync(native, [...argv, "-outfmt", "6"], { env });
      if (status === null) assert.notEqual(cli.status, 0, label);
      else assert.equal(cli.status, status, label);
      if (cli.status === 0) assert.match(cli.stderr.toString(), /^Warning: \[blastn\] Query_1 mq1 query: Could not calculate ungapped Karlin-Altschul/, label);
      const validated = reactor.call("losat_web2_validate", argv.join("\0"));
      assert.equal(validated.status, 0, `${label}: validate leaves the tables to run`);
    }
    // What validate rejects (the options' checks), the CLI rejects with the same message.
    const options = ["blastn", "-query", path.join(dir, "acgt.fa"), "-subject", path.join(dir, "s.fa"), "-reward", "7", "-evalue", "0"];
    const validated = reactor.call("losat_web2_validate", options.join("\0"));
    assert.equal(validated.status, -1);
    const cli = require("node:child_process").spawnSync(native, options, { env });
    assert.equal(validated.error, cli.stderr.toString(), "the CLI's message");
    fs.rmSync(dir, { recursive: true });
  }
  // register and the scan of the program's kind (2: the TBLASTN query) agree on a
  // multi-record input, in chunks of any size.
  const input = fs.readFileSync(path.join(ROOT, "docs/evidence/tlosan_stage_c/multi_query_20260924/query.faa"));
  const registered = reactor.call("losat_web2_register", "tblastn", 0, input);
  assert.ok(registered.status > 0, registered.error);
  const records = JSON.parse(registered.streams.get(2).toString()).records;
  const scanner = reactor.exports.losat_web2_scan_begin(2);
  assert.ok(scanner > 0);
  for (let at = 0; at < input.length; at += 37) {
    assert.equal(reactor.call("losat_web2_scan_chunk", scanner, input.subarray(at, at + 37)).status, 0);
  }
  const scanned = reactor.call("losat_web2_scan_end", scanner);
  assert.equal(scanned.status, 0, scanned.error);
  const index = JSON.parse(scanned.streams.get(2).toString()).records;
  assert.deepEqual(index.map((r) => [r.id, r.length]), records.map((r) => [r.id, r.length]), "scan and register agree");
  assert.equal(reactor.call("losat_web2_release", registered.status).status, 0);
  assert.equal(reactor.call("losat_web2_release", registered.status).status, -1, "a released handle is gone");
  // NCBI's protein reader reads a first record without a defline (its ID is empty), and
  // LOSAT Web rejects a gap line, which its index cannot locate (abi_v2.md §4, §9).
  const headless = reactor.call("losat_web2_register", "blastp", 0, Buffer.from("ACGT\n>x\nAC\n"));
  assert.ok(headless.status > 0, headless.error);
  assert.deepEqual(JSON.parse(headless.streams.get(2).toString()).records.map((r) => [r.id, r.length]),
                   [["", 4], ["x", 2]]);
  assert.equal(reactor.call("losat_web2_release", headless.status).status, 0);
  const gap = reactor.call("losat_web2_register", "blastn", 1, Buffer.from(">s\nACGT\n>?10\nACGT\n"));
  assert.equal(gap.status, -1);
  assert.match(gap.error, /^line 3 is a gap line .*not supported in the browser$/);
  checkFiles(reactor);
}

// The named input files (abi_v2.md §4 `file_begin`): the same *file* JSON whatever the chunks, the
// sizes kept, one registration per name, and the errors of the file's contents in the JSON.
function checkFiles(reactor) {
  const counts = path.join(ROOT, "LOSAT/tests/winmask/counts");
  const whole = registerFile(reactor, "u9.obinary", path.join(counts, "mg1655_u9.obinary"));
  assert.equal(whole.response.stat_type, "obinary");
  assert.equal(whole.response.unit_size, 9);
  assert.equal(whole.response.error, null);
  const twice = reactor.call("losat_web2_file_begin", 1, "u9.obinary", 1);
  assert.equal(twice.status, -1, "a registered name is not registered again");
  assert.match(twice.error, /^an input file named u9\.obinary is registered/);
  assert.equal(reactor.call("losat_web2_file_release", whole.handle).status, 0);
  assert.equal(reactor.call("losat_web2_file_release", whole.handle).status, -1, "a released file is gone");
  const small = registerFile(reactor, "u9.obinary", path.join(counts, "mg1655_u9.obinary"), 4093);
  const { handle, ...same } = small.response;
  const { handle: first, ...expected } = whole.response;
  assert.deepEqual(same, expected, "chunks of 4093 bytes: the same file");
  assert.equal(reactor.call("losat_web2_file_release", small.handle).status, 0);
  const binary = registerFile(reactor, "u9.binary", path.join(counts, "dm2L_u9.binary"));
  assert.equal(binary.response.stat_type, "binary");
  assert.equal(reactor.call("losat_web2_file_release", binary.handle).status, 0);
  for (const [kind, name, error] of [[2, "x", /^unknown input file kind 2$/], [1, "", /^an input file needs a name$/]]) {
    const refused = reactor.call("losat_web2_file_begin", kind, name, 1);
    assert.equal(refused.status, -1);
    assert.match(refused.error, error);
  }
  const over = reactor.call("losat_web2_file_begin", 1, "over", 2);
  assert.equal(reactor.call("losat_web2_file_chunk", over.status, Buffer.from("abc")).status, -1, "past the size");
  assert.equal(reactor.call("losat_web2_file_end", over.status).status, -1, "a file past its size is dropped");
  const short = reactor.call("losat_web2_file_begin", 1, "short", 2);
  assert.equal(reactor.call("losat_web2_file_chunk", short.status, Buffer.from("a")).status, 0);
  const ended = reactor.call("losat_web2_file_end", short.status);
  assert.equal(ended.status, -1);
  assert.equal(ended.error, "the input file short ended after 1 of its 2 bytes");
  const cut = fs.readFileSync(path.join(counts, "mg1655_u9.obinary"));
  const begun = reactor.call("losat_web2_file_begin", 1, "cut", cut.length >> 1);
  assert.equal(reactor.call("losat_web2_file_chunk", begun.status, cut.subarray(0, cut.length >> 1)).status, 0);
  const content = reactor.call("losat_web2_file_end", begun.status);
  assert.equal(content.status, 0, "an error of the contents comes with the run");
  assert.match(JSON.parse(content.streams.get(2).toString()).error, /^not enough data/);
  assert.equal(reactor.call("losat_web2_file_release", begun.status).status, 0);
}

async function main() {
  const options = parseArgs(process.argv.slice(2));
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), "losat-v-abi-"));
  const results = [];
  const selected = (options.engines || "serial,threads").split(",");
  const engines = [
    { name: "serial", open: () => openSerial(options.serial), threads: [1] },
    { name: "threads", open: () => openThreads(options.threads), threads: [1, 2, 4] },
  ].filter((engine) => selected.includes(engine.name));
  const all = JSON.parse(fs.readFileSync(options.cases, "utf8"));
  for (const engine of engines) {
    const { reactor, close, waitForWorkers } = await engine.open();
    assert.equal(reactor.exports.losat_web2_abi_version(), 2);
    checkSurface(reactor, options.native);
    console.log(`ok ${engine.name} describe, validate, register, scan`);
    const subjects = new Map(); // "program subject" -> handle, kept across runs
    const counts = new Map(); // the -window_masker_db value -> { handle, file }, kept across runs
    let countsRegistered = 0, countsKept = 0;
    for (const search of all) {
      const program = search.program;
      const query = inputPath(search, "-query"), subject = inputPath(search, "-subject");
      const expected = nativeOutputs(options.native, search, scratch);
      const subjectKey = `${program} ${subject}`;
      if (!subjects.has(subjectKey)) {
        const registered = reactor.call("losat_web2_register", program, 1, fs.readFileSync(subject));
        assert.ok(registered.status > 0, `register subject: ${registered.error}`);
        subjects.set(subjectKey, registered.status);
      }
      const q = reactor.call("losat_web2_register", program, 0, fs.readFileSync(query));
      assert.ok(q.status > 0, `register query: ${q.error}`);
      const queryIds = JSON.parse(q.streams.get(2).toString()).records.map((record) => record.id);
      const named = registerCounts(reactor, counts, search);
      if (named) named.registered ? countsRegistered++ : countsKept++;
      for (const threads of engine.threads) {
        const custom = search.outfmt === undefined ? [] : ["-outfmt", search.outfmt];
        const argv = [...search.argv, ...custom, "-num_threads", String(threads)].join("\0");
        const run = reactor.call("losat_web2_run", argv, q.status, subjects.get(subjectKey));
        if (waitForWorkers) await waitForWorkers();
        const label = `${engine.name} n${threads} ${search.cases.join(",")}`;
        if (search.expect === "failure") {
          checkFailure(run, expected, search, label);
          results.push({ engine: engine.name, threads, program, cases: search.cases, argv: search.argv, formats: {}, frozen: {},
                         failure: run.error });
          console.log(`ok ${label} (fails as the CLI)`);
          continue;
        }
        assert.equal(run.status, 0, `${label}: ${run.error}`);
        const row = { engine: engine.name, threads, program, cases: search.cases, argv: search.argv, outfmt: search.outfmt, formats: {}, frozen: {} };
        for (const format of FORMATS[program]) {
          const actual = run.streams.get(format) || Buffer.alloc(0);
          const want = expected.outputs.get(format);
          assert.equal(sha256(actual), sha256(want), `${label}: outfmt ${format}`);
          row.formats[format] = sha256(actual);
          if (search.frozen[format]) row.frozen[format] = sha256(toFrozen(actual, search.frozen_prefix)) === search.frozen[format];
        }
        if (search.files) {
          const sent = run.streams.get(CUSTOM_STREAM) || Buffer.alloc(0);
          const entries = storedZipEntries(sent);
          // A standard reader (Python's zipfile) reads it too, every CRC-32 checked.
          const archive = path.join(scratch, "archive.zip");
          fs.writeFileSync(archive, sent);
          const checked = require("node:child_process").spawnSync("python3", ["-I", "-c",
            "import sys, zipfile\nwith zipfile.ZipFile(sys.argv[1]) as z:\n    sys.exit(z.testzip() is not None)", archive]);
          assert.equal(checked.status, 0, `${label}: zipfile reads the archive: ${checked.stderr}`);
          fs.rmSync(archive, { force: true });
          const cli = expected.outputs.get(CUSTOM_STREAM);
          assert.deepEqual(entries.map(([name]) => name), cli.map(([name]) => name), `${label}: -outfmt ${search.outfmt} files`);
          for (const [index, [name, bytes]] of entries.entries()) {
            assert.equal(sha256(bytes), sha256(cli[index][1]), `${label}: -outfmt ${search.outfmt} ${name}`);
          }
          row.formats.custom = sha256(sent);
          const master = sha256(toFrozen(entries[entries.length - 1][1], search.frozen_prefix));
          const sides = entries.slice(0, -1).map(([, bytes]) => sha256(toFrozen(bytes, search.frozen_prefix))).join(",");
          if (search.frozen.custom) row.frozen.custom = master === search.frozen.custom;
          if (search.frozen.sides !== undefined) row.frozen.sides = sides === search.frozen.sides;
        } else if (search.outfmt !== undefined) {
          const sent = run.streams.get(CUSTOM_STREAM) || Buffer.alloc(0);
          const actual = search.sam ? withoutCl(sent) : sent;
          const cli = expected.outputs.get(CUSTOM_STREAM);
          assert.equal(sha256(actual), sha256(search.sam ? withoutCl(cli) : cli), `${label}: -outfmt ${search.outfmt}`);
          row.formats.custom = sha256(actual);
          if (search.frozen.custom) row.frozen.custom = sha256(toFrozen(actual, search.frozen_prefix)) === search.frozen.custom;
        } else {
          assert.ok(!run.streams.has(CUSTOM_STREAM), `${label}: no stream ${CUSTOM_STREAM} without -outfmt`);
        }
        assert.equal(sha256(run.streams.get(3) || Buffer.alloc(0)), sha256(expected.stderr), `${label}: diagnostics`);
        row.diagnostics = sha256(expected.stderr);
        row.hits = WITH_HITS.has(program) ? checkHits(run.streams, program, label, queryIds) : null;
        results.push(row);
        console.log(`ok ${label}`);
      }
      assert.equal(reactor.call("losat_web2_release", q.status).status, 0);
    }
    for (const handle of subjects.values()) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
    for (const { handle } of counts.values()) assert.equal(reactor.call("losat_web2_file_release", handle).status, 0);
    if (countsRegistered + countsKept) console.log(`ok ${engine.name} -window_masker_db: ${countsRegistered} counts registrations, ${countsKept} searches with counts kept or unregistered`);
    await close();
  }
  if (options.out) fs.writeFileSync(path.join(options.out, "v-abi-results.json"), JSON.stringify(results, null, 2) + "\n");
  const frozen = results.flatMap((row) => Object.entries(row.frozen).map(([format, same]) => ({ ...row, format, same })));
  const differing = [...new Set(frozen.filter((row) => !row.same).map((row) => `${row.cases.join(",")} outfmt ${row.format}`))];
  console.log(`${results.length} runs matched the native CLI; frozen hashes: ${frozen.length - frozen.filter((row) => !row.same).length}/${frozen.length} equal${differing.length ? `, differing: ${differing.join("; ")}` : ""}`);
  fs.rmSync(scratch, { recursive: true, force: true });
}

main().catch((error) => {
  console.error(error);
  process.exit(1);
});
