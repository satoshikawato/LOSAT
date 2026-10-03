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
// exists (Gate A, the Stage G matrix) the match is recorded as well.
//
// Usage: node v_abi.js --native LOSAT --serial losat-web-serial.wasm
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

const FORMATS = { blastp: [0, 6, 7], tblastn: [0, 6, 7], blastn: [0, 6, 7], tblastx: [0, 6, 7] };
const WITH_HITS = new Set(["blastp", "tblastn", "blastn", "tblastx"]);

function sha256(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

function parseArgs(argv) {
  const options = {};
  for (let i = 0; i < argv.length; i += 2) options[argv[i].replace(/^--/, "")] = argv[i + 1];
  for (const required of ["native", "serial", "threads", "cases"]) {
    if (!options[required]) throw new Error(`--${required} is required`);
  }
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

  emit(stream, ptr, len) {
    const chunk = Buffer.from(new Uint8Array(this.memory().buffer, ptr, len));
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
    return Buffer.from(new Uint8Array(this.memory().buffer, e.losat_web2_last_error_ptr(), e.losat_web2_last_error_len())).toString();
  }

  // Calls an export; Buffer and string arguments become (ptr, len) pairs.
  call(name, ...args) {
    const e = this.exports, allocations = [], flat = [];
    for (const arg of args) {
      if (typeof arg === "number") { flat.push(arg); continue; }
      const bytes = Buffer.isBuffer(arg) ? arg : Buffer.from(arg);
      const ptr = e.losat_web2_alloc(bytes.length);
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

async function openSerial(wasmPath) {
  const wasi = new WASI({ version: "preview1", returnOnExit: true, args: [], env: {} });
  let reactor;
  const imports = {
    wasi_snapshot_preview1: wasi.wasiImport,
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

function nativeOutputs(native, search, scratch) {
  const outputs = new Map();
  let stderr = null;
  for (const format of FORMATS[search.program]) {
    const out = path.join(scratch, `native.${format}.out`);
    const argv = [...search.argv, "-outfmt", String(format), "-num_threads", "1", "-out", out];
    const result = require("node:child_process").spawnSync(native, argv, { cwd: search.cwd, maxBuffer: 1 << 30 });
    if (result.status !== 0) throw new Error(`native ${argv.join(" ")} failed: ${result.stderr}`);
    outputs.set(format, fs.readFileSync(out));
    fs.rmSync(out);
    if (stderr === null) stderr = result.stderr;
    else assert.equal(sha256(result.stderr), sha256(stderr), "the CLI writes the same warnings for every format");
  }
  return { outputs, stderr };
}

function checkHits(streams, program, label) {
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
    if (program === "tblastn") assert.ok(description.subject_gencodes.includes(32), "tblastn: genetic code 32");
  }
  assert.equal(reactor.call("losat_web2_describe", "blastx").status, -1, "blastx joins in SX");
  // validate: accepted argv, and CLI errors with the CLI's words.
  assert.equal(reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s"].join("\0")).status, 0);
  for (const argv of [["blastp", "-query", "q", "-subject", "s", "-evalue", "abc"],
                      ["tblastn", "-query", "q", "-subject", "s", "-db_gencode", "7"],
                      ["blastn", "-query", "q", "-subject", "s", "-nosuchoption", "1"]]) {
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, argv.join(" "));
    const cli = require("node:child_process").spawnSync(native, argv);
    assert.equal(validated.error, cli.stderr.toString(), `${argv.join(" ")}: the CLI's message`);
  }
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
  assert.match(reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-outfmt", "6"].join("\0")).error, /not accepted/);
  // register and scan agree on a multi-record input, in chunks of any size.
  const input = fs.readFileSync(path.join(ROOT, "docs/evidence/tlosan_stage_c/multi_query_20260924/query.faa"));
  const registered = reactor.call("losat_web2_register", "tblastn", 0, input);
  assert.ok(registered.status > 0, registered.error);
  const records = JSON.parse(registered.streams.get(2).toString()).records;
  const scanner = reactor.exports.losat_web2_scan_begin(0);
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
  const bad = reactor.call("losat_web2_register", "blastp", 0, Buffer.from("ACGT\n>x\nAC\n"));
  assert.equal(bad.status, -1);
  assert.match(bad.error, /Expected > at record start/);
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
      for (const threads of engine.threads) {
        const argv = [...search.argv, "-num_threads", String(threads)].join("\0");
        const run = reactor.call("losat_web2_run", argv, q.status, subjects.get(subjectKey));
        if (waitForWorkers) await waitForWorkers();
        const label = `${engine.name} n${threads} ${search.cases.join(",")}`;
        assert.equal(run.status, 0, `${label}: ${run.error}`);
        const row = { engine: engine.name, threads, program, cases: search.cases, argv: search.argv, formats: {}, frozen: {} };
        for (const format of FORMATS[program]) {
          const actual = run.streams.get(format) || Buffer.alloc(0);
          const want = expected.outputs.get(format);
          assert.equal(sha256(actual), sha256(want), `${label}: outfmt ${format}`);
          row.formats[format] = sha256(actual);
          if (search.frozen[format]) row.frozen[format] = sha256(actual) === search.frozen[format];
        }
        assert.equal(sha256(run.streams.get(3) || Buffer.alloc(0)), sha256(expected.stderr), `${label}: diagnostics`);
        row.diagnostics = sha256(expected.stderr);
        row.hits = WITH_HITS.has(program) ? checkHits(run.streams, program, label) : null;
        results.push(row);
        console.log(`ok ${label}`);
      }
      assert.equal(reactor.call("losat_web2_release", q.status).status, 0);
    }
    for (const handle of subjects.values()) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
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
