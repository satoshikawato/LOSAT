"use strict";
// V-ABI (plan §6.2): runs the two adapter reactors under Node through ABI v2
// (docs/web/abi_v2.md) and compares every output format with the native CLI of the same
// commit run for that format alone, by SHA-256. The standard error of the CLI is
// compared with stream 3, and the HSP records (stream 1) are checked against the
// outfmt 6 rows and outfmt 0 sections they point to.
//
// One subject is registered once per program and reused by consecutive runs with other
// queries and options (plan §4.6). The serial reactor runs with 1 thread, the threaded
// reactor with 1, 2 and 4 (the host appends -num_threads, docs/web/abi_v2.md §7). A search whose
// *start* JSON says `register_inputs` (-parse_deflines) reads its query and subject together for
// each run's argv (`losat_web2_register_inputs`, abi_v2.md §4; line WQA's sweep, D7, 2026-11-19),
// and the IDs that the response shows for the records are the first two columns of outfmt 6.
//
// The searches come from web/adapter/tools/v_abi_cases.py; where a frozen SHA-256
// exists (Gate A, the Stage G matrix, NCBI's fixtures) the match is recorded as well, for a
// search with `"frozen_prefix"` after the frozen path prefix is put back (`toFrozen`). A
// search marked `"expect": "failure"` must fail in the CLI's run of the argv's own -outfmt (its
// `"outfmt"`, else outfmt 0) and in the ABI with the same error, and stream 3 then the error must
// be that run's standard error, byte for byte (abi_v2.md §5, §7; line WQA's sweep, D3 review,
// 2026-11-19). The CLI's runs of the other formats may complete: in a run of two searches
// (`-num_descriptions`, `-num_alignments`) the outcome is that of the argv's own search, and the
// other writes formats that the argv does not ask for (§7). A search that completes completes in
// the CLI's run of its own -outfmt, and the CLI's runs of its other formats may fail, whose streams
// are then empty (the run discards them, stream 9); its stream 3 is the CLI's standard error for
// the argv followed by what the other formats' standard error adds (§5, `nativeOutputs`; line
// WQA's sweep, D3 review r2, 2026-11-19). A search with `"outfmt"` 0, 6 or 7 alone runs with that
// -outfmt, which is one of the program's formats; one with another `"outfmt"` (a format besides
// 0, 6 and 7, or 6 or 7 with fields) runs
// with that -outfmt: its output, stream 4, is compared with the CLI's run of the same argv,
// whose standard error then starts stream 3 (abi_v2.md §5, §7); SAM's (`"sam": true`) without the
// CL tag of its @PG line, where each writes its own command line. The output of outfmt 13 and 14
// (`"files": true`) is a ZIP archive, whose entries are compared with the files of the CLI's run
// with `-out results`: the side files `results_<n>.json|xml`, and the master file as
// `results.json|xml`.
//
// A search with `"fails_alone"` lists the formats of 0, 6 and 7 whose CLI run of the argv fails
// while the argv's own format runs (a record that is both a query and a subject under
// -parse_deflines and has masks: LOSAT refuses outfmt 0 alone; line WQA's sweep, D7, 2026-11-19),
// and their CLI runs must fail. As for every format whose CLI run fails in a run that completes, the
// run writes no stream for them, stream 3 adds their error after the warnings (`nativeOutputs`; the
// adapter writes a refusal of a format that failed alone, abi_v2.md §7), and the run's signals of
// stream 9 name them, each with the place of its error, which holds it (abi_v2.md §5). Its `"outfmt"` may be 6 or 7 (the argv's own
// format then is that stream, and no stream 4 is written). Every run that succeeds says once where its messages go
// (stream 10): stream 3 without -logfile (line WQA's sweep, D7, 2026-11-20).
//
// A BLASTN search with `-window_masker_db` names a counts file: the host registers that file under
// the argv's value (`file_begin`, `file_chunk`, `file_end`; abi_v2.md §4) before the run, and keeps
// it for the next searches that name it, as the application keeps one (abi_v2.md §6); a value that
// is not a regular file is not registered, which the run reports as the CLI reports a file it
// cannot open. A search with `"window_masker_path"` (BLASTN `-window_masker_taxid`) runs the CLI
// with that WINDOW_MASKER_PATH, and the host registers the directory's counts files under LOSAT
// Web's directory (`registerTaxidCounts`).
//
// A search with `"database"` (`-db`) has no subject: the run takes subject handle 0 and the engine opens
// the database and the lists of its filters by the argv's names in the run's working directory, which
// the host preopens as `/work`, the reactor's current directory (abi_v2.md §3): the search's directory, as
// the CLI runs there (line WQA D12a, 2026-11-20). Each engine runs it on an instance whose `/work` is that
// directory; the threaded engine's
// threads share the table of the files it opens (LOSAT/tests/wasi_files.js), as the browser host's
// threads share theirs (web/packages/losat/src/infra/reactor/run-directory.ts). A search with `"run_directory"`
// (BLASTN `-filtering_db` with `-subject`, line RG item D4) has a subject as every search has and runs on the
// instance of its directory too, where the engine opens the repeat database by the option's name; its subject
// handle is held and released on that instance.
//
// The words with which the command line prints its help or its version (`checkPrints`; line WQA's D14, 2026-11-20) run
// for every program and word with the host's -num_threads: *start* says what the run prints, `validate` accepts the argv,
// and `run`, given no inputs (handles 0), writes the CLI's standard output on stream 11 and its standard error on stream
// 3, and no output format; the words around them, the registry file and the file of -logfile as the CLI reads them.
//
// A search with `"tool"` is a run of a tool of the command line (`losat windowmasker <argv>`; line NP WW3a, 2026-11-22):
// `"run"` calls `losat_web2_tool` on an instance whose `/work` is a folder of the search's `"files"`, and runs the
// native CLI in another folder of the same files with the standard input at the null device; the exit status
// (stream 2, `status`), the standard output (stream 11), the standard error (stream 3) and every file the run wrote
// (the frames of stream 12, abi_v2.md §5, by name) are compared with the CLI's, and the folder of the CLI's run has no
// file besides those that the frames name; a run of `"expect": "refusal"` is LOSAT's refusal of an output name that LOSAT
// Web does not write (exit status 1, the refusal on stream 3, no output, no file), which the CLI writes. `"validate"`
// calls `losat_web2_validate` on the argv: valid -> 0 and nothing emitted, else -1 and the CLI's standard error then the
// file of -logfile without the last line end (abi_v2.md §4); `"expect": "late_failure"` is an argv that the CLI's run
// refuses after it opened a file, which validate accepts (its checks stop before the files, D-WW3-10).
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
const { directoryReads } = require(path.join(ROOT, "LOSAT/tests/wasi_files"));

const FORMATS = { blastp: [0, 6, 7], tblastn: [0, 6, 7], blastn: [0, 6, 7], tblastx: [0, 6, 7], blastx: [0, 6, 7] };
// The stream of the output of the argv's -outfmt.
const CUSTOM_STREAM = 4;
// The signals that discard what the host received of an output stream (abi_v2.md §5).
const DISCARD_STREAM = 9;
const WITH_HITS = new Set(["blastp", "tblastn", "blastn", "tblastx", "blastx"]);

function sha256(bytes) {
  return crypto.createHash("sha256").update(bytes).digest("hex");
}

// The SHA-256 of LOSAT's output read back in NCBI's words, to compare with a frozen hash of NCBI's
// bytes: LOSAT states its own version and error prefixes (AGENTS.md rule 8), which
// LOSAT/tests/losat_wording.py `hash-form` (`ncbi_hash_form`) reads back exactly. An output that
// still writes NCBI's words has no such form and matches no frozen hash.
function ncbiHash(bytes) {
  const read = require("node:child_process").spawnSync("python3", ["-I", path.join(ROOT, "LOSAT/tests/losat_wording.py"), "hash-form", "-"],
                                                        { input: bytes, maxBuffer: 1 << 30 });
  if (read.error) throw read.error;
  return read.status === 0 ? sha256(read.stdout) : `ncbi-words:${read.stderr}`;
}

// A SAM document without the CL tag of its @PG line (LOSAT/tests/db_search_fixtures.py
// `strip_sam_cl`: from the `\tCL:` after the version of `\n@PG\tID:0\tVN:<version>` up to the
// next `\tPN:blastn\n`).
function withoutCl(document) {
  const head = Buffer.from("\n@PG\tID:0\tVN:");
  const at = document.indexOf(head);
  if (at < 0) return document;
  const tab = document.indexOf(Buffer.from("\t"), at + head.length);
  if (tab < 0 || !document.subarray(tab, tab + 4).equals(Buffer.from("\tCL:"))) return document;
  const start = tab;
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
  // A run that failed leaves the master file (`results`) empty, after the side files written before the error.
  if (fs.existsSync(path.join(directory, "results"))) files.push([`results.${extension}`, fs.readFileSync(path.join(directory, "results"))]);
  return files;
}


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
    this.chunks = new Map();
  }

  // Pointers and lengths are unsigned; JavaScript receives an i32 as a signed number
  // (docs/web/abi_v2.md §6), so a pointer of 2 GiB or more arrives negative.
  // A signal of stream 9 discards what the host received of the stream that it names (§5), as the
  // application discards it: the stream of a run is then empty. The signals are kept on stream 9.
  emit(stream, ptr, len) {
    const chunk = Buffer.from(new Uint8Array(this.memory().buffer, ptr >>> 0, len >>> 0));
    if (stream === DISCARD_STREAM) this.streams.delete(chunk.readUInt32LE(0));
    if (!this.streams.has(stream)) this.streams.set(stream, []);
    this.streams.get(stream).push(chunk);
  }

  take() {
    const result = new Map();
    for (const [stream, chunks] of this.streams) result.set(stream, Buffer.concat(chunks));
    this.chunks = new Map(this.streams); // the emits of the last call, one entry each (the frames of stream 12)
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

// The registry file of a search's -conffile (`start`'s `conffile`), which the host registers as kind 2
// before it validates and registers the inputs (abi_v2.md §4): the file of that name in the search's
// directory, read as the CLI reads it; none where there is no such file, whose run fails as the CLI's
// (`Registry file "x.ini" not found`). Returns its handle, or null.
function registerRegistry(reactor, search, name) {
  const file = path.resolve(search.cwd, name);
  if (!fs.statSync(file, { throwIfNoEntry: false })?.isFile()) return null;
  const bytes = fs.readFileSync(file);
  const begun = reactor.call("losat_web2_file_begin", 2, name, bytes.length);
  assert.ok(begun.status > 0, `file_begin ${name}: ${begun.error}`);
  assert.equal(reactor.call("losat_web2_file_chunk", begun.status, bytes).status, 0);
  const ended = reactor.call("losat_web2_file_end", begun.status);
  assert.equal(ended.status, 0, `file_end ${name}: ${ended.error}`);
  assert.equal(JSON.parse(ended.streams.get(2).toString()).kind, "registry");
  return begun.status;
}

// The counts files of the searches: the argv's `-window_masker_db` value -> its handle and path.
// A search that names a registered value with another file releases it first.
function registerCounts(reactor, counts, search) {
  const at = search.argv.indexOf("-window_masker_db");
  if (at < 0) return null;
  const name = search.argv[at + 1];
  const file = path.resolve(search.cwd, name);
  const known = counts.get(name);
  if (known && known.file === file) return { name, held: "kept" };
  if (known) {
    assert.equal(reactor.call("losat_web2_file_release", known.handle).status, 0);
    counts.delete(name);
  }
  if (name === "" || !fs.existsSync(file) || !fs.statSync(file).isFile()) return { name, held: "absent" };
  const { handle } = registerFile(reactor, name, file);
  counts.set(name, { handle, file });
  return { name, held: "registered" };
}

// LOSAT Web's directory of -window_masker_taxid (web/adapter/src/files.rs WINDOW_MASKER_PATH).
const WEB_WINDOW_MASKER_PATH = "windowmasker";

// The counts files of a search's -window_masker_taxid directory (`"window_masker_path"`, the CLI's
// WINDOW_MASKER_PATH): every readable regular file <path>/<n>/wmasker.obinary, registered as
// windowmasker/<n>/wmasker.obinary and kept across runs (a file the CLI cannot read, or a directory
// of that name, is not registered: LOSAT Web has neither). Returns null without the directory.
function registerTaxidCounts(reactor, counts, search) {
  if (search.window_masker_path === undefined) return null;
  const root = path.resolve(search.cwd, search.window_masker_path);
  let registered = 0, kept = 0;
  for (const entry of fs.readdirSync(root, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue;
    const file = path.join(root, entry.name, "wmasker.obinary");
    try {
      if (!fs.statSync(file).isFile()) continue;
      fs.accessSync(file, fs.constants.R_OK);
    } catch {
      continue;
    }
    const name = `${WEB_WINDOW_MASKER_PATH}/${entry.name}/wmasker.obinary`;
    const known = counts.get(name);
    if (known && known.file === file) {
      kept++;
      continue;
    }
    if (known) {
      assert.equal(reactor.call("losat_web2_file_release", known.handle).status, 0);
      counts.delete(name);
    }
    const { handle } = registerFile(reactor, name, file);
    counts.set(name, { handle, file });
    registered++;
  }
  return { name: WEB_WINDOW_MASKER_PATH, held: registered ? "registered" : "kept" };
}

// The serial reactor; `directory`, the run's working directory that it preopens as `/work` (a search of -db).
async function openSerial(wasmPath, directory) {
  const wasi = new WASI({ version: "preview1", returnOnExit: true, args: [], env: {}, ...(directory === undefined ? {} : { preopens: { "/work": directory } }) });
  let reactor, memory;
  // A folder of the run's directory opened for reading fails its reads with EISDIR, as on Linux and in the browser
  // host (directoryReads, LOSAT/tests/wasi_files.js).
  const guestMemory = { get buffer() { return memory.buffer; }, grow: (pages) => memory.grow(pages) };
  const imports = {
    wasi_snapshot_preview1: directoryReads(unsignedWasiImports(wasi.wasiImport), guestMemory, directory ?? "/"),
    losat_host: { emit: (stream, ptr, len) => reactor.emit(stream, ptr, len) },
  };
  const { instance } = await WebAssembly.instantiate(fs.readFileSync(wasmPath), imports);
  memory = instance.exports.memory;
  wasi.initialize(instance);
  reactor = new Reactor(instance.exports, () => instance.exports.memory);
  return { reactor, close: async () => {} };
}

// The threaded reactor; `directory` as openSerial's.
async function openThreads(wasmPath, directory) {
  let reactor;
  const host = await createThreadHost(wasmPath, "threaded-reactor", [], {
    directory,
    imports: { losat_host: { emit: (stream, ptr, len) => reactor.emit(stream, ptr, len) } },
  });
  reactor = new Reactor(host.instance.exports, () => host.memory);
  return { reactor, close: () => host.close(), waitForWorkers: () => host.waitForWorkers() };
}

// ---------------------------------------------------------------------------

function inputPath(search, option) {
  return path.resolve(search.cwd, search.argv[search.argv.indexOf(option) + 1]);
}

// The CLI's environment of a search: the harness's, with the search's directory of
// -window_masker_taxid as NCBI's WINDOW_MASKER_PATH (and without one inherited otherwise).
function nativeEnv(search) {
  const env = { ...process.env };
  delete env.WINDOW_MASKER_PATH;
  if (search.window_masker_path !== undefined) env.WINDOW_MASKER_PATH = search.window_masker_path;
  return env;
}

// The hit list size that the CLI derives for the output format `number` from the options of `argv`
// (`report_counts`, NCBI blast_args.cpp:2894-2973): the formats of a run with the same size are one
// search (abi_v2.md §7).
function hitListSize(argv, number) {
  const value = (flag) => {
    const at = argv.lastIndexOf(flag);
    return at < 0 ? undefined : Number(argv[at + 1]);
  };
  const [max, descriptions, alignments] = ["-max_target_seqs", "-num_descriptions", "-num_alignments"].map(value);
  if (max !== undefined) return max;
  if (number > 4) return alignments ?? 500;
  return descriptions === undefined && alignments === undefined ? 500 : Math.max(descriptions ?? 500, alignments ?? 250);
}

// The lines of `text` (each with its line end) that it has more times than `held` (line -> copies),
// in order; `held` then holds them (the adapter's `Diagnostics`). The adapter compares the error
// of a format that failed alone as one text (`Diagnostics::error`), and this compares its lines
// too: the two agree while no line of an error is a line that stream 3 already holds from another
// text. NCBI's exception reports and posts end with an empty line, so a case where two such texts
// meet (BLASTX's HSPs of length 0 with -outfmt 5: the post of the XML and the refusals of 6, 7 and
// 0, the adapter's test `formats_that_stop_after_a_query_warning_add_it_once`) is not one for this
// harness (reviews r2 L2 and r3 L5 of 2026-11-20).
function linesAfter(held, text) {
  const copies = new Map();
  let added = "";
  for (const line of text.split(/(?<=\n)/)) {
    if (line === "") continue;
    const copy = (copies.get(line) || 0) + 1;
    copies.set(line, copy);
    if (copy > (held.get(line) || 0)) added += line;
  }
  for (const [line, copy] of copies) held.set(line, Math.max(held.get(line) || 0, copy));
  return added;
}

// The formats `formats` in the order of a search's formats after its first: the tabular ones, then
// those up to 4 (abi_v2.md §7).
function searchOrder(formats) {
  return [...formats.filter((format) => format > 4), ...formats.filter((format) => format <= 4)];
}

// The CLI's outputs of a search, one run per format with -out. A search with `"expect": "failure"`
// fails in the CLI's run of the argv's own -outfmt (exit 3 or 1), whose files keep what was
// written before the failure (the outfmt 0 prolog before a first query batch fails); the runs of
// its other formats may complete (a run of two searches, abi_v2.md §7). A search that completes
// completes in the CLI's run of its own -outfmt, and the runs of its other formats may fail (the
// other search of a run of two searches, a format that fails alone, §7), whose outputs the run
// then discards: they are empty (line WQA's sweep, D3 review r2, 2026-11-19). `own` is the CLI's
// standard error for the run's argv: that of the search's -outfmt, else that of outfmt 0.
// `stderr` is stream 3 of a run that completes (§5): `own`, then the lines of the other formats'
// standard error that it does not hold, each line as many times as one of them has it more often
// than stream 3 holds it, the formats of the argv's own search first (the tabular ones, then those
// up to 4). A structured -outfmt's run does not write the queries' warnings, which the other
// formats' runs add (line WQA's sweep, D3 review r3 M2, 2026-11-19: this harness kept only the
// formatting handler's lines of the other formats there, and so took their loss), as they add the
// error of a format that fails alone (LOSAT's refusal of outfmt 0 of a record that is both a query
// and a subject and has masks, line WQA's sweep, D7, 2026-11-19).
function nativeOutputs(native, search, scratch) {
  const outputs = new Map();
  const stderrs = new Map(), failed = new Map();
  const failure = search.expect === "failure";
  // The argv's own format: 0, 6 or 7 alone is one of the program's formats, any other -outfmt is
  // the output of its own (`null`).
  const outfmt = search.outfmt === undefined ? "0" : search.outfmt.trim();
  const ownFormat = /^[067]$/.test(outfmt) ? Number(outfmt) : null;
  // The formats that the search says fail alone (`"fails_alone"`), whose CLI run must fail.
  const alone = new Set(search.fails_alone || []);
  let status = null;
  for (const format of FORMATS[search.program]) {
    const out = path.join(scratch, `native.${format}.out`);
    const argv = [...search.argv, "-outfmt", String(format), "-num_threads", "1", "-out", out];
    const result = require("node:child_process").spawnSync(native, argv, { cwd: search.cwd, env: nativeEnv(search), maxBuffer: 1 << 30 });
    const fails = result.status !== 0;
    if (format === ownFormat) {
      if (fails !== failure) {
        throw new Error(`native ${argv.join(" ")} ${fails ? "failed" : "succeeded"} (exit ${result.status}): ${result.error ?? result.stderr}`);
      }
      status = result.status;
    } else if (alone.has(format) && !fails) {
      throw new Error(`native ${argv.join(" ")} succeeded, where "fails_alone" says that it fails`);
    }
    const file = fs.existsSync(out) ? fs.readFileSync(out) : Buffer.alloc(0);
    outputs.set(format, fails && !failure ? Buffer.alloc(0) : file);
    fs.rmSync(out, { force: true });
    stderrs.set(format, result.stderr);
    failed.set(format, fails);
  }
  let own = ownFormat === null ? null : stderrs.get(ownFormat);
  if (ownFormat === null) {
    const directory = fs.mkdtempSync(path.join(scratch, "native.custom."));
    const out = path.join(directory, search.files ? "results" : "native.custom.out");
    const argv = [...search.argv, "-outfmt", search.outfmt, "-num_threads", "1", "-out", out];
    const result = require("node:child_process").spawnSync(native, argv, { cwd: search.cwd, env: nativeEnv(search), maxBuffer: 1 << 30 });
    if ((result.status !== 0) !== failure) {
      throw new Error(`native ${argv.join(" ")} ${failure ? "succeeded" : "failed"} (exit ${result.status}): ${result.error ?? result.stderr}`);
    }
    // The file of a CLI run that failed holds what it wrote before the error (line WQA's sweep, D15, 2026-11-20).
    if (!failure || fs.existsSync(out) || search.files) {
      outputs.set(CUSTOM_STREAM, search.files ? nativeFiles(directory, search.outfmt) : (fs.existsSync(out) ? fs.readFileSync(out) : Buffer.alloc(0)));
    }
    fs.rmSync(directory, { recursive: true, force: true });
    [own, status] = [result.stderr, result.status];
  }
  let stderr = null;
  // The streams of the formats whose CLI run fails while the run completes, which the run's signals
  // of stream 9 name (abi_v2.md §5): the other search's, or a format's that fails alone (§7).
  let failedStreams = [];
  if (!failure) {
    const number = parseInt(outfmt, 10);
    const size = hitListSize(search.argv, number);
    const others = FORMATS[search.program].filter((format) => format !== ownFormat);
    const same = searchOrder(others.filter((format) => hitListSize(search.argv, format) === size));
    const other = searchOrder(others.filter((format) => hitListSize(search.argv, format) !== size));
    const held = new Map();
    // The CLI's standard error as bytes (latin1 keeps each byte): a reader warning may hold bytes that are not
    // UTF-8 (`Unknown gap-type: caf\xe9`, a `>?` gap line's mods; line WQA's D13, 2026-11-20).
    let text = linesAfter(held, own.toString("latin1"));
    for (const format of [...same, ...other]) text += linesAfter(held, stderrs.get(format).toString("latin1"));
    stderr = Buffer.from(text, "latin1");
    failedStreams = others.filter((format) => failed.get(format)).sort((a, b) => a - b);
  }
  return { outputs, stderr, status, own, failedStreams };
}

// Whether the search's `"outfmt"` writes stream 4: a format besides 0, 6 and 7, or 6 or 7 with
// fields (an argv's own plain 0, 6 or 7 is that stream).
function customStream(search) {
  return search.outfmt !== undefined && !/^\s*[067]\s*$/.test(search.outfmt);
}

// An error of `run` or `validate` as the CLI writes it (abi_v2.md §4, §5): NCBI's messages as they
// are (with their line end), but for those of the application start, which come without their
// last line end, and which begin with NCBI's severity and error code (`Critical: (106.23)`,
// `Error: (CRegistryException::eSection)`); LOSAT's own errors after `Error: ` and with a line end
// (line WQA's sweep, D3 review r2 L2, 2026-11-19).
function cliError(error) {
  if (error.endsWith("\n")) return error;
  if (/^(Critical|Error): \(/.test(error)) return `${error}\n`;
  return `Error: ${error}\n`;
}

// A run that fails as the CLI fails: the CLI's error is the last message of its standard
// error, and the warnings that the CLI writes before it are stream 3, which a failed run sends
// (abi_v2.md §5): stream 3 then the error are the CLI's standard error for the run's argv, byte
// for byte (line WQA's sweep, D3, 2026-11-19), also in a run of two searches, whose outcome is
// that of the argv's own search (§7). The output of the argv's own -outfmt (stream 4 for one of its
// own, else the stream of its format, 0 without one) is exactly the CLI's file for the argv: what
// it wrote before the error (version 7; line WQA's sweep, D15, 2026-11-20). Every other output of
// the run is dropped, each with a signal of stream 9 whose place is 0xFFFFFFFE, the error that ended
// the run (the run failed before its search made them: none, and then they are empty). A search with
// `"prolog": false` fails before any output (an error of the options or of their files), so the CLI's
// outfmt 0 file is empty; the others fail after the CLI wrote the outfmt 0 prolog, or the CLI's
// outfmt 0 completes (the other search of a run whose own -outfmt is a tabular one).
function checkFailure(run, expected, search, label) {
  const program = search.program;
  assert.equal(run.status, -1, `${label}: the run fails as the CLI (exit ${expected.status})`);
  assert.ok(run.error.length > 0 && expected.own.toString().endsWith(cliError(run.error)), `${label}: the CLI's error (${run.error})`);
  const diagnostics = run.streams.get(3) || Buffer.alloc(0);
  assert.equal(diagnostics.toString() + cliError(run.error), expected.own.toString(),
               `${label}: stream 3 and the error are the CLI's standard error`);
  const outfmt = search.outfmt === undefined ? "0" : search.outfmt.trim();
  const own = customStream(search) ? CUSTOM_STREAM : Number(outfmt);
  const sent = run.streams.get(own) || Buffer.alloc(0);
  const cli = expected.outputs.get(own);
  if (search.files) {
    // The side files written before the error and the master file as the CLI leaves it (empty), with a
    // central directory (D15 review round 2, 2026-11-21).
    const entries = sent.length === 0 ? [] : storedZipEntries(sent);
    assert.deepEqual(entries.map(([name]) => name), cli.map(([name]) => name), `${label}: -outfmt ${search.outfmt} files before the failure`);
    for (const [index, [name, bytes]] of entries.entries()) assert.equal(sha256(bytes), sha256(cli[index][1]), `${label}: ${name}`);
  } else if (search.sam) {
    assert.equal(sha256(withoutCl(sent)), sha256(withoutCl(cli)), `${label}: -outfmt ${search.outfmt} before the failure`);
  } else {
    assert.equal(sha256(sent), sha256(cli), `${label}: the bytes of ${customStream(search) ? `-outfmt ${search.outfmt}` : `outfmt ${own}`} are the CLI's file (${sent.length} of ${cli.length} bytes)`);
  }
  const signals = run.streams.get(9) || Buffer.alloc(0);
  assert.equal(signals.length % 8, 0, `${label}: signals of 8 bytes (stream 9)`);
  const dropped = Array.from({ length: signals.length / 8 }, (_, i) => {
    assert.equal(signals.readUInt32LE(8 * i + 4), 0xfffffffe, `${label}: the place of the error of a dropped output is the run's`);
    return signals.readUInt32LE(8 * i);
  }).sort((a, b) => a - b);
  const others = FORMATS[program].filter((format) => format !== own);
  if (dropped.length > 0) assert.deepEqual(dropped, others, `${label}: the outputs that the failed run drops (stream 9)`);
  for (const format of others) assert.ok(!run.streams.has(format), `${label}: nothing of the dropped outfmt ${format}`);
  assert.ok(!run.streams.has(1), `${label}: no HSP records`);
  // `"sends"`: the formats whose report of the batches before the failure is over the adapter's chunk (1 MiB), which the
  // run sends before it fails (review L3, 2026-11-20).
  for (const format of search.sends || []) {
    assert.ok((run.streams.get(format) || Buffer.alloc(0)).length > 0, `${label}: outfmt ${format} sent before the failure`);
  }
  if (search.prolog === false) {
    assert.equal(expected.outputs.get(0).length, 0, `${label}: the CLI writes no outfmt 0 before the failure`);
  } else {
    assert.ok(expected.outputs.get(0).length > 0, `${label}: the CLI's outfmt 0 file has the prolog`);
  }
}

function checkHits(streams, program, label, queryIds, shownIds, database = false, html = false) {
  const records = (streams.get(1) || Buffer.alloc(0)).toString().split("\n").filter(Boolean).map((line) => JSON.parse(line));
  // A stream without bytes (a search without hits) is never emitted.
  const out6 = streams.get(6) || Buffer.alloc(0), out0 = streams.get(0) || Buffer.alloc(0);
  const rows = out6.toString().split("\n").filter(Boolean);
  assert.equal(records.length, rows.length, `${label}: one record per outfmt 6 row`);
  const ranks = new Map(), subjectOf = new Map();
  // BLASTN's matrix-only scoring (`-task rmblastn`, reward 0 and penalty 0) opens its report with
  // RMBLASTN and shows the raw score alone (` Score = 27204`; blast_format.cpp:168-177,
  // showalign.cpp:3586-3591,3605); the other searches show `<bits> bits (<raw>),`.
  const rawScoreOnly = out0.subarray(0, 9).toString() === "RMBLASTN ";
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
    if (shownIds) assert.equal(fields[0], shownIds.query[record.q_idx], `${label}: the shown query ID of HSP ${index}`);
    if (database) {
      // A subject of a database is no registered record: `s_idx` is its OID in the database (abi_v2.md §8),
      // and every row of one OID shows that subject's ID.
      const shown = subjectOf.get(record.s_idx) ?? fields[1];
      assert.equal(fields[1], shown, `${label}: the subject ID of OID ${record.s_idx} at HSP ${index}`);
      subjectOf.set(record.s_idx, shown);
    } else if (shownIds) {
      assert.equal(fields[1], shownIds.subject[record.s_idx], `${label}: the shown subject ID of HSP ${index}`);
    }
    // A subject beyond the alignments that outfmt 0 shows (for example BLASTN's 250) has
    // no section and no heading.
    if (record.out0 === null) {
      assert.equal(record.out0_subject, null, `${label}: out0_subject of HSP ${index} without a section`);
      return;
    }
    const section = out0.subarray(record.out0[0], record.out0[1]).toString();
    // A subject whose title NCBI cannot decode (pd.undec.*, under -parse_deflines) has an empty
    // heading, and NCBI's line stands in place of each alignment (showalign.cpp:1976-1981;
    // LOSAT's report/pairwise.rs write_skipped_alignment; line WQA's sweep, D7, 2026-11-19).
    if (/^Sequence with id \S+ no longer exists in database\.\.\.alignment skipped\n$/.test(section)) {
      assert.equal(record.out0_subject[0], record.out0_subject[1], `${label}: the empty heading of HSP ${index}`);
      assert.ok(record.out0_subject[1] <= record.out0[0], `${label}: out0_subject of HSP ${index}`);
      return;
    }
    assert.ok(section.startsWith(" Score ="), `${label}: out0 of HSP ${index}`);
    const score = rawScoreOnly ? section.startsWith(` Score = ${record.raw_score}\n`) : section.includes(`bits (${record.raw_score}),`);
    assert.ok(score, `${label}: raw score in out0 of HSP ${index}`);
    const heading = out0.subarray(record.out0_subject[0], record.out0_subject[1]).toString();
    // A local id that is not shown leaves "> " before the title; under -parse_deflines the
    // heading starts with the subject's id (`>BK000001.1 TPA: third party`; line WQA's sweep, D7, 2026-11-19),
    // as it does for a subject of a database, whose ids the database holds (line WQA D12a, 2026-11-20).
    // Under -html NCBI writes an anchor first in a heading whose id is not shown (`><a name=Subject_1></a> s1`,
    // the bytes of html.bp.s.f0): the raw prefix is "><a ", not "> " (line RG T4a2, 2026-11-23).
    const prefix = shownIds || database ? ">" : html ? "><a " : "> ";
    assert.ok(heading.startsWith(prefix) && record.out0_subject[1] <= record.out0[0], `${label}: out0_subject of HSP ${index}`);
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
      // The help is NCBI's comment of -task; the tasks are its choices (NCBI's set, in the order of NCBI's Web page).
      assert.equal(task.help, "Task to execute", "blastn: NCBI's help of -task");
      assert.deepEqual(task.choices, ["megablast", "dc-megablast", "blastn", "blastn-short", "rmblastn"], "blastn: tasks in -task's choices");
    }
  }
  // validate: accepted argv, and CLI errors with the CLI's words. An accepted argv says whether its run reads the inputs
  // with NCBI's data loaders on and with -parse_deflines on (the *validation* JSON, abi_v2.md §9; D16 and X7-5).
  const accepted = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s"].join("\0"));
  assert.equal(accepted.status, 0);
  assert.deepEqual(JSON.parse(accepted.streams.get(2).toString()), { data_loaders: true, parse_deflines: false });
  const acceptedPd = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-parse_deflines"].join("\0"));
  assert.equal(acceptedPd.status, 0);
  assert.deepEqual(JSON.parse(acceptedPd.streams.get(2).toString()), { data_loaders: true, parse_deflines: true });
  const acceptedPdDry = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-parse_deflines", "-dryrun"].join("\0"));
  assert.equal(acceptedPdDry.status, 0);
  assert.deepEqual(JSON.parse(acceptedPdDry.streams.get(2).toString()), { data_loaders: true, parse_deflines: false });
  // NCBI's application start (abi_v2.md §4 `losat_web2_start`, line RG E1w-3a): the *start* JSON, the
  // file of -conffile as a registry file (kind 2), and -dryrun.
  const started = reactor.call("losat_web2_start", ["tblastx", "-query", "q", "-subject", "s", "-conffile", "e1w3a.ini", "-logfile", "run.log", "-dryrun"].join("\0"));
  assert.equal(started.status, 0, started.error);
  assert.deepEqual(JSON.parse(started.streams.get(2).toString()), { conffile: "e1w3a.ini", log_file: "run.log", strategy_file: null, dry_run: true, register_inputs: false, prints: null, database: false });
  const parsed = reactor.call("losat_web2_start", ["blastn", "-query", "q", "-subject", "s", "-parse_deflines"].join("\0"));
  assert.equal(JSON.parse(parsed.streams.get(2).toString()).register_inputs, true, "-parse_deflines reads the inputs together");
  const conffile = ["blastp", "-query", "q", "-subject", "s", "-conffile", "e1w3a.ini", "-outfmt", "x", "-dryrun"].join("\0");
  assert.equal(reactor.call("losat_web2_validate", conffile).error, 'Critical: (106.23) Registry file "e1w3a.ini" not found, aborting.');
  const registry = Buffer.from("[BLAST]\n");
  const begunRegistry = reactor.call("losat_web2_file_begin", 2, "e1w3a.ini", registry.length);
  assert.ok(begunRegistry.status > 0, begunRegistry.error);
  assert.equal(reactor.call("losat_web2_file_chunk", begunRegistry.status, registry).status, 0);
  const endedRegistry = reactor.call("losat_web2_file_end", begunRegistry.status);
  assert.equal(JSON.parse(endedRegistry.streams.get(2).toString()).kind, "registry");
  assert.equal(reactor.call("losat_web2_validate", conffile).status, 0, "a dry run reads no -outfmt");
  assert.equal(reactor.call("losat_web2_file_release", begunRegistry.status).status, 0);
  // `[BLAST] DATA_LOADERS=none` of the registry turns the data loaders off (blast_scope_src.cpp:78-80).
  const noLoaders = Buffer.from("[BLAST]\nDATA_LOADERS=none\n");
  const begunNoLoaders = reactor.call("losat_web2_file_begin", 2, "nol.ini", noLoaders.length);
  assert.equal(reactor.call("losat_web2_file_chunk", begunNoLoaders.status, noLoaders).status, 0);
  assert.equal(reactor.call("losat_web2_file_end", begunNoLoaders.status).status, 0);
  const off = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-conffile", "nol.ini"].join("\0"));
  assert.equal(off.status, 0, off.error);
  assert.deepEqual(JSON.parse(off.streams.get(2).toString()), { data_loaders: false, parse_deflines: false });
  const offPd = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-conffile", "nol.ini", "-parse_deflines"].join("\0"));
  assert.equal(offPd.status, 0, offPd.error);
  assert.deepEqual(JSON.parse(offPd.streams.get(2).toString()), { data_loaders: false, parse_deflines: true });
  assert.equal(reactor.call("losat_web2_file_release", begunNoLoaders.status).status, 0);
  for (const argv of [["blastp", "-query", "q", "-subject", "s", "-evalue", "abc"],
                      ["tblastn", "-query", "q", "-subject", "s", "-db_gencode", "7"],
                      ["blastn", "-query", "q", "-subject", "s", "-nosuchoption", "1"],
                      ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast", "-template_type", "coding"],
                      ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast", "-template_type", "coding",
                       "-template_length", "17"]]) {
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, argv.join(" "));
    const cli = require("node:child_process").spawnSync(native, argv);
    assert.equal(validated.error, cli.stderr.toString(), `${argv.join(" ")}: the CLI's message`);
  }
  // BLASTX's ranges (line SX, DW-39): validate reads their syntax, the subject range first.
  const blastx = ["blastx", "-query", "q", "-subject", "s"];
  assert.equal(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "1-9", "-subject_loc", "2-5"].join("\0")).status, 0);
  assert.match(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "9-1", "-subject_loc", "5-2"].join("\0")).error,
               /^LOSAT engine error: Invalid specification of subject location/);
  assert.match(reactor.call("losat_web2_validate", [...blastx, "-query_loc", "9-1"].join("\0")).error,
               /^LOSAT engine error: Invalid specification of query location/);
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
  // The tasks of Session SD and line RG X6, and NCBI's option check of a template
  // (blast_options.c:1247-1261).
  for (const task of ["dc-megablast", "blastn-short", "rmblastn"]) {
    assert.equal(reactor.call("losat_web2_validate", ["blastn", "-query", "q", "-subject", "s", "-task", task].join("\0")).status, 0, task);
  }
  assert.match(reactor.call("losat_web2_validate", ["blastn", "-query", "q", "-subject", "s", "-task", "dc-megablast",
                                                    "-word_size", "13"].join("\0")).error,
               /^LOSAT query\/options error: Invalid discontiguous template parameters: word size must be either 11 or 12/);
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
  // register and the scan of the program's input (the TBLASTN query, NCBI's protein reader) agree on a
  // multi-record input, in chunks of any size.
  const input = fs.readFileSync(path.join(ROOT, "docs/evidence/tlosan_stage_c/multi_query_20260924/query.faa"));
  const registered = reactor.call("losat_web2_register", "tblastn", 0, 1, 0, input);
  assert.ok(registered.status > 0, registered.error);
  const records = JSON.parse(registered.streams.get(2).toString()).records;
  const scanner = reactor.call("losat_web2_scan_begin", "tblastn", 0, 1, 0).status;
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
  // NCBI's protein reader reads a first record without a defline (its ID is empty), and a gap line's letters
  // are its record's, as the CLI reads them (abi_v2.md §4, §9; line WQA's D13, 2026-11-20).
  const headless = reactor.call("losat_web2_register", "blastp", 0, 1, 0, Buffer.from("ACGT\n>x\nAC\n"));
  assert.ok(headless.status > 0, headless.error);
  assert.deepEqual(JSON.parse(headless.streams.get(2).toString()).records.map((r) => [r.id, r.length]),
                   [["", 4], ["x", 2]]);
  assert.equal(reactor.call("losat_web2_release", headless.status).status, 0);
  const gap = reactor.call("losat_web2_register", "blastn", 1, 1, 0, Buffer.from(">s\nACGT\n>?10\nACGT\n"));
  assert.ok(gap.status > 0, gap.error);
  assert.deepEqual(JSON.parse(gap.streams.get(2).toString()).records.map((r) => [r.id, r.length]), [["s", 18]]);
  assert.equal(reactor.call("losat_web2_release", gap.status).status, 0);
  // A first line that may be a Seq-id: the CLI's refusal with the data loaders on, data with them off (D13 (c), D16).
  const seqId = Buffer.from("2HHB\nMKTAYIAKQR\n");
  const refused = reactor.call("losat_web2_register", "blastp", 0, 1, 0, seqId);
  assert.equal(refused.status, -1);
  assert.match(refused.error, /^the first line of the query \("2HHB"\) is not a defline and may be a sequence identifier .* LOSAT's BLASTP /);
  const scanRefused = reactor.call("losat_web2_scan_begin", "blastp", 0, 1, 0).status;
  assert.equal(reactor.call("losat_web2_scan_chunk", scanRefused, seqId).status, 0);
  assert.equal(reactor.call("losat_web2_scan_end", scanRefused).error, refused.error, "the scan refuses in the same words");
  const data = reactor.call("losat_web2_register", "blastp", 0, 0, 0, seqId);
  assert.ok(data.status > 0, data.error);
  assert.deepEqual(JSON.parse(data.streams.get(2).toString()).records.map((r) => [r.id, r.length]), [["", 13]]);
  assert.equal(reactor.call("losat_web2_release", data.status).status, 0);
  // scan_begin and register with parse_deflines=1 read other-molecule Seq-id records with NCBI's alphabet and
  // encoding, matching register_inputs (line RG-XT item X7-5, ABI version 14).
  for (const [program, file, plainLens, pdLens] of [
    ["blastn", "LOSAT/tests/fasta/outfmt0/pd_mol/s_sp_nucx_bn.fa", [287], [300]],
    ["blastp", "LOSAT/tests/fasta/outfmt0/pd_mol/s_ti_prot.fa", [150, 110], [107, 110]],
  ]) {
    const bytes = fs.readFileSync(path.join(ROOT, file));
    for (const [pd, want] of [[0, plainLens], [1, pdLens]]) {
      const s = reactor.call("losat_web2_scan_begin", program, 1, 1, pd).status;
      assert.ok(s > 0, `${program} pd=${pd}`);
      assert.equal(reactor.call("losat_web2_scan_chunk", s, bytes).status, 0);
      const done = reactor.call("losat_web2_scan_end", s);
      assert.equal(done.status, 0, done.error);
      assert.deepEqual(JSON.parse(done.streams.get(2).toString()).records.map((r) => r.length), want, `${file} scan pd=${pd}`);
      const reg = reactor.call("losat_web2_register", program, 1, 1, pd, bytes);
      assert.ok(reg.status > 0, `${file} register pd=${pd}: ${reg.error}`);
      assert.deepEqual(JSON.parse(reg.streams.get(2).toString()).records.map((r) => r.length), want, `${file} register pd=${pd}`);
      assert.equal(reactor.call("losat_web2_release", reg.status).status, 0);
    }
  }
  // register with parse_deflines=1 matches register_inputs on CheckDataLine boundary inputs (Finding 2).
  const spJunk = Buffer.from(">sp|P99987.1|NAME_J\nEEEEEE\n%%%%%%\n");
  const tiJunk = Buffer.from(">gnl|TI|12400\nEEEEEE\n%%%%%%\n");
  assert.equal(reactor.call("losat_web2_register", "blastn", 1, 1, 0, spJunk).status, -1);
  const spJunkPd = reactor.call("losat_web2_register", "blastn", 1, 1, 1, spJunk);
  assert.ok(spJunkPd.status > 0, spJunkPd.error);
  assert.deepEqual(JSON.parse(spJunkPd.streams.get(2).toString()).records.map((r) => r.length), [6]);
  assert.equal(reactor.call("losat_web2_release", spJunkPd.status).status, 0);
  const tiJunkPlain = reactor.call("losat_web2_register", "blastp", 1, 1, 0, tiJunk);
  assert.ok(tiJunkPlain.status > 0, tiJunkPlain.error);
  assert.deepEqual(JSON.parse(tiJunkPlain.streams.get(2).toString()).records.map((r) => r.length), [6]);
  assert.equal(reactor.call("losat_web2_release", tiJunkPlain.status).status, 0);
  assert.equal(reactor.call("losat_web2_register", "blastp", 1, 1, 1, tiJunk).status, -1);
  assert.match(reactor.call("losat_web2_register", "blastp", 1, 1, 2, tiJunk).error, /^parse_deflines must be 0 or 1, not 2$/);
  // validate omits parse_deflines when an argv fails to parse (Finding 3).
  const unparsedVal = reactor.call("losat_web2_validate", ["blastp", "-query", "q", "-subject", "s", "-parse_deflines", "-evalue"].join("\0"));
  assert.equal(unparsedVal.status, -1);
  assert.deepEqual(JSON.parse(unparsedVal.streams.get(2).toString()), { data_loaders: true });
  for (const [args, error] of [[["blastp", 2, 1, 0], /^unknown input role 2$/], [["blastp", 0, 2, 0], /^data_loaders must be 0 or 1, not 2$/],
                               [["blastp", 0, 1, 2], /^parse_deflines must be 0 or 1, not 2$/], [["nosuch", 0, 1, 0], /./]]) {
    const begun = reactor.call("losat_web2_scan_begin", ...args);
    assert.equal(begun.status, -1, args.join(" "));
    assert.match(begun.error, error);
  }
  checkFiles(reactor);
}

// The words with which the command line prints its help or its version (line WQA's D14, 2026-11-20; abi_v2.md §5 stream
// 11, §9 `prints`), for every program and word, and the words, registry file and file of -logfile around them, each
// compared with the native CLI's run of the same argv from an empty directory.
// `--help` is no word of NCBI's reader (`Unknown argument: "-help"`, line RG item K4i-2b, 2026-11-23).
const HELP_WORDS = ["-h", "-help", "-help-full", "-xmlhelp"];
const VERSION_WORDS = ["-version", "-version-full", "-version-full-xml", "-version-full-json"];

function checkPrints(reactor, native) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "losat-v-abi-prints-"));
  const env = Object.fromEntries(Object.entries(process.env).filter(([name]) => !name.startsWith("LOSAT_")));
  const cli = (argv) => require("node:child_process").spawnSync(native, argv, { cwd: dir, env, maxBuffer: 1 << 26 });
  const start = (argv) => {
    const started = reactor.call("losat_web2_start", argv.join("\0"));
    assert.equal(started.status, 0, `${argv.join(" ")}: start: ${started.error}`);
    return JSON.parse(started.streams.get(2).toString());
  };
  const inputs = ["-query", "q.fa", "-subject", "s.fa"];
  // A run that prints `kind`: the CLI's standard output on stream 11, its standard error on stream 3, nothing else.
  const prints = (argv, kind) => {
    const label = argv.join(" ");
    const json = start(argv);
    assert.equal(json.prints, kind, `${label}: what start says the run prints`);
    assert.equal(json.dry_run, false, `${label}: a run that prints is not a dry run`);
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, 0, `${label}: validate: ${validated.error}`);
    const run = reactor.call("losat_web2_run", argv.join("\0"), 0, 0);
    assert.equal(run.status, 0, `${label}: run: ${run.error}`);
    const expected = cli(argv);
    assert.equal(expected.status, 0, `${label}: the CLI's exit status`);
    assert.ok(expected.stdout.length > 0, `${label}: the CLI prints`);
    const printed = run.streams.get(11) || Buffer.alloc(0);
    assert.equal(sha256(printed), sha256(expected.stdout), `${label}: stream 11 is the CLI's standard output`);
    assert.equal((run.streams.get(3) || Buffer.alloc(0)).toString(), expected.stderr.toString(), `${label}: stream 3 is the CLI's standard error`);
    for (const stream of [0, 1, 4, 5, 6, 7, 9]) assert.ok(!(run.streams.get(stream)?.length), `${label}: nothing on stream ${stream}`);
    assert.equal((run.streams.get(10) || Buffer.alloc(0)).length, 4, `${label}: one signal of where the messages go (stream 10)`);
    return { json, run };
  };
  let count = 0;
  for (const program of Object.keys(FORMATS)) {
    for (const [kind, words] of [["help", HELP_WORDS], ["version", VERSION_WORDS]]) {
      for (const word of words) {
        prints([program, ...inputs, word, "-num_threads", "1"], kind);
        count++;
      }
    }
  }
  // After a help word the CLI reads no value; a version word comes before -dryrun, the registry file (none is read) and
  // the file of -logfile, which the CLI makes empty.
  prints(["blastn", ...inputs, "-h", "-evalue", "abc", "-num_threads", "1"], "help");
  // `-out`, which a run refuses, does not stop the help (line WQA's D14 review, 2026-11-20).
  prints(["blastp", ...inputs, "-h", "-out", "x", "-num_threads", "1"], "help");
  prints(["blastn", ...inputs, "-out", "x", "-xmlhelp", "-num_threads", "1"], "help");
  prints(["blastx", ...inputs, "-dryrun", "-help", "-num_threads", "1"], "help");
  assert.equal(prints(["blastp", ...inputs, "-version", "-conffile", "missing.ini", "-num_threads", "1"], "version").json.conffile, null);
  const logged = prints(["tblastx", ...inputs, "-version", "-logfile", "v.log", "-num_threads", "1"], "version");
  assert.equal(logged.json.log_file, "v.log");
  const made = fs.readFileSync(path.join(dir, "v.log"));
  assert.deepEqual(storedZipEntries(logged.run.streams.get(8)).map(([name, bytes]) => [name, sha256(bytes)]), [["v.log", sha256(made)]],
                   "the file of -logfile (stream 8) is the CLI's");
  fs.rmSync(path.join(dir, "v.log"));
  // A -logfile name with a directory: the CLI opens the file of -logfile before it prints its help and its version, in
  // a working directory where the name's directories exist (`checkLogfileNames`), and the run gives that file, empty
  // (line WQA's D14, 2026-11-21).
  for (const [word, kind] of [["-h", "help"], ["-version", "version"]]) {
    fs.mkdirSync(path.join(dir, "a"));
    const run = prints(["blastp", ...inputs, word, "-logfile", "a/b.log", "-num_threads", "1"], kind).run;
    const file = fs.readFileSync(path.join(dir, "a", "b.log"));
    assert.deepEqual(storedZipEntries(run.streams.get(8)).map(([name, bytes]) => [name, sha256(bytes)]), [["a/b.log", sha256(file)]],
                     `${word} -logfile a/b.log: the file of -logfile (stream 8) is the CLI's`);
    fs.rmSync(path.join(dir, "a"), { recursive: true });
  }
  // A help or version word as the value of an option that takes one is that value, as the CLI's parser reads it; only
  // NCBI's scan of the versions before the arguments reads every word (`-subject -version`).
  prints(["blastp", "-query", "q.fa", "-subject", "-version", "-num_threads", "1"], "version");
  prints(["tblastn", ...inputs, "-outfmt", "-version-full", "-num_threads", "1"], "version");
  count += 10;
  // What stops the CLI before its help: a value before the help word, and the registry file of -conffile, which the
  // command line reads before it reads its arguments.
  for (const argv of [["blastn", ...inputs, "-evalue", "abc", "-h", "-num_threads", "1"],
                      ["blastp", ...inputs, "-h", "-conffile", "missing.ini", "-num_threads", "1"]]) {
    const label = argv.join(" ");
    start(argv);
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    const expected = cli(argv);
    assert.notEqual(expected.status, 0, `${label}: the CLI fails`);
    assert.equal(validated.status, -1, `${label}: validate fails`);
    assert.equal(validated.error.trimEnd(), expected.stderr.toString().trimEnd(), `${label}: the CLI's message`);
    count++;
  }
  // The help words that are an option's value (line WQA's D14, 2026-11-21): the run is the CLI's, which fails, given the
  // registered inputs: `-query -h` names no subject (NCBI's error before the query is read), `-outfmt -h` and
  // `-evalue -h` are values that the CLI rejects.
  const sequences = { blastp: ["MKVLAAGIVG", "MKVLAAGIVG"], blastx: ["ACGTACGTAC", "MKVLAAGIVG"], tblastn: ["MKVLAAGIVG", "ACGTACGTAC"],
                      tblastx: ["ACGTACGTAC", "ACGTACGTAC"], blastn: ["ACGTACGTAC", "ACGTACGTAC"] };
  const fails = (argv) => {
    const label = argv.join(" ");
    assert.equal(start(argv).prints, null, `${label}: not a run that prints`);
    const expected = cli(argv);
    assert.notEqual(expected.status, 0, `${label}: the CLI fails`);
    const message = expected.stderr.toString().trimEnd();
    const validated = reactor.call("losat_web2_validate", argv.join("\0"));
    assert.equal(validated.status, -1, `${label}: validate fails`);
    assert.equal(validated.error.trimEnd(), message, `${label}: validate gives the CLI's message`);
    const [query, subject] = sequences[argv[0]].map((residues, kind) =>
      reactor.call("losat_web2_register", argv[0], kind, 1, 0, Buffer.from(`>r${kind}\n${residues}\n`)).status);
    assert.ok(query > 0 && subject > 0, `${label}: the inputs are registered`);
    const run = reactor.call("losat_web2_run", argv.join("\0"), query, subject);
    for (const handle of [query, subject]) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
    assert.equal(run.status, -1, `${label}: run fails`);
    assert.equal(run.error.trimEnd(), message, `${label}: run gives the CLI's message`);
    assert.ok(!run.streams.get(3)?.length && !run.streams.get(11)?.length, `${label}: nothing on streams 3 and 11`);
    count++;
  };
  for (const program of Object.keys(FORMATS)) fails([program, "-query", "-h", "-num_threads", "1"]);
  fails(["blastp", ...inputs, "-outfmt", "-h", "-num_threads", "1"]);
  fails(["blastn", ...inputs, "-evalue", "-xmlhelp", "-num_threads", "1"]);
  // `-out -h` is the file `-h` of -out, which a run refuses whatever its name (abi_v2.md §4): not the help.
  const out = ["blastp", ...inputs, "-out", "-h", "-num_threads", "1"];
  assert.equal(start(out).prints, null, "-out -h: not a run that prints");
  assert.equal(reactor.call("losat_web2_validate", out.join("\0")).error.trimEnd(), "-out is not accepted: a run writes every output format of the program");
  count++;
  fs.rmSync(dir, { recursive: true });
  return count;
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
  for (const [kind, name, error] of [[3, "x", /^unknown input file kind 3$/], [1, "", /^an input file needs a name$/]]) {
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

// -logfile names (abi_v2.md §4 *start*; line WQA's sweep, D11, 2026-11-20): a name is read as the CLI on Linux reads it
// in a working directory where the name's directories exist. BLASTN and TBLASTX with -max_target_seqs 1 (TBLASTX
// writes the warning where -logfile sends the diagnostics) on the outfmt11 inputs, without and with a `.ncbirc` that
// NCBI cannot read (the one byte 0xEF: the usage report's error, which the CLI writes to its standard error at once and
// to a -logfile file when it opens it), each against the CLI run in a new directory that has the name's directories, one
// run per format with -out: outfmt 0, 6 and 7; stream 3 is the CLI's standard error of outfmt 6; stream 8 holds the file
// that the name makes, under the name, or, for the names of the standard output, what the CLI writes there (its whole
// standard output, the reports being in -out); stream 10 says where the warning went. `<dir>` is the run's directory.
const LOG_NAMES = ["run.log", "logs/run.log", "a/b/run.log", "./run.log", "<dir>/abs/run.log", "/dev/stdout", "/dev/fd/1",
  "-", "/dev/stderr", "/dev/fd/2", "/proc/self/fd/2", "//dev/stderr", ".", "..", "logs/", "logs/.", "", "/dev/null", "//dev/null",
  // Line WQA's D11 round 2 (2026-11-20): names that NCBI reads otherwise than the system: a trailing '/' names the file
  // without it, the spaces at the ends of a relative name are cut, and a file system that reports no free space
  // (procfs: /dev/fd, /proc) gets no message; and devices whose writes are lost.
  "run.log/", "logs//", "a/b/", "<dir>/abs/run.log/", " run.log", "run.log ", " logs/run.log", "  ", "/proc/thread-self/fd/2",
  "/dev/fd/../stderr", "/dev/full", "/dev/zero",
  // What the system does not open (a file descriptor that the run has not opened, a file that the user may not write, a
  // name in a device), and what it opens otherwise than by the name: the standard input (the null device), a file of
  // procfs (no free space) and a file of the working directory through procfs.
  "/dev/fd/3", "/proc/self/status", "/dev/null/x", "/dev/stdout/../stderr", "/dev/stdin", "/dev/fd/0", "/proc/self/comm",
  "/proc/self/cwd/run.log",
  // Round 3: a directory of '/'s only is none (`CDirEntry::Reset` empties `//`), a user makes no file in a directory of
  // root's and opens none of its devices but those whose writes are lost (a new terminal's master).
  "//run.log", "/run.log", "/etc/run.log", "/var/log/run.log", "/home/run.log", "/dev/run.log", "/dev/kmsg", "/dev/console",
  "/dev/ptmx",
  // Round 4: one rule of owners (a user makes files in the working directory, a home and the directories where every user
  // does; every other directory is root's, where only root makes one), the links of Linux (`/var/lock` is `/run/lock`,
  // `/bin` is `/usr/bin`), and the system's limits of a path (a component of 255 bytes). `<id>` is a name of this run's
  // own, for the files that the CLI makes outside the run's directory.
  "/var/lock/<id>.log", "/var/run/../tmp/<id>.log", "/data/run.log", "/lost+found/run.log", "/snap/run.log",
  "/bin/../tmp/<id>.log", "/var/lock/../tmp/<id>.log", "/dev/fd/../../tmp/<id>.log", "a".repeat(256),
  `logs/${"b".repeat(256)}`,
  // Round 5: a name climbs out of directories of root's that the model does not list, as the CLI passes them.
  "/etc/default/../../tmp/<id>.log", "/usr/share/doc/../../../tmp/<id>.log",
  // Round 6: a name looks nothing up in a directory that a user cannot search (mode 0700, root's), as the CLI's EACCES.
  "/root/../tmp/<id>.log", "/var/cache/ldconfig/../../../tmp/<id>.log",
  // Round 7: the system follows 40 links in a lookup (`/var/lock` is one) and refuses the 41st, as the CLI's ELOOP.
  `${"/var/lock/../..".repeat(40)}/tmp/<id>.log`, `${"/var/lock/../..".repeat(41)}/tmp/<id>.log`];

// The rows of `/`-names (LOG_NAMES) compare the Web with the CLI of the host, which must be a Linux like LOSAT Web's
// (Ubuntu 24.04, as CI's runner): a user that is not root, `/var/lock` and `/var/run` links to `/run/lock` and `/run`,
// `/run/lock` and `/tmp` where every user makes files (mode 1777), `/etc/default` and `/usr/share/doc`, and `/root` and
// `/var/cache/ldconfig` that the user cannot search. Returns why this host differs, or undefined.
function hostDiffers() {
  const sticky = (directory) => {
    try {
      return (fs.statSync(directory).mode & 0o1777) === 0o1777;
    } catch {
      return false;
    }
  };
  const resolves = (link, target) => {
    try {
      return fs.realpathSync(link) === target;
    } catch {
      return false;
    }
  };
  const isDirectory = (directory) => {
    try {
      return fs.statSync(directory).isDirectory();
    } catch {
      return false;
    }
  };
  const closed = (directory) => {
    if (!isDirectory(directory)) return false;
    try {
      fs.accessSync(directory, fs.constants.X_OK);
      return false;
    } catch {
      return true;
    }
  };
  const differences = [
    [process.getuid?.() === 0, "the run is root's"],
    [!resolves("/var/lock", "/run/lock"), "/var/lock is not /run/lock"],
    [!resolves("/var/run", "/run"), "/var/run is not /run"],
    [!sticky("/run/lock"), "/run/lock is not mode 1777"],
    [!sticky("/tmp"), "/tmp is not mode 1777"],
    [!isDirectory("/etc/default"), "no /etc/default"],
    [!isDirectory("/usr/share/doc"), "no /usr/share/doc"],
    [!closed("/root"), "the user searches /root, or it is not there"],
    [!closed("/var/cache/ldconfig"), "the user searches /var/cache/ldconfig, or it is not there"],
  ].filter(([differs]) => differs).map(([, why]) => why);
  return differences.length > 0 ? differences.join(", ") : undefined;
}

// The CLI's run with its standard output and error on pipes, as a shell's `|` or a terminal gives them: node's child
// processes get sockets, which `/dev/stdout` and `/dev/stderr` cannot open (ENXIO), so that NCBI and the CLI send the
// log to `cerr` (python's subprocess makes pipes). Its standard input is the null device, as LOSAT Web's run has none.
function nativeOnPipes(native, argv, directory) {
  const outFile = path.join(directory, ".stdout"), errFile = path.join(directory, ".stderr");
  const script = "import subprocess, sys\n" +
    "p = subprocess.run(sys.argv[3:], capture_output=True, stdin=subprocess.DEVNULL)\n" +
    "open(sys.argv[1], 'wb').write(p.stdout); open(sys.argv[2], 'wb').write(p.stderr); sys.exit(p.returncode)\n";
  const result = require("node:child_process").spawnSync("python3", ["-I", "-c", script, outFile, errFile, native, ...argv],
                                                         { cwd: directory, env: { PATH: process.env.PATH, HOME: directory } });
  assert.equal(result.stderr.length, 0, `python3: ${result.stderr}`);
  const stdout = fs.readFileSync(outFile), stderr = fs.readFileSync(errFile);
  fs.rmSync(outFile);
  fs.rmSync(errFile);
  return { status: result.status, stdout, stderr };
}

// The regular files under `directory`, as paths.
function filesIn(directory) {
  return fs.readdirSync(directory, { recursive: true, withFileTypes: true })
    .filter((entry) => entry.isFile())
    .map((entry) => path.join(entry.parentPath ?? entry.path, entry.name));
}

async function checkLogfileNames(reactor, native, scratch, engine, waitForWorkers) {
  const differs = hostDiffers();
  assert.ok(differs === undefined, `-logfile names: this host is not a Linux like LOSAT Web's (${differs}); ` +
            "run V-ABI as a user on Ubuntu 24.04 or a Linux like it");
  const inputs = path.join(ROOT, "LOSAT/tests/fasta/outfmt11");
  const query = path.join(inputs, "arch_tbx_q.fna"), subject = path.join(inputs, "arch_tbx_s.fna");
  const standardOutput = new Set(["/dev/stdout"]);
  const device = (name) => name === "-" || name === "" || /^\/+(dev|proc|sys)\//.test(name);
  let rows = 0;
  for (const program of ["blastn", "tblastx"]) {
    const handles = [0, 1].map((role) => {
      const registered = reactor.call("losat_web2_register", program, role, 1, 0, fs.readFileSync(role === 0 ? query : subject));
      assert.ok(registered.status > 0, registered.error);
      return registered.status;
    });
    for (const ncbirc of [false, true]) {
      let registry;
      if (ncbirc) {
        const begun = reactor.call("losat_web2_file_begin", 2, ".ncbirc", 1);
        assert.ok(begun.status > 0, begun.error);
        assert.equal(reactor.call("losat_web2_file_chunk", begun.status, Buffer.from([0xef])).status, 0);
        assert.equal(reactor.call("losat_web2_file_end", begun.status).status, 0);
        registry = begun.status;
      }
      for (const template of LOG_NAMES) {
        const directory = fs.mkdtempSync(path.join(scratch, "logfile."));
        const name = template.replace("<dir>", directory).replace("<id>", `losat-v-abi-${process.pid}-${engine}-${rows}`);
        // The file that the CLI makes outside the run's directory: the path that NCBI opens, the name read lexically.
        const outside = template.includes("<id>") ? path.posix.normalize(name) : undefined;
        const label = `${engine} ${program}${ncbirc ? " .ncbirc 0xEF" : ""} -logfile ${JSON.stringify(name)}`;
        // The directories that the name gives: the one that NCBI checks (the name's, without its trailing '/'s) and that
        // of the file it opens (the name without the spaces at its ends).
        // Not outside it: a user makes no directory of root's (`/data`), which the CLI then does not find.
        if (!device(name) && (!name.trim().startsWith("/") || name.startsWith(directory))) {
          for (const given of [name.replace(/\/+$/, ""), name.trim().replace(/\/+$/, "")]) {
            const slash = given.lastIndexOf("/");
            if (slash > 0) fs.mkdirSync(path.resolve(directory, given.slice(0, slash)), { recursive: true });
          }
        }
        if (ncbirc) fs.writeFileSync(path.join(directory, ".ncbirc"), Buffer.from([0xef]));
        const before = new Set(filesIn(directory));
        const words = [program, "-query", query, "-subject", subject, "-max_target_seqs", "1"];
        const outputs = new Map();
        let cli;
        for (const format of FORMATS[program]) {
          const out = path.join(directory, `native.${format}.out`);
          const result = nativeOnPipes(native, [...words, "-outfmt", String(format), "-num_threads", "1", "-out", out, "-logfile", name], directory);
          assert.equal(result.status, 0, `${label}: native outfmt ${format}: ${result.stderr}`);
          outputs.set(format, fs.readFileSync(out));
          fs.rmSync(out);
          // The file that the run made (NCBI's name for it: `logs` for `logs/`, `run.log` for ` run.log`).
          const madeFiles = filesIn(directory).filter((file) => !before.has(file));
          if (outside !== undefined && fs.existsSync(outside)) madeFiles.push(outside);
          assert.ok(madeFiles.length <= 1, `${label}: native outfmt ${format} made ${madeFiles}`);
          const made = madeFiles.length === 1 ? fs.readFileSync(madeFiles[0]) : undefined;
          // The next format's run appends to the file: each run starts without it.
          if (made !== undefined) fs.rmSync(madeFiles[0]);
          if (format === 6) cli = { stderr: result.stderr, stdout: result.stdout, made };
        }
        const files = standardOutput.has(name) ? [[name, cli.stdout]] : cli.made === undefined ? [] : [[name, cli.made]];
        if (!standardOutput.has(name)) assert.equal(cli.stdout.length, 0, `${label}: the CLI's standard output with -out`);
        const run = reactor.call("losat_web2_run", [...words, "-logfile", name, "-num_threads", "1"].join("\0"), ...handles);
        if (waitForWorkers) await waitForWorkers();
        assert.equal(run.status, 0, `${label}: ${run.error}`);
        for (const format of FORMATS[program]) {
          assert.equal(sha256(run.streams.get(format) || Buffer.alloc(0)), sha256(outputs.get(format)), `${label}: outfmt ${format}`);
        }
        assert.equal((run.streams.get(3) || Buffer.alloc(0)).toString(), cli.stderr.toString(), `${label}: stream 3 is the CLI's standard error`);
        const sent = run.streams.has(8) ? storedZipEntries(run.streams.get(8)).map(([entry, bytes]) => [entry, bytes.toString()]) : [];
        assert.deepEqual(sent, files.map(([entry, bytes]) => [entry, bytes.toString()]), `${label}: the files of the run (stream 8)`);
        const warning = `Warning: [${program}] Examining 5 or more matches is recommended\n`;
        const place = run.streams.get(10).readUInt32LE(0);
        const where = cli.stderr.toString().includes(warning) ? 3 : files.some(([, bytes]) => bytes.toString().includes(warning)) ? 8 : 0xffffffff;
        assert.equal(place, where, `${label}: where the run's messages go (stream 10)`);
        fs.rmSync(directory, { recursive: true, force: true });
        rows++;
      }
      if (registry !== undefined) assert.equal(reactor.call("losat_web2_file_release", registry).status, 0);
    }
    for (const handle of handles) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
  }
  return rows;
}

// -export_search_strategy (line RG item F1w, review 1 S3): the file of a run of two or more query batches at NCBI's
// default batch sizes (the browser gives no BATCH_SIZE), for each program, on each engine and its thread counts (the
// threaded engine's parallel search, -num_threads 4 too), against the CLI's file of the same argv: stream 8 holds the
// CLI's file under the name, and the argv's output (outfmt 6, and outfmt 11's archive, whose options handle the
// strategy shares) is the CLI's. The queries and subjects are seeded letters, written once into `scratch`; each
// subject is pieces of the queries, so that the searches find hits.
function seededFasta(seed, records, length, letters, prefix) {
  let state = seed >>> 0;
  const next = () => {
    state = (state + 0x6d2b79f5) >>> 0;
    let t = state;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
  const sequences = [];
  for (let record = 0; record < records; record++) {
    let sequence = "";
    for (let at = 0; at < length; at++) sequence += letters[Math.floor(next() * letters.length)];
    sequences.push(sequence);
  }
  const text = sequences.map((sequence, record) => `>${prefix}${record} record ${record}\n${sequence.match(/.{1,60}/g).join("\n")}\n`).join("");
  return { text, sequences };
}

const STRATEGY_BATCHES = [
  // [program, words, query: [records, length, alphabet], subject alphabet]: query letters past the program's batch size
  // (blastp 10000, tblastn 20000 residues, blastx and tblastx 10002, BLASTN's -task blastn 100000 bases).
  ["blastp", [], [30, 400, "protein"], "protein"],
  ["tblastn", [], [60, 400, "protein"], "nucleotide"],
  ["blastx", [], [30, 400, "nucleotide"], "protein"],
  ["tblastx", [], [30, 400, "nucleotide"], "nucleotide"],
  ["blastn", ["-task", "blastn"], [260, 400, "nucleotide"], "nucleotide"],
];
const ALPHABETS = { protein: "ACDEFGHIKLMNPQRSTVWY", nucleotide: "ACGT" };
// The CLI's translation of a piece of a nucleotide query, for a protein subject (the standard code).
function translated(bases) {
  const code = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
  const index = { T: 0, C: 1, A: 2, G: 3 };
  let protein = "";
  for (let at = 0; at + 3 <= bases.length; at += 3) {
    protein += code[index[bases[at]] * 16 + index[bases[at + 1]] * 4 + index[bases[at + 2]]];
  }
  return protein.replaceAll("*", "W");
}

async function checkStrategyFiles(reactor, native, scratch, engine, waitForWorkers) {
  let rows = 0;
  for (const [index, [program, words, [records, length, queryAlphabet], subjectAlphabet]] of STRATEGY_BATCHES.entries()) {
    const query = seededFasta(0xf1 + index, records, length, ALPHABETS[queryAlphabet], "q");
    // Three subjects, each of pieces of three queries: a protein subject of a nucleotide query takes their translation,
    // a nucleotide subject of a protein query takes seeded bases.
    const pieces = (from) => [0, 1, 2].map((shift) => {
      const sequence = query.sequences[(from + shift * 7) % records].slice(60, 360);
      if (subjectAlphabet === queryAlphabet) return sequence;
      return subjectAlphabet === "protein" ? translated(sequence) : seededFasta(0x5e + from + shift, 1, 900, ALPHABETS.nucleotide, "x").sequences[0];
    }).join("");
    const subjectText = [0, 1, 2].map((record) => `>s${record} subject ${record}\n${pieces(record * 3).match(/.{1,60}/g).join("\n")}\n`).join("");
    const directory = fs.mkdtempSync(path.join(scratch, "strategy."));
    fs.writeFileSync(path.join(directory, "q.fa"), query.text);
    fs.writeFileSync(path.join(directory, "s.fa"), subjectText);
    const handles = [query.text, subjectText].map((text, role) => {
      const registered = reactor.call("losat_web2_register", program, role, 1, 0, Buffer.from(text));
      assert.ok(registered.status > 0, registered.error);
      return registered.status;
    });
    for (const format of program === "blastp" || program === "blastx" ? [6, 11] : [6]) {
      const argv = [program, "-query", "q.fa", "-subject", "s.fa", ...words, "-outfmt", String(format), "-export_search_strategy", "s.asn"];
      const cli = nativeOnPipes(native, [...argv, "-num_threads", "1", "-out", "o.out"], directory);
      const label = `${engine.name} ${program} outfmt ${format}`;
      assert.equal(cli.status, 0, `${label}: native: ${cli.stderr}`);
      const file = fs.readFileSync(path.join(directory, "s.asn")), out = fs.readFileSync(path.join(directory, "o.out"));
      const requests = file.toString().split("Blast4-request ::=").length - 1;
      assert.ok(requests >= 2, `${label}: the CLI wrote ${requests} request(s), not a run of two or more batches`);
      for (const threads of engine.threads) {
        const run = reactor.call("losat_web2_run", [...argv, "-num_threads", String(threads)].join("\0"), ...handles);
        if (waitForWorkers) await waitForWorkers();
        assert.equal(run.status, 0, `${label} threads ${threads}: ${run.error}`);
        const sent = run.streams.has(8) ? storedZipEntries(run.streams.get(8)) : [];
        assert.deepEqual(sent.map(([name, bytes]) => [name, sha256(bytes)]), [["s.asn", sha256(file)]],
                         `${label} threads ${threads}: the file of the run (stream 8) is the CLI's (${requests} requests)`);
        const stream = format === 6 ? 6 : CUSTOM_STREAM;
        assert.equal(sha256(run.streams.get(stream) || Buffer.alloc(0)), sha256(out), `${label} threads ${threads}: outfmt ${format}`);
        rows++;
      }
      fs.rmSync(path.join(directory, "s.asn"));
      fs.rmSync(path.join(directory, "o.out"));
    }
    for (const handle of handles) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
    fs.rmSync(directory, { recursive: true, force: true });
  }
  return rows;
}

// The frames of stream 12 (abi_v2.md §5): `chunks`, the run's emits, as the files they write in the order of their
// numbers: [{ name, bytes, writes, ended }], `bytes` the file as the frames leave it and `writes` the number of its
// frames of bytes. Each frame is one emit of at most 1 MiB: 8 bytes (the file's number from 0, in the order the run
// opens its files, and the kind, each a u32 LE) and a payload; a file's frames come as beginning (0, its name), bytes
// (1, at a u64 LE offset, written there with zero bytes before them past the end), truncation (3, a u64 LE length) and
// end (2), and every file that the run opened ends when the run does, in the order of their numbers.
function toolFiles(chunks) {
  const files = [];
  let ends = 0;
  for (const chunk of chunks) {
    assert.ok(chunk.length >= 8 && chunk.length <= 1 << 20, `a frame of ${chunk.length} bytes`);
    const [number, kind] = [chunk.readUInt32LE(0), chunk.readUInt32LE(4)];
    const payload = chunk.subarray(8);
    if (kind === 0) {
      assert.equal(number, files.length, "files are numbered from 0 in the order they open");
      assert.equal(ends, 0, "no file opens after the first end (the ends come when the run ends)");
      files.push({ name: payload.toString("utf8"), bytes: Buffer.alloc(0), writes: 0, ended: false });
      continue;
    }
    assert.ok(number < files.length && !files[number].ended, `frame of kind ${kind} of file ${number}: it is open`);
    const file = files[number];
    if (kind === 1 || kind === 3) {
      assert.ok(payload.length >= 8 && (kind === 1 || payload.length === 8), `frame of kind ${kind} of ${payload.length} bytes`);
      const at = Number(payload.readBigUInt64LE(0));
      assert.ok(Number.isSafeInteger(at) && at < 1 << 30, `an offset or length of ${payload.readBigUInt64LE(0)}`);
      if (kind === 1) {
        const bytes = payload.subarray(8);
        assert.ok(bytes.length > 0, "a frame of bytes has bytes");
        const end = at + bytes.length;
        if (end > file.bytes.length) file.bytes = Buffer.concat([file.bytes, Buffer.alloc(end - file.bytes.length)]);
        bytes.copy(file.bytes, at);
        file.writes++;
      } else {
        file.bytes = at <= file.bytes.length ? file.bytes.subarray(0, at) : Buffer.concat([file.bytes, Buffer.alloc(at - file.bytes.length)]);
      }
    } else {
      assert.equal(kind, 2, `the kind of a frame (${kind})`);
      assert.equal(payload.length, 0, "the end of a file has no payload");
      assert.equal(number, ends, "the files end in the order of their numbers");
      file.ended = true;
      ends++;
    }
  }
  for (const [number, file] of files.entries()) assert.ok(file.ended, `file ${number} (${file.name}) ends with the run`);
  return files;
}

// The folder of a tool search's files: `files` (published path -> source file) copied under `directory`.
function writeToolFiles(directory, files) {
  for (const [name, source] of Object.entries(files)) {
    fs.mkdirSync(path.dirname(path.join(directory, name)), { recursive: true });
    fs.copyFileSync(source, path.join(directory, name));
  }
}

// One search of `"tool"` (see the header): `instance` is the engine's instance whose "/work" is `directory`, a folder of
// the search's files, which the engine does not write to. Returns the results row.
async function checkTool(engine, instance, directory, search, options, scratch) {
  const label = `${engine.name} ${search.cases.join(",")}`;
  const words = ["windowmasker", ...search.argv];
  const before = new Set(filesIn(directory).map((file) => path.relative(directory, file)));
  // The CLI's run, in a folder of the same files.
  const folder = fs.mkdtempSync(path.join(scratch, "tool."));
  writeToolFiles(folder, search.files);
  const known = new Set(filesIn(folder).map((file) => path.relative(folder, file)));
  const cli = nativeOnPipes(options.native, words, folder);
  // The CLI writes the folder's path where the run's directory is named (NCBI's search path of a database that is not
  // found, `[db/none] in search path [<folder>::]`); the engine's directory is `/work` (the host's preopen, abi_v2.md §3).
  const place = (bytes) => Buffer.from(bytes.toString("latin1").split(folder).join("/work"), "latin1");
  cli.stdout = place(cli.stdout);
  cli.stderr = place(cli.stderr);
  const made = filesIn(folder).map((file) => path.relative(folder, file)).filter((name) => !known.has(name)).sort();
  const created = made.map((name) => [name, place(fs.readFileSync(path.join(folder, name)))]);
  fs.rmSync(folder, { recursive: true, force: true });
  const reactor = instance.reactor;
  let row;
  if (search.tool === "validate") {
    const validated = reactor.call("losat_web2_validate", words.join("\0"));
    assert.equal(reactor.chunks.size, 0, `${label}: validate emits nothing`);
    if (search.expect === "late_failure") {
      // The command line's run fails after it opened a file, where validation stops (abi_v2.md §4).
      assert.notEqual(cli.status, 0, `${label}: the CLI's run fails`);
      assert.equal(validated.status, 0, `${label}: validate accepts what the command line checks before its files: ${validated.error}`);
      assert.equal(validated.streams.size, 0, `${label}: nothing emitted`);
    } else if (search.expect !== "failure") {
      assert.equal(cli.status, 0, `${label}: the CLI's run of a valid argv (${cli.stderr})`);
      assert.equal(validated.status, 0, `${label}: validate: ${validated.error}`);
      assert.equal(validated.streams.size, 0, `${label}: nothing emitted`);
    } else {
      // The CLI's standard error, then the file of -logfile, without the last line end.
      const log = created.filter(([name]) => name === search.argv[search.argv.indexOf("-logfile") + 1]).map(([, bytes]) => bytes);
      let text = Buffer.concat([cli.stderr, ...log]).toString();
      if (text.endsWith("\n")) text = text.slice(0, -1);
      assert.notEqual(cli.status, 0, `${label}: the CLI's run of an invalid argv`);
      assert.equal(validated.status, -1, `${label}: validate refuses the argv (the CLI's exit ${cli.status})`);
      assert.equal(validated.error, text, `${label}: validate gives the CLI's standard error`);
      assert.equal(validated.streams.size, 0, `${label}: nothing emitted`);
    }
    row = { engine: engine.name, threads: 0, program: search.program, cases: search.cases, argv: search.argv, formats: {}, frozen: {},
            tool: "validate", status: validated.status };
  } else {
    const run = reactor.call("losat_web2_tool", words.join("\0"));
    const frames = (reactor.chunks.get(12) || []).map((chunk) => Buffer.from(chunk));
    if (instance.waitForWorkers) await instance.waitForWorkers();
    assert.equal(run.status, 0, `${label}: ${run.error}`);
    const status = JSON.parse(run.streams.get(2).toString()).status;
    if (search.expect === "refusal") {
      // LOSAT's refusal of an output name that LOSAT Web does not write (abi_v2.md §4), where the command line writes
      // the file: exit status 1, the refusal on stream 3, no output and no file.
      const named = search.argv[search.argv.indexOf("-out") + 1];
      assert.equal(status, 1, `${label}: the refusal's exit status`);
      assert.match((run.streams.get(3) || "").toString(), /^Error: an output file '.*' .* is not supported by LOSAT's WINDOWMASKER\n$/, `${label}: the refusal`);
      assert.ok((run.streams.get(3) || "").toString().includes(`'${named}'`), `${label}: the refusal names ${named}`);
      assert.ok(!run.streams.has(11) && toolFiles(frames).length === 0, `${label}: no output and no file`);
      assert.deepEqual(filesIn(directory).map((file) => path.relative(directory, file)).sort(), [...before].sort(), `${label}: the run's folder`);
      console.log(`ok ${label}`);
      return { engine: engine.name, threads: 0, program: search.program, cases: search.cases, argv: search.argv, formats: {}, frozen: {},
               tool: "run", status, files: [] };
    }
    assert.equal(status, cli.status, `${label}: the exit status (${run.streams.get(3) || ""})`);
    assert.equal(status !== 0, search.expect === "failure", `${label}: the run ${status !== 0 ? "fails" : "completes"} (exit ${status})`);
    assert.equal(sha256(run.streams.get(11) || Buffer.alloc(0)), sha256(cli.stdout), `${label}: stream 11 is the CLI's standard output`);
    assert.equal((run.streams.get(3) || Buffer.alloc(0)).toString(), cli.stderr.toString(), `${label}: stream 3 is the CLI's standard error`);
    for (const stream of [0, 1, 4, 5, 6, 7, 8, 9, 10]) assert.ok(!run.streams.has(stream), `${label}: no stream ${stream}`);
    const files = toolFiles(frames);
    assert.deepEqual(files.map((file) => file.name).sort(), created.map(([name]) => name), `${label}: the files of the run`);
    for (const file of files) {
      const want = created.find(([name]) => name === file.name)[1];
      assert.equal(sha256(file.bytes), sha256(want), `${label}: the bytes of ${file.name} (${want.length} bytes)`);
    }
    const logfile = search.argv.indexOf("-logfile");
    if (logfile >= 0 && files.length > 0) assert.equal(files[0].name, search.argv[logfile + 1], `${label}: the file of -logfile opens first`);
    const largest = Math.max(0, ...files.map((file) => file.writes));
    if (search.frames !== undefined) assert.ok(largest >= search.frames, `${label}: ${largest} frames of bytes, at least ${search.frames}`);
    row = { engine: engine.name, threads: 0, program: search.program, cases: search.cases, argv: search.argv, formats: {}, frozen: {},
            tool: "run", status, files: files.map((file) => [file.name, file.bytes.length, file.writes]) };
  }
  // The engine does not write to the run's directory.
  assert.deepEqual(filesIn(directory).map((file) => path.relative(directory, file)).sort(), [...before].sort(), `${label}: the run's folder`);
  console.log(`ok ${label}`);
  return row;
}

async function main() {
  const options = parseArgs(process.argv.slice(2));
  const scratch = fs.mkdtempSync(path.join(os.tmpdir(), "losat-v-abi-"));
  const results = [];
  const selected = (options.engines || "serial,threads").split(",");
  const engines = [
    { name: "serial", open: (directory) => openSerial(options.serial, directory), threads: [1] },
    { name: "threads", open: (directory) => openThreads(options.threads, directory), threads: [1, 2, 4] },
  ].filter((engine) => selected.includes(engine.name));
  const all = JSON.parse(fs.readFileSync(options.cases, "utf8"));
  for (const engine of engines) {
    const { reactor: own, close, waitForWorkers } = await engine.open();
    assert.equal(own.exports.losat_web2_abi_version(), 14);
    checkSurface(own, options.native);
    console.log(`ok ${engine.name} describe, validate, register, scan`);
    const logRows = await checkLogfileNames(own, options.native, scratch, engine.name, waitForWorkers);
    console.log(`ok ${engine.name} -logfile names: ${logRows} runs as the CLI's`);
    const strategyRows = await checkStrategyFiles(own, options.native, scratch, engine, waitForWorkers);
    console.log(`ok ${engine.name} -export_search_strategy: ${strategyRows} runs of two or more batches as the CLI's`);
    console.log(`ok ${engine.name} help and version: ${checkPrints(own, options.native)} argvs as the CLI`);
    const subjects = new Map(); // "program subject data_loaders" -> handle, kept across runs
    const subjectReactors = new Map(); // the same key -> the reactor that holds the handle (a run directory's instance has its own)
    // Per reactor, the -window_masker_db value or taxid file name -> { handle, file }, kept across runs: a file registered
    // on one instance (the engine's own or a run directory's) is known there only, and released there.
    const counts = new Map();
    const countsOn = (reactor) => counts.get(reactor) ?? counts.set(reactor, new Map()).get(reactor);
    let countsRegistered = 0, countsKept = 0;
    const databases = new Map(); // a -db search's directory -> the engine's instance whose "/work" it is
    const tools = new Map(); // a tool search's files -> [folder, the engine's instance whose "/work" it is]
    for (const search of all) {
      const program = search.program;
      if (search.tool) {
        const key = JSON.stringify(search.files);
        if (!tools.has(key)) {
          const directory = fs.mkdtempSync(path.join(scratch, "toolfiles."));
          writeToolFiles(directory, search.files);
          tools.set(key, [directory, await engine.open(directory)]);
        }
        const [directory, instance] = tools.get(key);
        results.push(await checkTool(engine, instance, directory, search, options, scratch));
        continue;
      }
      // A search with `"run_directory"` has a subject like any other and files in its directory too (the repeat
      // database of BLASTN `-filtering_db`, line RG item D4): it runs on the instance of that directory.
      const inDirectory = search.database || search.run_directory;
      if (inDirectory && !databases.has(search.cwd)) databases.set(search.cwd, await engine.open(search.cwd));
      const instance = inDirectory ? databases.get(search.cwd) : { reactor: own, waitForWorkers };
      const reactor = instance.reactor;
      const query = inputPath(search, "-query"), subject = search.database ? undefined : inputPath(search, "-subject");
      const expected = nativeOutputs(options.native, search, scratch);
      const started = reactor.call("losat_web2_start", search.argv.join("\0"));
      const start = started.status === 0 ? JSON.parse(started.streams.get(2).toString()) : undefined;
      const together = start?.register_inputs === true;
      // The host registers the file of -conffile, then validates the argv, whose *validation* says whether the
      // run reads its inputs with NCBI's data loaders on, and registers the inputs so (abi_v2.md §4, §9; line
      // WQA's D16, 2026-11-20). An argv that validate refuses is read with them on, as before any registry.
      const registry = start?.conffile ? registerRegistry(reactor, search, start.conffile) : null;
      const custom = search.outfmt === undefined ? [] : ["-outfmt", search.outfmt];
      const first = reactor.call("losat_web2_validate", [...search.argv, ...custom, "-num_threads", "1"].join("\0"));
      const dataLoaders = first.status === 0 ? JSON.parse(first.streams.get(2).toString()).data_loaders : true;
      const subjectKey = `${program} ${subject} ${dataLoaders} ${inDirectory ? search.cwd : ""}`;
      // An input that the CLI's reader refuses (a first line that may be a Seq-id with the data loaders on) fails
      // the search as the CLI fails, before any hit: the subject first, as the CLI reads it first. The application
      // shows the error on the input and runs nothing; the CLI has written no output for a subject (read with the
      // arguments) and at most the outfmt 0 prolog for a query (read after it).
      const refusedInput = (registered, role) => {
        assert.equal(search.expect, "failure", `${engine.name} ${search.cases.join(",")}: register ${role}: ${registered.error}`);
        assert.ok(registered.error.length > 0 && expected.own.toString().endsWith(cliError(registered.error)),
                  `${engine.name} ${search.cases.join(",")}: the CLI's error (${registered.error})`);
        for (const format of FORMATS[program]) {
          if (role === "query" && format === 0) continue;
          assert.equal(expected.outputs.get(format).length, 0, `${search.cases.join(",")}: the CLI writes no outfmt ${format}`);
        }
        results.push({ engine: engine.name, threads: 0, program, cases: search.cases, argv: search.argv, formats: {}, frozen: {}, failure: registered.error });
        console.log(`ok ${engine.name} ${search.cases.join(",")} (fails as the CLI, reading the ${role})`);
        if (registry !== null) assert.equal(reactor.call("losat_web2_file_release", registry).status, 0);
      };
      if (!together && !search.database && !subjects.has(subjectKey)) {
        const registered = reactor.call("losat_web2_register", program, 1, dataLoaders ? 1 : 0, 0, fs.readFileSync(subject));
        if (registered.status < 0) {
          refusedInput(registered, "subject");
          continue;
        }
        subjects.set(subjectKey, registered.status);
        subjectReactors.set(subjectKey, reactor);
      }
      const q = together ? undefined : reactor.call("losat_web2_register", program, 0, dataLoaders ? 1 : 0, 0, fs.readFileSync(query));
      if (q && q.status < 0) {
        refusedInput(q, "query");
        continue;
      }
      let queryIds = q && JSON.parse(q.streams.get(2).toString()).records.map((record) => record.id);
      const taxid = registerTaxidCounts(reactor, countsOn(reactor), search);
      const named = registerCounts(reactor, countsOn(reactor), search) ?? taxid;
      if (named) named.held === "registered" ? countsRegistered++ : countsKept++;
      for (const threads of engine.threads) {
        const argv = [...search.argv, ...custom, "-num_threads", String(threads)].join("\0");
        // The application validates the argv on the instance that holds its counts file, then runs it
        // (abi_v2.md §4): an error of validate is the CLI's first error. So is that of a search that
        // fails without a counts file, whose validate makes the checks of the argv's own search (§7).
        if (named ? named.held !== "absent" : search.expect === "failure") {
          const validated = reactor.call("losat_web2_validate", argv);
          if (search.expect === "failure" && validated.status !== 0) {
            assert.ok(expected.own.toString().endsWith(cliError(validated.error)), `${engine.name} n${threads} ${search.cases.join(",")}: validate gives the CLI's error (${validated.error})`);
          } else {
            assert.equal(validated.status, 0, `${engine.name} n${threads} ${search.cases.join(",")}: validate: ${validated.error}`);
          }
        }
        const label = `${engine.name} n${threads} ${search.cases.join(",")}`;
        let handles = [q?.status, search.database ? 0 : subjects.get(subjectKey)], shownIds;
        if (together) {
          const inputs = reactor.call("losat_web2_register_inputs", argv, fs.readFileSync(query), subject === undefined ? Buffer.alloc(0) : fs.readFileSync(subject));
          if (inputs.status < 0) {
            // The CLI fails reading the inputs, before any output.
            assert.equal(search.expect, "failure", `${label}: register_inputs: ${inputs.error}`);
            assert.ok(inputs.error.length > 0 && expected.own.toString().endsWith(`${inputs.error}\n`), `${label}: the CLI's error (${inputs.error})`);
            for (const format of FORMATS[program]) assert.equal(expected.outputs.get(format).length, 0, `${label}: the CLI writes no outfmt ${format}`);
            results.push({ engine: engine.name, threads, program, cases: search.cases, argv: search.argv, formats: {}, frozen: {}, failure: inputs.error });
            console.log(`ok ${label} (fails as the CLI, reading the inputs)`);
            continue;
          }
          const response = JSON.parse(inputs.streams.get(2).toString());
          handles = [response.query.handle, response.subject === null ? 0 : response.subject.handle];
          queryIds = response.query.records.map((record) => record.id);
          shownIds = { query: response.query.records.map((record) => record.shown_id), subject: (response.subject?.records ?? []).map((record) => record.shown_id) };
        }
        const run = reactor.call("losat_web2_run", argv, ...handles);
        if (instance.waitForWorkers) await instance.waitForWorkers();
        if (together) for (const handle of handles.filter((handle) => handle !== 0)) assert.equal(reactor.call("losat_web2_release", handle).status, 0);
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
          if (search.frozen[format]) row.frozen[format] = ncbiHash(toFrozen(actual, search.frozen_prefix)) === search.frozen[format];
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
          const master = ncbiHash(toFrozen(entries[entries.length - 1][1], search.frozen_prefix));
          const sides = entries.slice(0, -1).map(([, bytes]) => ncbiHash(toFrozen(bytes, search.frozen_prefix))).join(",");
          if (search.frozen.custom) row.frozen.custom = master === search.frozen.custom;
          if (search.frozen.sides !== undefined) row.frozen.sides = sides === search.frozen.sides;
        } else if (customStream(search)) {
          const sent = run.streams.get(CUSTOM_STREAM) || Buffer.alloc(0);
          const actual = search.sam ? withoutCl(sent) : sent;
          const cli = expected.outputs.get(CUSTOM_STREAM);
          assert.equal(sha256(actual), sha256(search.sam ? withoutCl(cli) : cli), `${label}: -outfmt ${search.outfmt}`);
          row.formats.custom = sha256(actual);
          if (search.frozen.custom) row.frozen.custom = ncbiHash(toFrozen(actual, search.frozen_prefix)) === search.frozen.custom;
        } else {
          assert.ok(!run.streams.has(CUSTOM_STREAM), `${label}: no stream ${CUSTOM_STREAM} without an -outfmt of its own`);
        }
        const diagnostics = run.streams.get(3) || Buffer.alloc(0);
        assert.equal(sha256(diagnostics), sha256(expected.stderr), `${label}: diagnostics\n${diagnostics}`);
        row.diagnostics = sha256(expected.stderr);
        // The run's signal of stream 10 says where its messages go (abi_v2.md §5; line WQA's sweep, D7, 2026-11-20): one
        // signal of 4 bytes; stream 3 without -logfile; the file of -logfile only where the run gave it.
        const messages = run.streams.get(10) || Buffer.alloc(0);
        assert.equal(messages.length, 4, `${label}: one signal of 4 bytes (stream 10)`);
        const messagesPlace = messages.readUInt32LE(0);
        assert.ok([3, 8, 0xffffffff].includes(messagesPlace), `${label}: the place of the run's messages (${messagesPlace})`);
        if (!search.argv.some((word) => word === "-logfile" || word.startsWith("-logfile="))) {
          assert.equal(messagesPlace, 3, `${label}: the run's messages go to stream 3 without -logfile`);
        }
        // The run's signals of stream 9 name the streams that it could not write: those whose CLI run of the argv
        // fails (abi_v2.md §5; line WQA's sweep, D7, 2026-11-19), each with where its error is (2026-11-20): stream 3,
        // which then holds it, or nowhere where the run's messages go nowhere; never the file of -logfile, which holds
        // what the CLI's file holds for the argv (line WQA's sweep, D3 review r2, 2026-11-19).
        const signals = run.streams.get(9) || Buffer.alloc(0);
        assert.equal(signals.length % 8, 0, `${label}: signals of 8 bytes (stream 9)`);
        const failed = Array.from({ length: signals.length / 8 }, (_, i) => [signals.readUInt32LE(8 * i), signals.readUInt32LE(8 * i + 4)]);
        const failedStreams = failed.map(([stream]) => stream).sort((a, b) => a - b);
        assert.deepEqual(failedStreams, expected.failedStreams, `${label}: the streams that the run says failed (stream 9)`);
        for (const [stream, place] of failed) {
          if (place === 3) assert.ok(diagnostics.length > 0, `${label}: the error of stream ${stream} is on stream 3, which is empty`);
          else {
            assert.equal(place, 0xffffffff, `${label}: the place of the error of stream ${stream}`);
            assert.equal(messagesPlace, 0xffffffff, `${label}: the error of stream ${stream} is nowhere, where the run's messages go`);
          }
        }
        if (messagesPlace === 8) assert.ok(run.streams.has(8), `${label}: the run's messages are in the file of -logfile, which the run did not give`);
        if (search.fails_alone) row.fails_alone = search.fails_alone;
        row.hits = WITH_HITS.has(program) ? checkHits(run.streams, program, label, queryIds, shownIds, search.database, search.argv.includes("-html")) : null;
        if (together) row.inputs = "register_inputs";
        results.push(row);
        console.log(`ok ${label}`);
      }
      if (q) assert.equal(reactor.call("losat_web2_release", q.status).status, 0);
      if (registry !== null) assert.equal(reactor.call("losat_web2_file_release", registry).status, 0);
    }
    for (const [key, handle] of subjects) assert.equal(subjectReactors.get(key).call("losat_web2_release", handle).status, 0);
    for (const [holder, held] of counts) for (const { handle } of held.values()) assert.equal(holder.call("losat_web2_file_release", handle).status, 0);
    for (const instance of databases.values()) await instance.close();
    for (const [, instance] of tools.values()) await instance.close();
    if (countsRegistered + countsKept) console.log(`ok ${engine.name} -window_masker_db and -window_masker_taxid: ${countsRegistered} searches that registered counts, ${countsKept} searches with counts kept or unregistered`);
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
