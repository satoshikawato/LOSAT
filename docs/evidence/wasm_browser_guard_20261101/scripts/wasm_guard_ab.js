"use strict";
// wg-wasm-guard (2026-11-01): in-process A/B of the shared-memory guard on the threaded Wasm command.
// ONE Node process, two arms alternating (ABBA: round r even = guarded then raw, odd = raw then guarded), per case
// per round. Each search = createThreadHost + start + waitForWorkers + close of the command (host start-up, guard,
// compile and the search itself; that is what a user of run_losat_wasi_threads.js pays).
//   guarded = the build-time guarded artifact (custom section losat.shared_memory_guard; the host uses it as it is)
//   raw     = the unguarded cargo output; the host's guardSharedMemory is replaced by a pass-through (no guard at all)
// usage: node wasm_guard_ab.js --rounds N --threads T --out FILE.jsonl [--cases bp,bpff,...] [--outdir DIR]
//          [--host DIR] [--guarded WASM] [--raw WASM] [--debug 1] [--warmup 1|0] [--first-round R]
// One JSON line per search: case, arm, round, wall_s, exit, error, out_sha256, out_bytes, load1/5/15, guard_s, compile_s,
// compiled_sha256, compiled_ok (the bytes handed to WebAssembly.compile are those of the intended arm's file).
// A first line {"header":...} records the file hashes. A trap/failure is a failed sample (exit != 0, error set), the run goes on.
const fs = require("node:fs");
const path = require("node:path");
const crypto = require("node:crypto");
const { performance } = require("node:perf_hooks");

const opt = {};
for (let i = 2; i < process.argv.length; i += 2) opt[process.argv[i].replace(/^--/, "")] = process.argv[i + 1];
const HOST = opt.host || "/home/kawato/losat-work/.worktrees/wg-guard/LOSAT/tests";
const BASE = "/home/kawato/losat-baselines/wg-wasm-guard-20261101";
const ARMS = {
  guarded: opt.guarded || `${BASE}/art/losat-threaded-command.wasm`,
  raw: opt.raw || "/home/kawato/.cache/losat-work/wt-wg-guard/threaded-command/wasm32-wasip1-threads/release/losat.wasm",
};
const ROUNDS = Number(opt.rounds || 1), THREADS = Number(opt.threads || 4);
const OUT = opt.out; if (!OUT) throw new Error("--out FILE.jsonl required");
const OUTDIR = opt.outdir || process.env.TMPDIR || `${BASE}/tmp`;
if (opt.debug === "1") process.env.LOSAT_WASI_THREADS_DEBUG = "1";

const inherited = Object.keys(process.env).filter(n => n.startsWith("LOSAT_") && n !== "LOSAT_WASI_THREADS_DEBUG");
if (inherited.length) throw new Error(`LOSAT_* must not come from outside: ${inherited}`);
process.env.LOSAT_WASI_GUARD_CACHE = "0";   // no on-disk guard cache: the raw arm guards nothing, the guarded arm needs none

// ---- cases (line W's 9 daily cases; argv of v4-progress-20261025/ws/vp_cases.sh largs, inputs mapped to local copies) ----
const L = "/home/kawato/losat-baselines";
const T = `${BASE}/timing/data`;   // re-wrapped copies, sha256 identical to the workstation files
const I = {
  sp4_first3000: `${T}/sp4_first3000.fa`,
  sp4: `${L}/perf-of-20261013/ws/full7_a3/Species4.fa`, sp6: `${L}/perf-of-20261013/ws/full7_a3/Species6.fa`,
  q100: `${L}/v4-W/tmp/s18in/q100.fna`, s16: `${L}/v4-W/tmp/s18in/16s_all.fna`,
  cdna: `${L}/v4-A10/wasm_in/dmel_cdna_first31168.fa`, dna3m: `${L}/v4-A10/wasm_in/dmel_dna_3M.fa`,
  short: `${L}/v4-W/tmp/s18in/short_q30000.fa`,
  g2l: `${L}/v4-W/tmp/bxdata/g2l50k.fa`, pepdb: `${L}/v4-W/tmp/bxdata/dmel_pep`,
  txq: `${T}/txB_q.fa`, txs: `${T}/txB_s.fa`,
  // tblastn -db: hsap_p1.faa + the hsap_dna BLAST DB exist only on the workstation (work/tbn/in, work/tbn/db)
  tbq: process.env.WG_TBN_QUERY || null, tbdb: process.env.WG_TBN_DB || null,
};
const CASES = {
  bp:   ["blastp", "-query", I.sp4_first3000, "-subject", I.sp6, "-outfmt", "6", "-evalue", "0.001"],
  bpff: ["blastp", "-task", "blastp-fast", "-query", I.sp4, "-subject", I.sp6, "-outfmt", "6", "-evalue", "0.001"],
  mb:   ["blastn", "-task", "megablast", "-query", I.q100, "-subject", I.s16, "-outfmt", "6"],
  bn:   ["blastn", "-task", "blastn", "-query", I.cdna, "-subject", I.dna3m, "-outfmt", "6"],
  dc:   ["blastn", "-task", "dc-megablast", "-query", I.cdna, "-subject", I.dna3m, "-outfmt", "6"],
  bns:  ["blastn", "-task", "blastn-short", "-query", I.short, "-subject", I.dna3m, "-outfmt", "6", "-evalue", "0.001"],
  tbn:  ["tblastn", "-query", I.tbq, "-db", I.tbdb, "-evalue", "1e-5", "-max_target_seqs", "5", "-outfmt", "6"],
  bx:   ["blastx", "-query", I.g2l, "-db", I.pepdb, "-evalue", "1e-5", "-outfmt", "6"],
  tbx:  ["tblastx", "-query", I.txq, "-subject", I.txs, "-outfmt", "6", "-query_gencode", "4", "-db_gencode", "4"],
};
const ORDER = ["bp", "bpff", "mb", "bn", "dc", "bns", "tbn", "bx", "tbx"];
const wanted = (opt.cases || ORDER.join(",")).split(",");
const cases = [], skipped = [];
for (const c of wanted) {
  if (!CASES[c]) throw new Error(`unknown case ${c}`);
  const missing = CASES[c].filter(a => a === null || (a.startsWith("/") && !fs.existsSync(a) && !fs.existsSync(a + ".pin") && !fs.existsSync(a + ".phr") && !fs.existsSync(a + ".nin")));
  if (CASES[c].includes(null) || missing.length) skipped.push(c); else cases.push(c);
}

// ---- host: replace guardSharedMemory BEFORE wasi_thread_host.js is required (it destructures at require time) ----
const sm = require(path.join(HOST, "wasi_shared_memory.js"));
const origGuard = sm.guardSharedMemory;
let armNow = null, stat = null;
sm.guardSharedMemory = (bytes) => {
  const t = performance.now();
  const r = armNow === "raw" ? { bytes, fill: 0, copy: 0, record: null } : origGuard(bytes);
  if (stat) stat.guard_s += (performance.now() - t) / 1000;
  return r;
};
const origCompile = WebAssembly.compile;
WebAssembly.compile = function (bytes, ...rest) {
  const t = performance.now(); if (stat) stat.compiled.push(bytes);
  const p = origCompile.call(this, bytes, ...rest);
  return p.finally(() => { if (stat) stat.compile_s += (performance.now() - t) / 1000; });
};
const { createThreadHost } = require(path.join(HOST, "wasi_thread_host.js"));

// The host ends the process by SIGTERM when a thread does not start within 30 s (awaitThreadStart), and a worker ends
// it the same way when its thread traps (wasi_thread_worker.js failWorker). Do not catch it: the main thread may be
// blocked in Wasm on a latch of the dead thread, so a caught SIGTERM left the process spinning (2026-11-01, 7 h at
// ~700 % CPU after a trap of the raw arm). run_ab.sh records the search of OUT.current as failed and goes on with the
// next round in a new process (--first-round).
const sigterms = 0;
const CURRENT = OUT + ".current";
const sha = (b) => crypto.createHash("sha256").update(b).digest("hex");
const load = () => fs.readFileSync("/proc/loadavg", "utf8").split(" ").slice(0, 3).map(Number);
const emit = (row) => fs.appendFileSync(OUT, JSON.stringify(row) + "\n");

async function search(arm, round, c) {
  const outFile = path.join(OUTDIR, `wgab.${process.pid}.${c}.${arm}.out`);
  fs.rmSync(outFile, { force: true });
  const argv = [...CASES[c], "-num_threads", String(THREADS), "-out", outFile];
  const [l1, l5, l15] = load();
  armNow = arm; stat = { guard_s: 0, compile_s: 0, compiled: [] };
  fs.writeFileSync(CURRENT, JSON.stringify({ case: c, arm, round, threads: THREADS, load1: l1, load5: l5, load15: l15 }) + "\n");
  let exit = null, error = "", host = null;
  const sig0 = sigterms;
  const t0 = performance.now();
  try {
    host = await createThreadHost(ARMS[arm], "threaded-command", argv);
    exit = host.start();
    await host.waitForWorkers();
  } catch (err) { error = String(err && err.message || err).slice(0, 300); if (exit === null) exit = -1; }
  const wall = (performance.now() - t0) / 1000;
  if (sigterms > sig0) { error = `host SIGTERM (thread start timeout) ${error}`.trim(); if (exit === 0 || exit === null) exit = -2; }
  try { if (host) await Promise.race([host.close(), new Promise(r => setTimeout(r, 20000))]); } catch (err) { error += ` close:${err.message}`; }
  const s = stat; stat = null; armNow = null;
  // untimed: what the host compiled, and the output
  const shas = s.compiled.map(b => sha(b));
  const intended = headerShas[arm];
  let outSha = null, outBytes = null;
  if (fs.existsSync(outFile)) { const b = fs.readFileSync(outFile); outSha = sha(b); outBytes = b.length; fs.rmSync(outFile, { force: true }); }
  emit({ case: c, arm, round, wall_s: Number(wall.toFixed(4)), exit, error, out_sha256: outSha, out_bytes: outBytes,
         load1: l1, load5: l5, load15: l15, threads: THREADS, guard_s: Number(s.guard_s.toFixed(4)), compile_s: Number(s.compile_s.toFixed(4)),
         compiled_sha256: shas, compiled_ok: shas.length === 1 && shas[0] === intended });
  console.error(`r${round} ${c} ${arm}: ${wall.toFixed(3)} s exit ${exit}${error ? " ERR " + error : ""} out ${outSha ? outSha.slice(0, 12) : "none"} guard ${s.guard_s.toFixed(3)} compile ${s.compile_s.toFixed(3)} compiled_ok ${shas.length === 1 && shas[0] === intended}`);
}

const headerShas = {};
(async () => {
  const gBytes = fs.readFileSync(ARMS.guarded), rBytes = fs.readFileSync(ARMS.raw);
  headerShas.guarded = sha(gBytes); headerShas.raw = sha(rBytes);
  const regen = origGuard(rBytes);   // untimed: the guard applied to the raw bytes (same cargo output?)
  emit({ header: { node: process.version, host: HOST, rounds: ROUNDS, threads: THREADS, cases, skipped, arms: ARMS,
    file_sha256: headerShas, guarded_has_record: sm.guardRecord(gBytes) !== null, raw_has_record: sm.guardRecord(rBytes) !== null,
    guard_of_raw_equals_guarded_file: sha(regen.bytes) === headerShas.guarded, guard_of_raw_counts: { fill: regen.fill, copy: regen.copy },
    started: new Date().toISOString(), load: load() } });
  if (skipped.length) console.error(`skipped (inputs missing): ${skipped}`);
  // --warmup 1 (default): one discarded pair per case first (round -1; the first search of a process reads the inputs cold, x4 seen on dc)
  if (Number(opt.warmup ?? 1) > 0) for (const c of cases) for (const arm of ["guarded", "raw"]) await search(arm, -1, c);
  for (let r = Number(opt["first-round"] || 0); r < ROUNDS; r++) {
    const order = r % 2 === 0 ? ["guarded", "raw"] : ["raw", "guarded"];
    for (const c of cases) for (const arm of order) await search(arm, r, c);
  }
  fs.rmSync(CURRENT, { force: true });
})().catch(err => { console.error(err); process.exitCode = 1; });
