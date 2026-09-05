"use strict";
const fs = require("fs");
const { WASI } = require("wasi");
const { performance } = require("perf_hooks");
const crypto = require("crypto");

// NCBI c++/src/algo/blast/core/blast_engine.c:1410-1413:
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF)
// NCBI c++/src/algo/blast/format/blast_format.cpp:828-832:
// ITERATE(CSeq_align_set::Tdata, itr, copy_aln_set.Get()) { ... tabinfo.Print(); }
// Each command performs a complete unchanged search. Only the compiled module
// is shared: every job gets a fresh instance, WASI object, and linear memory.
// Never call _start twice on an instance or label this a resident search engine.
(async () => {
  const [, , wasmPath, jobsPath] = process.argv;
  const jobs = JSON.parse(fs.readFileSync(jobsPath, "utf8"));
  const compileStart = performance.now();
  const module = await WebAssembly.compile(fs.readFileSync(wasmPath));
  const compileSeconds = (performance.now() - compileStart) / 1000;
  if (WebAssembly.Module.imports(module).some((x) => x.module === "wasi")) {
    throw new Error("warm serial runner rejects a threaded command");
  }
  const samples = [];
  for (const job of jobs) {
    const start = performance.now();
    const cpu = process.cpuUsage();
    const wasi = new WASI({
      version: "preview1", args: job.argv, env: process.env,
      preopens: { "/": "/" }, returnOnExit: true,
    });
    const instance = await WebAssembly.instantiate(module, {
      wasi_snapshot_preview1: wasi.wasiImport,
    });
    const instanceSeconds = (performance.now() - start) / 1000;
    if (typeof instance.exports._start !== "function") {
      throw new Error("command _start missing");
    }
    const commandStart = performance.now();
    const exit = wasi.start(instance);
    const elapsed = (performance.now() - start) / 1000;
    const cpuElapsed = process.cpuUsage(cpu);
    const sample = {
      repeat: job.repeat, ordered_argv: job.argv, exit_status: exit,
      boundary: "cached-module-new-instance-through-command-return",
      wall_seconds: elapsed, instantiate_seconds: instanceSeconds,
      command_seconds: (performance.now() - commandStart) / 1000,
      cpu_user_seconds: cpuElapsed.user / 1e6,
      cpu_system_seconds: cpuElapsed.system / 1e6,
      rss_after_bytes: process.memoryUsage().rss,
      process_lifetime_peak_rss_bytes: process.resourceUsage().maxRSS * 1024,
      linear_memory_after_bytes: instance.exports.memory.buffer.byteLength,
      output: job.output,
    };
    samples.push(sample);
    fs.writeFileSync(process.argv.at(-1), JSON.stringify({ compile_seconds: compileSeconds, samples }));
    if (exit !== 0) throw new Error(`command exit ${exit}`);
    const hash = crypto.createHash("sha256").update(fs.readFileSync(job.output)).digest("hex");
    if (hash !== job.expected_sha256) throw new Error("raw gate failed; refusing subsequent repeats");
  }
})().catch((error) => {
  console.error(error.stack || String(error));
  process.exitCode = 1;
});
