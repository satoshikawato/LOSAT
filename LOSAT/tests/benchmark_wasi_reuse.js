"use strict";
const fs = require("node:fs");
const assert = require("node:assert/strict");
const crypto = require("node:crypto");
const { performance } = require("node:perf_hooks");
const { WASI } = require("node:wasi");
const { inspectArtifact } = require("./wasi_artifact");
const { prepareThreadHost } = require("./wasi_thread_host");
const { runPair } = require("./check_wasi_reactor");
const hash = bytes => crypto.createHash("sha256").update(bytes).digest("hex");

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1410-1475
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) { ... }
// Each sample completes one search with the same inputs. Compile reuse and
// reactor instance reuse are separate measurements with explicit boundaries.
(async () => {
  const [artifact, kind, mode, jobsFile, resultFile] = process.argv.slice(2);
  const jobs = JSON.parse(fs.readFileSync(jobsFile, "utf8"));
  let prepared, module, timings, identity;
  if (kind.startsWith("threaded")) {
    prepared = await prepareThreadHost(artifact, kind); ({ timings, identity } = prepared);
  } else {
    const start = performance.now(), bytes = fs.readFileSync(artifact);
    identity = inspectArtifact(bytes, kind);
    const inspected = performance.now(); module = await WebAssembly.compile(bytes);
    timings = { validation_including_raw_compile_seconds:(inspected-start)/1000, guard_seconds:0, second_compile_seconds:(performance.now()-inspected)/1000 };
  }
  // NCBI reference: c++/src/algo/blast/api/seqsrc_multiseq.cpp:175-180
  // m_iTotalLength += (Int8) (*iter)->length;
  // Explicit API inputs are prepared once per case, outside invocation timing.
  // Command jobs retain file paths; repeated job metadata never duplicates FASTA text.
  const inputStart = performance.now(), inputs = new Map();
  if (kind.endsWith("reactor")) for (const job of jobs) {
    if (!inputs.has(job.case_id)) inputs.set(job.case_id, {
      query:fs.readFileSync(job.query_file,"utf8"), subject:fs.readFileSync(job.subject_file,"utf8"),
    });
  }
  timings.input_load_seconds=(performance.now()-inputStart)/1000;
  timings.resident_input_bytes=[...inputs.values()].reduce((n,v)=>n+Buffer.byteLength(v.query)+Buffer.byteLength(v.subject),0);
  const instances = new Map(), samples = [];
  const save = () => fs.writeFileSync(resultFile, JSON.stringify({ identity, timings, mode, samples }, null, 2));
  async function create(job) {
    if (prepared) return prepared.create(job.argv || []);
    const wasi = new WASI({version:"preview1", args:[artifact,...(job.argv||[])], env:process.env, preopens:{"/":"/"}, returnOnExit:true});
    const instance = await WebAssembly.instantiate(module, {wasi_snapshot_preview1:wasi.wasiImport});
    if (kind.endsWith("reactor")) wasi.initialize(instance);
    return {instance, memory:instance.exports.memory, events:[], start:()=>wasi.start(instance), waitForWorkers:async()=>{}, close:async()=>{}};
  }
  try {
    for (const job of jobs) {
      const key = job.case_id;
      let host = mode === "same-instance" ? instances.get(key) : null;
      const start = performance.now(), cpu = process.cpuUsage();
      if (!host) {
        host = await create(job);
        if (mode === "same-instance") instances.set(key, host);
      }
      const instantiated = performance.now();
      let data, status;
      if (kind.endsWith("reactor")) {
        const response = runPair(host, job.program, inputs.get(key).query, inputs.get(key).subject, job.format || "6", job.extra);
        status = response.status; data = response.result;
        assert.equal(status, 0, response.error);
        fs.writeFileSync(job.output, data);
      } else {
        status = host.start(); data = fs.readFileSync(job.output);
      }
      const returned = performance.now();
      const events = await host.waitForWorkers() || [];
      const exited = performance.now();
      const n = job.threads;
      const counts = Object.fromEntries(["spawn_attempt","spawned","ready","exited"].map(e=>[e,events.filter(x=>x.event===e).length]));
      if (prepared) for (const count of Object.values(counts)) assert.equal(count,n===1?0:n);
      assert.ok(events.filter(e=>e.event==="exited").every(e=>e.code===0));
      const usage = process.cpuUsage(cpu);
      const sample = { case_id:key, repeat:job.repeat, timed:job.timed, threads:n, status,
        boundary:mode, wall_seconds:(exited-start)/1000, instantiate_seconds:(instantiated-start)/1000,
        invocation_including_io_and_api_copy_seconds:(returned-instantiated)/1000, host_exit_wait_seconds:(exited-returned)/1000,
        cpu_user_seconds:usage.user/1e6, cpu_system_seconds:usage.system/1e6,
        process_lifetime_peak_rss_bytes:process.resourceUsage().maxRSS*1024,
        rss_after_bytes:process.memoryUsage().rss, memory_bytes:host.memory.buffer.byteLength,
        raw_output_sha256:hash(data), output:job.output, events, host_worker_counts:counts };
      samples.push(sample); save();
      assert.equal(status,0); assert.equal(sample.raw_output_sha256,job.expected_sha256,`raw mismatch: ${key}`);
      if (mode !== "same-instance") await host.close();
    }
  } finally { for (const host of instances.values()) await host.close(); save(); }
})().catch(error=>{console.error(error);process.exitCode=1;});
