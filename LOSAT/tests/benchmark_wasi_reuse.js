"use strict";
const fs = require("node:fs");
const assert = require("node:assert/strict");
const crypto = require("node:crypto");
const { performance } = require("node:perf_hooks");
const { WASI } = require("node:wasi");
const { prepareArtifact } = require("./wasi_artifact");
const { prepareThreadHost } = require("./wasi_thread_host");
const { runPair } = require("./check_wasi_reactor");
const hash = bytes => crypto.createHash("sha256").update(bytes).digest("hex");

// NCBI reference: c++/src/algo/blast/core/blast_engine.c:1410-1475
// while ((seq_arg.oid = BlastSeqSrcIteratorNext(seq_src, itr)) != BLAST_SEQSRC_EOF) { ... }
// Each sample completes one search with the same inputs. Compile reuse and
// reactor instance reuse are separate measurements with explicit boundaries.
(async () => {
  const [artifact, kind, mode, jobsFile, outputFlag, outputPath] = process.argv.slice(2);
  const resultFile = outputFlag === "-out" ? outputPath : outputFlag;
  assert.ok(["compiled-module", "same-instance"].includes(mode));
  assert.ok(mode !== "same-instance" || kind.endsWith("reactor"), "commands cannot reuse an exited instance");
  const jobs = JSON.parse(fs.readFileSync(jobsFile, "utf8"));
  // NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
  // x_PrintField(*iter); ... m_Ostream << "\\n";
  // Never let a stale result supply bytes for a successful command invocation.
  fs.closeSync(fs.openSync(resultFile, "wx"));
  const paths = new Set();
  for (const job of jobs) {
    assert.ok(!paths.has(job.output) && !fs.existsSync(job.output), `existing output: ${job.output}`);
    paths.add(job.output);
  }
  let prepared, module, timings, identity;
  if (kind.startsWith("threaded")) {
    prepared = await prepareThreadHost(artifact, kind); ({ timings, identity } = prepared);
  } else {
    const start = performance.now(), bytes = fs.readFileSync(artifact);
    // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
    // CATCH_ALL(status); ... return status;
    // Serial execution owns the inspected module; no second compile is requested.
    ({ module, identity } = prepareArtifact(bytes, kind));
    timings = { validation_including_raw_compile_seconds:(performance.now()-start)/1000, guard_seconds:0, second_compile_seconds:0 };
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
  const instances = new Map();
  const samples = jobs.map(job => ({ case_id:job.case_id, repeat:job.repeat, timed:job.timed,
    threads:job.threads, status:"NOT_RUN", eligible_for_timing:false, output:job.output }));
  const save = () => fs.writeFileSync(resultFile, JSON.stringify({ identity, timings, mode,
    node_argv:[process.execPath, ...process.execArgv], instance_scope:mode === "same-instance" ? "per-case" : "per-job",
    samples }, null, 2));
  save();
  async function create(job) {
    if (prepared) return prepared.create(job.argv || []);
    const wasi = new WASI({version:"preview1", args:[artifact,...(job.argv||[])], env:process.env, preopens:{"/":"/"}, returnOnExit:true});
    const instance = await WebAssembly.instantiate(module, {wasi_snapshot_preview1:wasi.wasiImport});
    if (kind.endsWith("reactor")) wasi.initialize(instance);
    return {instance, memory:instance.exports.memory, events:[], start:()=>wasi.start(instance), waitForWorkers:async()=>{}, close:async()=>{}};
  }
  // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:172-188
  // (*thread)->Run(); ... (*thread)->Join(&result);
  // Complete lifecycle, exact bytes, and success status must all hold for timing.
  try {
    for (const [index, job] of jobs.entries()) {
      const sample = samples[index], key = job.case_id;
      let host = mode === "same-instance" ? instances.get(key) : null;
      sample.status = "RUNNING"; save();
      // Evidence serialization is outside the measured runtime boundary.
      // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:178-188
      // (*thread)->Join(&result);
      // Capture wall-clock evidence at the same endpoints as the monotonic timer.
      const realtimeStart = Date.now(), start = performance.now(), cpu = process.cpuUsage();
      let runtimeEnd, realtimeEnd, runtimeCpu, failure = null;
      try {
        if (!host) {
          host = await create(job);
          if (mode === "same-instance") instances.set(key, host);
        }
        const instantiated = performance.now();
        let data, exitStatus;
        if (kind.endsWith("reactor")) {
          const response = runPair(host, job.program, inputs.get(key).query, inputs.get(key).subject, job.format || "6", job.extra);
          exitStatus = response.status; data = response.result;
          sample.exit_status = exitStatus;
          assert.equal(exitStatus, 0, response.error);
          fs.writeFileSync(job.output, data, {flag:"wx"});
        } else {
          exitStatus = host.start(); sample.exit_status = exitStatus;
          assert.equal(exitStatus, 0);
          data = fs.readFileSync(job.output);
        }
        const returned = performance.now();
        const events = await host.waitForWorkers() || [];
        const exited = performance.now();
        runtimeEnd = exited; realtimeEnd = Date.now(); runtimeCpu = process.cpuUsage(cpu);
        const n = job.threads;
        const counts = Object.fromEntries(["spawn_attempt","spawned","ready","exited"].map(e=>[e,events.filter(x=>x.event===e).length]));
        Object.assign(sample, { instantiate_seconds:(instantiated-start)/1000,
          invocation_including_io_and_api_copy_seconds:(returned-instantiated)/1000,
          host_exit_wait_seconds:(exited-returned)/1000, raw_output_sha256:hash(data),
          expected_sha256:job.expected_sha256, raw_equal:hash(data) === job.expected_sha256,
          events, host_worker_counts:counts, memory_bytes:host.memory.buffer.byteLength });
        if (prepared) {
          // NCBI reference: c++/src/algo/blast/core/blast_kappa.c:3429-3459
          // #pragma omp parallel ... num_threads(actual_num_threads)
          // LOSAT's existing total-thread contract includes the caller.
          const expected = n - 1;
          const tids = events.filter(e=>e.event === "spawn_attempt").map(e=>e.tid).sort((a,b)=>a-b);
          for (const event of Object.keys(counts)) {
            assert.equal(counts[event], expected);
            assert.deepEqual(events.filter(e=>e.event===event).map(e=>e.tid).sort((a,b)=>a-b), tids);
          }
          assert.equal(new Set(tids).size, expected);
        }
        assert.ok(events.filter(e=>e.event==="exited").every(e=>e.code===0));
        sample.thread_contract = "PASS";
        assert.ok(sample.raw_equal, `raw mismatch: ${key}`);
        sample.status = "PASS";
      } catch (error) {
        failure = error; sample.status = sample.raw_equal === false ? "PARITY_FAIL" : "FAIL";
        sample.reason = String(error.stack || error);
      } finally {
        runtimeEnd ??= performance.now(); realtimeEnd ??= Date.now(); runtimeCpu ??= process.cpuUsage(cpu);
        const closingRealtime = Date.now(), closing = performance.now(), closingCpu = process.cpuUsage();
        if (host && mode !== "same-instance") {
          try { await host.close(); } catch (error) {
            failure ||= error; sample.status = "FAIL"; sample.close_error = String(error);
          }
        }
        const closeSeconds = (performance.now()-closing)/1000, closeCpu = process.cpuUsage(closingCpu);
        const wallSeconds = (runtimeEnd-start)/1000 + closeSeconds;
        const realtimeSeconds = (realtimeEnd-realtimeStart + Date.now()-closingRealtime)/1000;
        const clockAgreement = Math.abs(realtimeSeconds-wallSeconds) <= Math.max(0.1, 0.05*wallSeconds);
        Object.assign(sample, { boundary:mode, wall_seconds:wallSeconds,
          realtime_seconds:realtimeSeconds, realtime_clock_agreement:clockAgreement,
          elapsed_clock:"performance.now() (CLOCK_MONOTONIC); Date.now() is an adjustable-clock diagnostic",
          close_seconds:closeSeconds, cpu_user_seconds:(runtimeCpu.user+closeCpu.user)/1e6,
          cpu_system_seconds:(runtimeCpu.system+closeCpu.system)/1e6,
          process_lifetime_peak_rss_bytes:process.resourceUsage().maxRSS*1024, rss_after_bytes:process.memoryUsage().rss,
          eligible_for_timing:sample.status === "PASS" && job.timed });
        save();
      }
      if (failure) throw failure;
    }
  } finally {
    const closing = performance.now();
    for (const host of instances.values()) await host.close();
    timings.final_instance_close_seconds = (performance.now()-closing)/1000;
    save();
  }

})().catch(error=>{console.error(error);process.exitCode=1;});
