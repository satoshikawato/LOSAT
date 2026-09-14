"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { createThreadHost } = require("./wasi_thread_host");

// NCBI reference: c++/src/algo/blast/format/blast_format.cpp:68-96,770-832
// m_Outfile(ostr); CBlastTabularInfo tabinfo(m_Outfile, ...);
// Preserve raw API bytes, including each format's existing labels.
function runPair(host, program, query, subject, format, args) {
  const api = host.instance.exports;
  const values = [program, query, subject, format, args.join("\0")];
  const allocations = values.map(value => {
    const bytes = Buffer.from(value), ptr = api.losat_web_alloc(bytes.length);
    new Uint8Array(host.memory.buffer, ptr, bytes.length).set(bytes);
    return [ptr, bytes.length];
  });
  try {
    const status = api.losat_web_run_pair(...allocations.flat());
    const result = Buffer.from(new Uint8Array(host.memory.buffer, api.losat_web_result_ptr(), api.losat_web_result_len()));
    const error = Buffer.from(new Uint8Array(host.memory.buffer, api.losat_web_error_ptr(), api.losat_web_error_len())).toString();
    return { status, result, error };
  } finally { for (const [ptr, length] of allocations) api.losat_web_dealloc(ptr, length); }
}

// NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
// (*thread)->Run(); (*thread)->Join(&result);
// Guest joins precede the API return; host exit callbacks are awaited explicitly.
async function checkReactor(artifact, fixtures, output) {
  fs.mkdirSync(output, { recursive: true });
  const nuc = fs.readFileSync(path.join(fixtures, "nuc3.fasta"), "utf8");
  const single = fs.readFileSync(path.join(fixtures, "nuc1.fasta"), "utf8");
  const aa = fs.readFileSync(path.join(fixtures, "aa3.fasta"), "utf8");
  const cases = [
    { id: "blastn-dp", program: "blastn", seq: nuc, args: ["-task", "blastn", "-word_size", "11"] },
    { id: "megablast", program: "blastn", seq: nuc, args: ["-task", "megablast"] },
    { id: "tblastx-linking", program: "tblastx", seq: single, args: [] },
    { id: "blastn-single-dp", program: "blastn", seq: single, args: ["-task", "blastn"] },
    { id: "tblastx", program: "tblastx", seq: nuc, args: [] },
    { id: "blastp", program: "blastp", seq: aa, args: [] },
  ];
  const host = await createThreadHost(artifact, "threaded-reactor");
  const records = [];
  const save = (id, response, events) => {
    fs.writeFileSync(path.join(output, `${id}.out`), response.result);
    fs.writeFileSync(path.join(output, `${id}.error`), response.error);
    records.push({ id, status: response.status, events, memory_bytes: host.memory.buffer.byteLength });
    fs.writeFileSync(path.join(output, "runs.json"), JSON.stringify(records, null, 2) + "\n");
  };
  try {
    for (const c of cases) {
      let reference;
      for (const [i, n] of [1, 2, 4, 8, 2, 1, 2, 2].entries()) {
        const response = runPair(host, c.program, c.seq, c.seq, "6", [...c.args, "-num_threads", String(n)]);
        const events = await host.waitForWorkers();
        save(`${c.id}-${i}-n${n}`, response, events);
        assert.equal(response.status, 0, response.error);
        if (!reference) reference = response.result;
        assert.deepEqual(response.result, reference, c.id);
        for (const event of ["spawn_attempt", "spawned", "ready", "exited"]) {
          assert.equal(events.filter(e => e.event === event).length, n - 1, `${c.id} n${n} ${event}`);
        }
        assert.ok(events.filter(e => e.event === "exited").every(e => e.code === 0));
      }
      // Failure must clear previous result bytes and retain the engine cause.
      for (const bad of ["0", "invalid", "256"]) {
        const response = runPair(host, c.program, c.seq, c.seq, "6", [...c.args, "-num_threads", bad]);
        const events = await host.waitForWorkers(); save(`${c.id}-reject-${bad}`, response, events);
        assert.equal(response.status, -1); assert.equal(response.result.length, 0);
        assert.ok(response.error.length > 0); assert.equal(events.length, 0);
      }
      for (const failAt of [1, 2]) {
        let attempt = 0; host.setRejectSpawn(() => ++attempt === failAt);
        const response = runPair(host, c.program, c.seq, c.seq, "6", [...c.args, "-num_threads", "4"]);
        const events = await host.waitForWorkers();
        save(`${c.id}-spawn-fail-${failAt}`, response, events);
        assert.equal(response.status, -1); assert.equal(response.result.length, 0);
        assert.match(response.error, /failed to build.*pool.*threads.*:/);
        assert.equal(events.filter(e => e.event === "ready").length, failAt - 1);
        assert.equal(events.filter(e => e.event === "exited").length, failAt - 1);
        host.setRejectSpawn(null);
        const recovery = runPair(host, c.program, c.seq, c.seq, "6", [...c.args, "-num_threads", "2"]);
        const recoveryEvents = await host.waitForWorkers(); save(`${c.id}-recovery-${failAt}`, recovery, recoveryEvents);
        assert.equal(recovery.status, 0, recovery.error); assert.deepEqual(recovery.result, reference);
      }
      // NCBI reference: c++/src/algo/blast/format/blast_format.cpp:770-832
      // CBlastTabularInfo tabinfo(m_Outfile, ...);
      // API labels are a separate existing contract; compare unchanged bytes.
      for (const format of c.program === "blastp" ? ["0", "7", "6 qseqid sseqid score bitscore qstart qend sstart send"] : c.program === "blastn" ? ["7"] : []) {
        let expected;
        for (const n of [1, 2, 4, 8]) {
          const response = runPair(host, c.program, c.seq, c.seq, format, [...c.args, "-num_threads", String(n)]);
          const events = await host.waitForWorkers(); save(`${c.id}-format-${format.split(" ")[0]}-n${n}`, response, events);
          assert.equal(response.status, 0, response.error);
          if (!expected) expected = response.result;
          assert.deepEqual(response.result, expected);
          // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
          // (*thread)->Run(); (*thread)->Join(&result);
          for (const event of ["spawn_attempt", "spawned", "ready", "exited"]) {
            assert.equal(events.filter(e => e.event === event).length, n - 1);
          }
          assert.ok(events.filter(e => e.event === "exited").every(e => e.code === 0));
        }
      }
      const invalid = runPair(host, c.program, "invalid FASTA", c.seq, "6", [...c.args, "-num_threads", "1"]);
      assert.equal(invalid.status, -1); assert.equal(invalid.result.length, 0);
      assert.match(invalid.error, /FASTA|fasta/);
      save(`${c.id}-invalid-input`, invalid, []);
    }
    // NCBI reference: c++/src/algo/blast/api/prelim_stage.cpp:177-188
    // (*thread)->Run(); (*thread)->Join(&result);
    // Repetition checks actual worker reclamation and bounded guest memory;
    // evidence is drained from the host after each call, not retained forever.
    const memorySizes = [];
    for (let i = 0; i < 24; i++) {
      const response = runPair(host, "tblastx", nuc, nuc, "6", ["-num_threads", "2"]);
      const events = await host.waitForWorkers();
      assert.equal(response.status, 0, response.error);
      assert.equal(events.filter(e=>e.event==="exited").length,1);
      assert.equal(host.events.length,0);
      memorySizes.push(host.memory.buffer.byteLength);
      save(`repeat-stress-${i}`, response, events);
    }
    assert.equal(new Set(memorySizes.slice(-12)).size,1,"guest linear memory kept growing on identical repeated searches");
  } finally { await host.close(); }
  return records;
}
module.exports = { runPair, checkReactor };
if (require.main === module) {
  checkReactor(...process.argv.slice(2)).then(rows => console.log(`reactor: ${rows.length} recorded runs passed`))
    .catch(error => { console.error(error); process.exitCode = 1; });
}
