"use strict";
// NCBI reference: c++/src/objtools/align_format/tabular.cpp:1100-1108
// x_PrintField(*iter); m_Ostream << "\n";
// Real serial WASI calls verify evidence retention. NCBI is a test oracle only.
const fs = require("node:fs");
const os = require("node:os");
const path = require("node:path");
const crypto = require("node:crypto");
const assert = require("node:assert/strict");
const { spawnSync } = require("node:child_process");
const artifact = path.resolve(process.argv[2]);
const directory = fs.mkdtempSync(path.join(os.tmpdir(), "losat-reuse-evidence-test-"));
const fixture = path.join(__dirname, "fasta/small_test.fasta");
const args = ["-query", fixture, "-subject", fixture, "-task", "blastn", "-outfmt", "6", "-num_threads", "1"];
const oracle = path.join(directory, "oracle.out");
const result = spawnSync(process.env.BLASTN_BIN || "blastn", [...args, "-out", oracle], {timeout:30000});
assert.equal(result.status, 0, result.stderr?.toString());
const expected = crypto.createHash("sha256").update(fs.readFileSync(oracle)).digest("hex");
const runner = path.join(__dirname, "benchmark_wasi_reuse.js");
const env = {...process.env, NODE_NO_WARNINGS:"1"};
delete env.NODE_OPTIONS;
for (const scenario of ["success", "mismatch", "bad-command", "stale-output"]) {
  const output = path.join(directory, `${scenario}.out`);
  const jobsFile = path.join(directory, `${scenario}-jobs.json`);
  const resultFile = path.join(directory, `${scenario}-result.json`);
  const jobs = [{case_id:"small", threads:1, repeat:0, timed:true, output,
    argv:["blastn", ...args, ...(scenario === "bad-command" ? ["--invalid-option"] : []), "-out", output],
    expected_sha256:scenario === "mismatch" ? "0".repeat(64) : expected}];
  fs.writeFileSync(jobsFile, JSON.stringify(jobs));
  if (scenario === "stale-output") fs.writeFileSync(output, "old result");
  const child = spawnSync(process.execPath, [runner, artifact, "serial-command", "compiled-module", jobsFile, "-out", resultFile], {env, timeout:30000});
  if (scenario === "stale-output") {
    assert.notEqual(child.status, 0);
    assert.equal(fs.readFileSync(output, "utf8"), "old result");
    continue;
  }
  const sample = JSON.parse(fs.readFileSync(resultFile)).samples[0];
  assert.equal(sample.status, scenario === "success" ? "PASS" : scenario === "mismatch" ? "PARITY_FAIL" : "FAIL");
  assert.equal(sample.eligible_for_timing, scenario === "success");
  assert.equal(child.status === 0, scenario === "success", child.stderr?.toString());
  assert.ok(sample.close_seconds >= 0);
  if (scenario === "success") assert.ok(sample.raw_equal);
}
console.log(`PASS: reuse success, raw mismatch, bad command, stale output; evidence ${directory}`);
