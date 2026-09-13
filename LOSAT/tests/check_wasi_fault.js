"use strict";
const fs = require("node:fs");
const { createThreadHost } = require("./wasi_thread_host");
const { runPair } = require("./check_wasi_reactor");
// NCBI reference: c++/src/app/blast/blast_app_util.hpp:260-266
// LOG_POST(Error << "Error: " << e.what()); exit_code = BLAST_UNKNOWN_ERROR;
// The parent integration process requires explicit failure, never a hung guest.
(async () => {
  const [artifact, query, fault] = process.argv.slice(2);
  const host = await createThreadHost(artifact, "threaded-reactor", [], { fault, startTimeoutMs: 2000 });
  try {
    const seq = fs.readFileSync(query, "utf8");
    const result = runPair(host, "blastn", seq, seq, "6", ["-task", "blastn", "-num_threads", "2"]);
    await host.waitForWorkers(2000);
    throw new Error(`injected ${fault} unexpectedly returned ${result.status}`);
  } finally { await host.close(); }
})().catch(error => { console.error(error); process.exitCode = 1; });
