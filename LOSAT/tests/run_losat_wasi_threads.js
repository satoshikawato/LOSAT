"use strict";
const { createThreadHost } = require("./wasi_thread_host");
// NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
// CATCH_ALL(status); return status;
(async () => {
  const wasm = process.argv[2];
  if (!wasm) throw new Error("usage: node run_losat_wasi_threads.js <threaded-command.wasm> [args...]");
  const host = await createThreadHost(wasm, "threaded-command", process.argv.slice(3));
  try {
    process.exitCode = host.start();
    await host.waitForWorkers();
  } finally { await host.close(); }
})().catch(err => { console.error(err); process.exitCode = 1; });
