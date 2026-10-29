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
  // Every thread worker has exited, the files are closed and the command's output was written
  // synchronously (WASI fd_write): exit with the command's status now. Node's own end of the
  // process (freeing the environment and the main isolate after the 'exit' event) took 5-10 ms
  // more on the workstation and does nothing for the command. A failure above still ends the
  // usual way, after its message.
  process.exit();
})().catch(err => { console.error(err); process.exitCode = 1; });
