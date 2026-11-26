"use strict";
const { createThreadHost } = require("./wasi_thread_host");
const { endOnAbortTrap } = require("./wasi_artifact");
// NCBI reference (598d8ae6): c++/src/app/blast/blastn_app.cpp:172-176
(async () => {
  const wasm = process.argv[2];
  if (!wasm) throw new Error("usage: node run_losat_wasi_threads.js <threaded-command.wasm> [args...]");
  const host = await createThreadHost(wasm, "threaded-command", process.argv.slice(3));
  try {
    process.exitCode = host.start();
    await host.waitForWorkers();
  } catch (err) {
    // The module's abort on the main thread ends the process at once by SIGABRT, as the native
    // command's abort ends every thread (wasi_artifact.js); the workers are not waited for.
    endOnAbortTrap(err);
    throw err;
  } finally { await host.close(); }
  // The process ends when its event loop is empty, with process.exitCode, and never with
  // process.exit(): Node's exit handler shuts down V8's worker threads before it frees the main
  // isolate, and a background compile there that waits for a garbage collection on this thread
  // then waits forever (one 4-thread tblastn command of a local gate, 2026-10-29: the main
  // thread joined a node-V8Worker in ConcurrentBaselineCompiler -> CollectionBarrier::
  // AwaitCollectionBackground for 3 hours). Freeing the isolate first ends such a wait.
})().catch(err => { console.error(err); process.exitCode = 1; });
