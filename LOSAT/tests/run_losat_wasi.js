#!/usr/bin/env node

"use strict";

const fs = require("fs");
const { WASI } = require("wasi");
const { inspectArtifact } = require("./wasi_artifact");

// NCBI reference:
// ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
// The command-Wasm certification passes the same query, subject, output,
// genetic-code, and thread arguments as native LOSAT. This runner only supplies
// WASI preview1; it does not rewrite output or implement BLAST behavior.
const wasmPath = process.argv[2];
if (!wasmPath) {
  console.error("usage: node tests/run_losat_wasi.js <LOSAT.wasm> [args...]");
  process.exit(2);
}

const args = [wasmPath, ...process.argv.slice(3)];
const wasi = new WASI({
  version: "preview1",
  // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
  // CATCH_ALL(status)
  // return status;
  returnOnExit: true,
  args,
  env: process.env,
  preopens: { "/": "/" },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
  // CATCH_ALL(status); return status;
  inspectArtifact(bytes, "serial-command");
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  if (typeof instance.exports._start !== "function") {
    throw new Error(
      `${wasmPath} does not export _start; use the LOSAT command-Wasm artifact.`,
    );
  }
  // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
  // CATCH_ALL(status)
  // return status;
  process.exitCode = wasi.start(instance);
})().catch((error) => {
  console.error(error && error.stack ? error.stack : String(error));
  process.exit(1);
});
