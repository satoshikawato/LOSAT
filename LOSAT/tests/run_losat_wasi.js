#!/usr/bin/env node

"use strict";

const fs = require("fs");
const { WASI } = require("wasi");

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
  args,
  env: process.env,
  preopens: { "/": "/" },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  if (typeof instance.exports._start !== "function") {
    throw new Error(
      `${wasmPath} does not export _start; use the LOSAT command-Wasm artifact.`,
    );
  }
  wasi.start(instance);
})().catch((error) => {
  console.error(error && error.stack ? error.stack : String(error));
  process.exit(1);
});
