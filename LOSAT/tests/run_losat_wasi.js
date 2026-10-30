#!/usr/bin/env node

"use strict";

const fs = require("fs");
const { WASI } = require("wasi");
const { prepareArtifact, unsignedWasiImports } = require("./wasi_artifact");
const { sharedFileTable, sharedFileImports } = require("./wasi_files");

// NCBI reference:
// ncbi-blast/c++/src/algo/blast/blastinput/cmdline_flags.cpp:46-75
// The command-Wasm certification passes the same query, subject, output,
// genetic-code, and thread arguments as native LOSAT. This runner only supplies
// WASI preview1; it does not rewrite output or implement BLAST behavior.
const wasmPath = process.argv[2];
if (!wasmPath) {
  console.error("usage: node tests/run_losat_wasi.js <losat.wasm> [args...]");
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
  env: { ...process.env, PWD: process.cwd() },
  // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729
// ```c++
//             fstrm->open(AsString().c_str(),IOS_BASE::in | mode);
//             if ( !fstrm->is_open() ) {
//                 delete fstrm;
// ```
// Preserve lexical filenames at the recorded host cwd.
// Root preopen permits absolute paths; PWD initializes the guest cwd for
// relative filenames through std::env::set_current_dir, without rewriting argv.
    preopens: { "/": "/" },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  // NCBI reference: c++/src/app/blast/blastn_app.cpp:172-176
  // CATCH_ALL(status); return status;
  // Module reuse is host preparation only; preserve validation and exit status.
  const { module } = prepareArtifact(bytes, "serial-command");
  // NCBI reference (598d8ae6): c++/src/corelib/ncbiargs.cpp:727-729 (see createWasi above)
  // The files the command opens from "/" are opened by name as the OS opens them, through the
  // table of files that the threaded host uses too (wasi_files.js, the comment above
  // SHARED_FD_BASE): uvwasi cannot open /dev/stdin or /dev/stdout when they are pipes. The module
  // exports its memory, so the table reads it through the instance.
  let memory = null;
  const guestMemory = { get buffer() { return memory.buffer; }, grow: (pages) => memory.grow(pages) };
  const local = unsignedWasiImports(wasi.wasiImport);
  const instance = await WebAssembly.instantiate(module, {
    // The serial memory has no declared maximum, so it can pass 2 GiB.
    wasi_snapshot_preview1: { ...local, ...sharedFileImports(sharedFileTable(), guestMemory, local) },
  });
  memory = instance.exports.memory;
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
