# WebAssembly and Browser

LOSAT runs as WebAssembly in two forms, built from the same code as the native `losat`:

1. **WASI command modules (`losat.wasm`).** The command line for WASI runtimes such as Node.js, in a serial (`wasm32-wasip1`) and a threaded (`wasm32-wasip1-threads`) build. Both are release files on [GitHub Releases](https://github.com/satoshikawato/LOSAT/releases).
2. **The npm package `losat`, for web pages.** It runs `blastn`, `blastp`, `blastx`, `tblastn` and `tblastx` in the user's browser, so sequences never leave the user's machine. LOSAT Web runs its searches through this package, and the browser version of [gbdraw](https://github.com/satoshikawato/gbdraw) uses the same WebAssembly engine.

## WASI command modules (`losat.wasm`)

The modules take the same arguments as the native `losat`. The repository's Node.js runners in `LOSAT/tests` give a module what it needs from its host (files, the working directory, the time zone, and memory past 2 GiB):

```bash
# serial module
node LOSAT/tests/run_losat_wasi.js losat.wasm blastn -query query.fna -subject subject.fna -outfmt 6

# threaded module
node LOSAT/tests/run_losat_wasi_threads.js losat-threads.wasm blastn -query query.fna -subject subject.fna -num_threads 4 -outfmt 6
```

[Platforms and WebAssembly](../platforms.md#wasi-losatwasm) lists what a WASI host must pass, the limits of WASI (such as compressed inputs), and other hosts. [Threads, speed and memory](../performance.md#webassembly) compares their speed with the native binary.

## The npm package `losat`

The package holds two WebAssembly builds of the engine and the TypeScript host that runs them:

- `losat-web-serial.wasm` runs one thread and works in any browser and in Node.js.
- `losat-web-threads.wasm` runs a search on several threads. It needs a cross-origin isolated page.

A search takes an argv, as on the command line, and FASTA sequences that the page gives it. It writes `-outfmt 0`, `6` and `7` byte for byte as the command line does, and an HSP record for each hit. For `-db` and BLASTN `-filtering_db`, the WebAssembly engine and the npm package read BLAST databases from the run's `files` (`EngineRunRequest.files`), while LOSAT Web's browser form refuses both.

`engine.runTool(argv, options)` runs `windowmasker`, `blast_formatter`, `makeblastdb`, `dustmasker`, `segmasker` and `convert2blastmask` with the command line's words on the files of `options.directory`, and gives the exit status, standard output, standard error and the files that the command line writes (for `makeblastdb`, the database files that it leaves in the run's directory, dated in the browser's time zone as the command line dates them in the user's).

### Getting the package

The package is not on npm. Build it from the repository:

1. Build the two WebAssembly modules: `python3 web/adapter/tools/build_reactors.py --target-dir <build directory> --output-dir <modules directory>`.
2. Install the web workspace's dependencies with `npm ci` in `web/`.
3. In `web/packages/losat`, run `LOSAT_WEB_REACTORS=<modules directory> npm run build`. The package is written to `dist/`; `npm pack` there makes a file you can install.

The build needs Node.js 22.18 or later, Python 3.11 or later, and Rust 1.92 with the `wasm32-wasip1` and `wasm32-wasip1-threads` targets.

### What a page must do

- **Serve the modules** as `application/wasm`. The host checks their SHA-256 before it runs them.
- **Isolate the page for threads.** The threaded module needs `SharedArrayBuffer`, which a browser gives only to a cross-origin isolated page. Serve the page and its scripts with:

  ```text
  Cross-Origin-Opener-Policy: same-origin
  Cross-Origin-Embedder-Policy: require-corp
  ```

  On a page that is not isolated, a search that asks for more than one thread runs on the serial module with one thread, and `RuntimeInfo.fallbackReason` says why.
- **Drop the outputs that the engine discards.** When an output cannot be written, the engine tells the host to drop what arrived of it (stream 9). The outputs are the command line's bytes only for a host that does this. `RunOutputReceiver` reports each such signal through its `onDiscard` callback.

The package's [README](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/web/packages/losat/README.md) shows how to start `WasmEngine`, how the outputs arrive, how memory is reserved, and how to call the serial module directly. The [ABI v2 specification](https://github.com/satoshikawato/LOSAT/blob/v0.2.0/docs/web/abi_v2.md) describes the interface between the host and the modules.
