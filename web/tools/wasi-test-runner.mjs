// Cargo runner for wasm32-wasip1 test binaries, so engine tests that compile only for
// wasm32 (for example `LOSAT/src/web_api.rs`) can run in CI:
//   CARGO_TARGET_WASM32_WASIP1_RUNNER="node web/tools/wasi-test-runner.mjs" \
//     cargo test --lib --target wasm32-wasip1 --no-default-features -- web_api::tests
import { readFile } from 'node:fs/promises';
import { createRequire } from 'node:module';
import { WASI } from 'node:wasi';
import { argv, env, exit } from 'node:process';

// The test binary's output to a full standard output or error pipe waits as the native write
// blocks (blockingStdioWrite, LOSAT/tests/wasi_artifact.js), instead of ending the run with EAGAIN.
const { blockingStdioWrite } = createRequire(import.meta.url)('../../LOSAT/tests/wasi_artifact.js');

const [wasmPath, ...args] = argv.slice(2);
const wasi = new WASI({
  version: 'preview1',
  args: [wasmPath, ...args],
  env,
  preopens: { '/': '/' },
  returnOnExit: true,
});
const module = await WebAssembly.compile(await readFile(wasmPath));
const { wasi_snapshot_preview1: wasiImport } = wasi.getImportObject();
const instance = await WebAssembly.instantiate(module, {
  wasi_snapshot_preview1: { ...wasiImport, fd_write: blockingStdioWrite(wasiImport.fd_write) },
});
exit(wasi.start(instance));
