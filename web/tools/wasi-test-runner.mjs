// Cargo runner for wasm32-wasip1 test binaries, so engine tests that compile only for
// wasm32 (for example `LOSAT/src/web_api.rs`) can run in CI:
//   CARGO_TARGET_WASM32_WASIP1_RUNNER="node web/tools/wasi-test-runner.mjs" \
//     cargo test --lib --target wasm32-wasip1 --no-default-features -- web_api::tests
import { readFile } from 'node:fs/promises';
import { WASI } from 'node:wasi';
import { argv, env, exit } from 'node:process';

const [wasmPath, ...args] = argv.slice(2);
const wasi = new WASI({
  version: 'preview1',
  args: [wasmPath, ...args],
  env,
  preopens: { '/': '/' },
  returnOnExit: true,
});
const module = await WebAssembly.compile(await readFile(wasmPath));
const instance = await WebAssembly.instantiate(module, wasi.getImportObject());
exit(wasi.start(instance));
