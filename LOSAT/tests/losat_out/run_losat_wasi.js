const { WASI } = require('wasi');
const fs = require('fs');

const wasmPath = process.argv[2];
const args = process.argv.slice(3);
const wasi = new WASI({
  version: 'preview1',
  args: [wasmPath, ...args],
  env: process.env,
  preopens: { '/': '/' },
});

(async () => {
  const bytes = fs.readFileSync(wasmPath);
  const module = await WebAssembly.compile(bytes);
  const instance = await WebAssembly.instantiate(module, {
    wasi_snapshot_preview1: wasi.wasiImport,
  });
  if (typeof instance.exports._start !== 'function') {
    console.error(`${wasmPath} does not export _start; use the LOSAT command Wasm artifact.`);
    process.exit(1);
  }
  wasi.start(instance);
})().catch((err) => {
  console.error(err);
  process.exit(1);
});
