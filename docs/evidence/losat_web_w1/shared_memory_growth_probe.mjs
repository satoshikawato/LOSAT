// S09 task 6 (plan §8): does a bulk memory operation (memory.fill / memory.copy) in one
// browser worker trap after another worker grows the shared memory, as V8 can in Node
// (LOSAT/tests/wasi_shared_memory.js guards against it there)? The module is the consumer
// fixture of LOSAT/tests/test_wasi_shared_memory.js: it signals a grower worker, waits
// inside Wasm until the memory has grown, then fills or copies in the new page. Every
// operation runs TRIALS times, unguarded and with the guard of wasi_shared_memory.js.
// Usage (repository root): node docs/evidence/losat_web_w1/shared_memory_growth_probe.mjs
//   TRIALS=20 (default); LOSAT_WEB_WEBKIT_EXECUTABLE=<launcher> for WebKit on a host whose
//   libraries Playwright cannot install.
import http from 'node:http';
import { createRequire } from 'node:module';
const repo = new URL('../../../', import.meta.url);
const require = createRequire(new URL('web/app/package.json', repo));
const { chromium, firefox, webkit } = require('@playwright/test');
const { guardSharedMemory } = require(new URL('LOSAT/tests/wasi_shared_memory.js', repo).pathname);
const TRIALS = Number(process.env.TRIALS || 20);
const H = { 'Cross-Origin-Opener-Policy': 'same-origin', 'Cross-Origin-Embedder-Policy': 'require-corp', 'Cross-Origin-Resource-Policy': 'same-origin' };
const fixture = `
const uleb = (value) => { const out = []; do { const b = value & 127; value >>>= 7; out.push(b | (value ? 128 : 0)); } while (value); return out; };
const vector = (items) => [...uleb(items.length), ...items.flat()];
const name = (s) => { const b = [...new TextEncoder().encode(s)]; return [...uleb(b.length), ...b]; };
const section = (id, bytes) => [id, ...uleb(bytes.length), ...bytes];
const body = (ops) => { const b = [0, ...ops, 0x0b]; return [...uleb(b.length), ...b]; };
const args = [0x20, 0, 0x20, 1, 0x20, 2];
function fixture(extra) {
  return new Uint8Array([0, 97, 115, 109, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 3, 0x7f, 0x7f, 0x7f, 0]),
    ...section(2, vector([[...name("env"), ...name("memory"), 2, 3, 1, 4]])),
    ...section(3, vector(extra.map(() => [0]))),
    ...section(7, vector(extra.map((_, i) => [...name("extra" + i), 0, i]))),
    ...section(10, vector(extra.map(body))),
  ]);
}
function consumer(operation) {
  return [0x41, 0, 0x41, 1, 0xfe, 0x17, 2, 0,
    0x03, 0x40, 0x41, 0, 0xfe, 0x10, 2, 0, 0x41, 2, 0x47, 0x0d, 0, 0x0b,
    ...args, 0xfc, ...(operation === "fill" ? [11, 0] : [10, 0, 0])];
}`;
const files = {
  '/index.html': ['text/html', '<!doctype html><title>g</title><script type="module" src="/main.js"></script>'],
  '/main.js': ['text/javascript', `${fixture}
    async function one(operation, parameters) {
      const memory = new WebAssembly.Memory({ initial: 1, maximum: 4, shared: true });
      new Uint8Array(memory.buffer, 83, 8).fill(83);
      let bytes = fixture([consumer(operation)]);
      if (window.guarded) bytes = new Uint8Array(await (await fetch('/guard', { method: 'POST', body: bytes })).arrayBuffer());
      const module = await WebAssembly.compile(bytes);
      const grower = new Worker('/grower.js');
      const runner = new Worker('/consumer.js');
      grower.postMessage({ memory });
      const done = new Promise((resolve) => { runner.onmessage = (e) => resolve(e.data); runner.onerror = (e) => resolve({ error: 'worker error ' + e.message }); });
      runner.postMessage({ memory, module, parameters });
      const result = await Promise.race([done, new Promise((r) => setTimeout(() => r({ error: 'timeout' }), 10000))]);
      grower.terminate(); runner.terminate();
      return { operation, ...result, byteLength: memory.buffer.byteLength };
    }
    window.run = async (TRIALS) => {
      const out = [];
      for (const [operation, parameters] of [['fill', [65536, 83, 8]], ['copy destination', [65536, 83, 8]], ['copy source', [64, 65536, 8]]]) {
        for (let i = 0; i < TRIALS; i++) out.push(await one(operation.startsWith('copy') ? 'copy' : 'fill', parameters).then((r) => ({ ...r, operation })));
      }
      return out;
    };`],
  '/grower.js': ['text/javascript', `onmessage = (e) => { const m = e.data.memory; const view = new Int32Array(m.buffer);
      while (Atomics.load(view, 0) !== 1) Atomics.wait(view, 0, 0, 10);
      m.grow(1); new Uint8Array(m.buffer, 65536, 8).fill(83); Atomics.store(view, 0, 2); };`],
  '/consumer.js': ['text/javascript', `onmessage = async (e) => {
      const { memory, module, parameters } = e.data;
      const instance = await WebAssembly.instantiate(module, { env: { memory } });
      try { instance.exports.extra0(...parameters);
        postMessage({ ok: true, bytes: [...new Uint8Array(memory.buffer, parameters[0], 8)] });
      } catch (error) { postMessage({ ok: false, trap: String(error) }); } };`],
};
const server = http.createServer((req, res) => {
  if (req.url === '/guard') { const chunks = []; req.on('data', (c) => chunks.push(c)); req.on('end', () => { const g = guardSharedMemory(Buffer.concat(chunks)); res.writeHead(200, { ...H, 'Content-Type': 'application/wasm' }); res.end(g.bytes); }); return; }
  const f = files[req.url.split('?')[0]]; if (!f) { res.writeHead(404); res.end(); return; } res.writeHead(200, { ...H, 'Content-Type': f[0] }); res.end(f[1]); }).listen(0);
const port = server.address().port;
const webkitExecutable = process.env.LOSAT_WEB_WEBKIT_EXECUTABLE;
for (const [label, type, opts] of [['chromium', chromium, {}], ['firefox', firefox, {}], ['webkit', webkit, webkitExecutable ? { executablePath: webkitExecutable } : {}]]) {
  const browser = await type.launch(opts);
  const page = await browser.newPage();
  page.on('pageerror', (e) => console.log(label, 'pageerror', e.message));
  await page.goto(`http://localhost:${port}/index.html`);
  for (const guarded of [false, true]) {
    const r = await page.evaluate(([g, n]) => { window.guarded = g; return window.run(n); }, [guarded, TRIALS]);
    const summary = {};
    for (const x of r) { const k = `${x.operation}: ${x.ok ? (x.bytes.every((b) => b === 83) ? 'ok' : 'WRONG BYTES') : (x.trap || x.error)}`; summary[k] = (summary[k] || 0) + 1; }
    console.log(label, browser.version(), guarded ? 'guarded' : 'unguarded', JSON.stringify(summary));
  }
  await browser.close();
}
server.close();
