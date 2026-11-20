onmessage = ({ data: { module, memory } }) => {
  const view = new Int32Array(memory.buffer);
  const deadline = performance.now() + 25000;
  while (Atomics.load(view, 0) === 0 && performance.now() < deadline) Atomics.wait(view, 0, 0, 100);
  if (Atomics.load(view, 0) === 0) return;
  const old = new WebAssembly.Instance(module, { env: { memory } }).exports.grow(1);
  Atomics.store(view, 1, old * 65536 + 64);
};
postMessage("ready");
