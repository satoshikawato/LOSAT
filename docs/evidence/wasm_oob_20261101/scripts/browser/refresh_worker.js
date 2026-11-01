onmessage = ({ data: { module, memory } }) => {
  let outcome;
  try { outcome = String(new WebAssembly.Instance(module, { env: { memory } }).exports.probe(0)); } catch (err) { outcome = err.message; }
  postMessage(outcome);
};
postMessage("ready");
