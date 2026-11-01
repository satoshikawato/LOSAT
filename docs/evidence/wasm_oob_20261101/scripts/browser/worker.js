onmessage = ({ data: { module, memory, kind } }) => {
  let outcome = "ok";
  try {
    const instance = new WebAssembly.Instance(module, { env: { memory } });
    instance.exports.probe(kind);
  } catch (err) { outcome = err.message; }
  postMessage(outcome);
};
postMessage("ready");
