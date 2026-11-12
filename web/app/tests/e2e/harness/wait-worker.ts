// A worker of worker-end.spec.ts: a shared memory of its own, of the given maximum, and an
// instance of the test's module on it. It says so, then waits on the memory without a timeout
// (the module's export `park`) until the page ends it.
const scope = self as unknown as {
  onmessage: ((event: MessageEvent<{ module: WebAssembly.Module; maximumPages: number }>) => void) | null;
  postMessage(message: string): void;
};

scope.onmessage = (event) => {
  let memory: WebAssembly.Memory;
  try {
    memory = new WebAssembly.Memory({ initial: 1, maximum: event.data.maximumPages, shared: true });
  } catch (error) {
    scope.postMessage(`no memory: ${String(error)}`);
    return;
  }
  const instance = new WebAssembly.Instance(event.data.module, {
    env: { memory },
    wasi_snapshot_preview1: { sched_yield: () => 0 },
  });
  scope.postMessage('waiting');
  (instance.exports['park'] as () => void)();
};
