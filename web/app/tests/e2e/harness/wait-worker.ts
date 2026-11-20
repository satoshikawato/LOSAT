// A worker of worker-end.spec.ts, as the Engine worker of a threaded runtime: a shared memory of
// its own, of the given maximum, with the application's keeper (memory-keeper-worker.ts) linked to
// the page (link.ts), and an instance of the test's module on it. Once the keeper is ready it says
// so, then waits on the memory without a timeout (the module's export `park`) until the page
// ends it.
import type { KeeperCommand } from '../../../../packages/losat/src/infra/engine-worker/protocol';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<{ module: WebAssembly.Module; maximumPages: number }>) => void) | null;
  postMessage(message: string | { readonly type: 'link'; readonly port: MessagePort }, transfer?: Transferable[]): void;
};

scope.onmessage = (event) => {
  const { module, maximumPages } = event.data;
  let memory: WebAssembly.Memory;
  try {
    memory = new WebAssembly.Memory({ initial: 1, maximum: maximumPages, shared: true });
  } catch (error) {
    scope.postMessage(`no memory: ${String(error)}`);
    return;
  }
  const keeper = new Worker(new URL('../../../../packages/losat/src/infra/engine-worker/memory-keeper-worker.ts', import.meta.url), { type: 'module' });
  keeper.onmessage = () => {
    const instance = new WebAssembly.Instance(module, {
      env: { memory },
      wasi_snapshot_preview1: { sched_yield: () => 0 },
    });
    scope.postMessage('waiting');
    (instance.exports['park'] as () => void)();
  };
  // As the Engine worker links the keeper to the page (engine-worker.ts).
  const link = new MessageChannel();
  scope.postMessage({ type: 'link', port: link.port2 }, [link.port2]);
  const keep: KeeperCommand = { type: 'keep', memory, maximumPages, port: link.port1 };
  keeper.postMessage(keep, [link.port1]);
};
