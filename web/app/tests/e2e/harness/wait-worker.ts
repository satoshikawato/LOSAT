// A worker of worker-end.spec.ts, as the Engine worker of a threaded runtime: a shared memory of
// its own, of the given maximum, with the application's keeper (memory-keeper-worker.ts) on the
// given channel, and an instance of the test's module on it. Once the keeper is ready it says
// so, then waits on the memory without a timeout (the module's export `park`) until the page
// ends it.
import type { KeeperCommand } from '../../../../packages/losat/src/infra/engine-worker/protocol';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<{ module: WebAssembly.Module; maximumPages: number; keeperChannel: string }>) => void) | null;
  postMessage(message: string): void;
};

scope.onmessage = (event) => {
  const { module, maximumPages, keeperChannel } = event.data;
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
  const keep: KeeperCommand = { type: 'keep', memory, maximumPages, channel: keeperChannel };
  keeper.postMessage(keep);
};
