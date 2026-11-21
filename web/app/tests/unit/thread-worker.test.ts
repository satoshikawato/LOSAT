// A thread worker (web/packages/losat/src/infra/engine-worker/thread-worker.ts) with a fake instance of the threaded
// module: it is ready only once its link to the main thread carries messages, and it reports a thread that stops on
// that link. Line WQA's sweep, D8 (2026-11-20): it opened a BroadcastChannel for those reports, which Firefox refuses
// with SecurityError where the site's data is blocked, so no thread worker became ready and a search of more than one
// thread could not run threaded. The browser here refuses BroadcastChannel. The real thread workers run in the browsers
// (tests/e2e/engine.spec.ts, tests/e2e/private-window.spec.ts).
import { afterEach, describe, expect, it, vi } from 'vitest';
import { InstanceUse } from '../../../packages/losat/src/infra/engine-worker/instance-use';
import {
  SLOT_FAILED,
  SLOT_PREPARING,
  SLOT_READY,
  SLOT_RUNNING,
  type ThreadCommand,
} from '../../../packages/losat/src/infra/engine-worker/protocol';

vi.mock('../../../packages/losat/src/infra/reactor/instance', () => ({
  instantiateThread: () =>
    Promise.resolve({
      start: () => {
        throw new Error('unreachable');
      },
    }),
}));

/** Waits until `condition` holds. */
async function until(condition: () => boolean): Promise<void> {
  for (let i = 0; i < 1000 && !condition(); i++) await new Promise((resolve) => setTimeout(resolve, 1));
  expect(condition()).toBe(true);
}

const ports: MessagePort[] = [];
afterEach(() => {
  for (const port of ports.splice(0)) port.close();
  vi.unstubAllGlobals();
});

describe("a thread worker where the site's data is blocked", () => {
  it('is ready once the answer to its hello came back on its link, and reports a thread that stops there', async () => {
    vi.stubGlobal(
      'BroadcastChannel',
      class {
        constructor() {
          throw new DOMException('The operation is insecure.', 'SecurityError');
        }
      },
    );
    const scope = { onmessage: null as ((event: MessageEvent<ThreadCommand>) => void) | null, postMessage: vi.fn() };
    vi.stubGlobal('self', scope);
    vi.resetModules();
    await import('../../../packages/losat/src/infra/engine-worker/thread-worker');

    const slot = new Int32Array(new SharedArrayBuffer(4));
    const use = InstanceUse.create();
    const link = new MessageChannel();
    ports.push(link.port1, link.port2);
    const received: unknown[] = [];
    link.port2.onmessage = (event: MessageEvent) => void received.push(event.data);
    const prepare: ThreadCommand = {
      type: 'prepare',
      module: {} as WebAssembly.Module,
      memory: {} as WebAssembly.Memory,
      slot: slot.buffer as SharedArrayBuffer,
      faults: link.port1,
      use: use.buffer,
    };
    scope.onmessage?.({ data: prepare } as MessageEvent<ThreadCommand>);
    await until(() => received.length === 1);
    expect(received).toEqual([{ type: 'hello' }]);
    // Its instance is ready, but no answer came: the link may not carry a fault yet (link.ts).
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(Atomics.load(slot, 0)).toBe(SLOT_PREPARING);
    link.port2.postMessage({ type: 'hello' });
    await until(() => Atomics.load(slot, 0) === SLOT_READY);

    // As ThreadHost.spawn claims the slot, counts the thread and starts it.
    Atomics.store(slot, 0, SLOT_RUNNING);
    use.startThread();
    const started = new Int32Array(new SharedArrayBuffer(4));
    scope.onmessage?.({ data: { type: 'start', tid: 3, startArg: 0, started: started.buffer as SharedArrayBuffer } } as MessageEvent<ThreadCommand>);
    await until(() => received.length === 2);
    expect(received[1]).toEqual({ type: 'fault', tid: 3, message: 'The engine stopped: unreachable' });
    expect(Atomics.load(slot, 0)).toBe(SLOT_FAILED);
    expect(use.uses()).toBe(0);
    expect(scope.postMessage).not.toHaveBeenCalled();
  });
});
