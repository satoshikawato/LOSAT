// The order in which the Engine worker (web/packages/losat/src/infra/engine-worker/engine-worker.ts) answers `init` and
// takes the searches that come with it, with a fake runtime and fake Web Locks. Line WQA's sweep, D2, its review, H1 and
// M1 (2026-11-19): the worker said that it started only once the browser granted its lock, and a search posted together
// with `init` (the first search of a new worker: after a cancel, before the warm-up, after a renewal) reached its
// synchronous call into the engine before that. The page therefore did not wait for the lock during that search, so a
// worker ended then left the run "Searching" for good, and an error then read as a worker that could not start. The
// worker's searches run in the browsers in tests/e2e/unavailable.spec.ts and tests/e2e/engine.spec.ts.
import { afterEach, describe, expect, it, vi } from 'vitest';
import type { LifeLocks } from 'losat';
import type { EngineCommand, EngineEvent, EngineInit, WorkerPrepare, WorkerRun } from '../../../packages/losat/src/infra/engine-worker/protocol';

/** What the worker did, in order: the type of each message that it posted, and each call into its runtime. */
const done = vi.hoisted(() => [] as string[]);
/** The options of the latest runtime that the worker made. */
const made = vi.hoisted(() => ({ options: undefined as unknown }));

vi.mock('../../../packages/losat/src/infra/engine-worker/runtime', () => ({
  EngineRuntime: class {
    constructor(options: unknown) {
      made.options = options;
    }
    run(request: { readonly runId: string }): Promise<unknown> {
      done.push(`search ${request.runId}`);
      return Promise.resolve({ path: 'serial', threads: 1, engineBuild: 'serial', memory: {}, subjectRetained: false });
    }
    prepare(): Promise<void> {
      done.push('prepare');
      return Promise.resolve();
    }
  },
}));

/** Web Locks whose requests are granted only when the test says so, as a browser grants them a moment later. */
class LateLocks implements LifeLocks {
  private readonly waiting: Array<() => void> = [];

  request(_name: string, _options: { readonly signal?: AbortSignal }, callback: () => unknown): Promise<unknown> {
    return new Promise((resolve) => this.waiting.push(() => resolve(callback())));
  }

  grant(): void {
    for (const grant of this.waiting.splice(0)) grant();
  }
}

const INIT: EngineInit = {
  type: 'init',
  generation: 1,
  modules: { serial: { asset: { url: '/engine/serial.wasm', sha256: '' } }, threads: { asset: { url: '/engine/threads.wasm', sha256: '' } } },
  builds: { serial: 'serial', threads: 'threads' },
  lifeLock: 'losat-web:engine-life:1',
} as unknown as EngineInit;

const PREPARE: WorkerPrepare = { type: 'prepare', generation: 1, threads: 1 };

function run(runId: string): WorkerRun {
  return { type: 'run', generation: 1, runId, output: { close: () => undefined } } as unknown as WorkerRun;
}

/** Starts the worker's module with `locks` as its Web Locks; returns how to post a command to it. */
async function startWorker(
  locks: LifeLocks | undefined,
): Promise<{ post: (command: EngineCommand) => void; events: EngineEvent[]; transfers: Array<Transferable[] | undefined> }> {
  const events: EngineEvent[] = [];
  const transfers: Array<Transferable[] | undefined> = [];
  const scope = {
    onmessage: null as ((event: MessageEvent<EngineCommand>) => void) | null,
    postMessage(event: EngineEvent, transfer?: Transferable[]) {
      events.push(event);
      transfers.push(transfer);
      done.push(event.type);
    },
  };
  vi.stubGlobal('self', scope);
  vi.stubGlobal('navigator', { locks });
  vi.resetModules();
  await import('../../../packages/losat/src/infra/engine-worker/engine-worker');
  return { post: (command) => scope.onmessage?.({ data: command } as MessageEvent<EngineCommand>), events, transfers };
}

/** Lets every promise and timer that is due settle. */
async function settle(): Promise<void> {
  for (let i = 0; i < 5; i++) await new Promise((resolve) => setTimeout(resolve, 0));
}

afterEach(() => {
  done.length = 0;
  vi.unstubAllGlobals();
});

describe('Engine worker', () => {
  it('says at once that it started, and starts no search before it holds its lock and has said so', async () => {
    const locks = new LateLocks();
    const { post, events } = await startWorker(locks);
    post(INIT);
    post(run('r1'));
    await settle();
    // Before, the search ran here, and the worker said nothing until the browser granted its lock.
    expect(done).toEqual(['started']);
    locks.grant();
    await settle();
    expect(done).toEqual(['started', 'life', 'search r1', 'done']);
    expect(events.find((event) => event.type === 'life')).toEqual({ type: 'life', generation: 1, held: true });
  });

  it('starts no preparation before it holds its lock and has said so', async () => {
    const locks = new LateLocks();
    const { post } = await startWorker(locks);
    post(INIT);
    post(PREPARE);
    await settle();
    expect(done).toEqual(['started']);
    locks.grant();
    await settle();
    expect(done).toEqual(['started', 'life', 'prepare']);
  });

  it('without Web Locks, says that it holds no lock, and then searches', async () => {
    const { post, events } = await startWorker(undefined);
    post(INIT);
    post(run('r1'));
    await settle();
    expect(done).toEqual(['started', 'life', 'search r1', 'done']);
    expect(events.find((event) => event.type === 'life')).toEqual({ type: 'life', generation: 1, held: false });
  });

  // Line WQA's sweep, D8 (2026-11-20): the thread workers and the keeper reached the page on BroadcastChannels, which
  // Firefox refuses where the site's data is blocked. The worker makes their links and posts the page's end of each.
  it("posts the page's end of each link that its runtime makes, and the faults of its thread host", async () => {
    vi.stubGlobal(
      'BroadcastChannel',
      class {
        constructor() {
          throw new DOMException('The operation is insecure.', 'SecurityError');
        }
      },
    );
    const { post, events, transfers } = await startWorker(undefined);
    post({ ...INIT, generation: 3 });
    const { main } = made.options as { main: { link(role: 'thread' | 'keeper'): MessagePort; fault(fault: unknown): void } };
    const keeper = main.link('keeper');
    const linked = events.at(-1) as { type: string; generation: number; role: string; port: MessagePort };
    expect(linked).toMatchObject({ type: 'link', generation: 3, role: 'keeper' });
    expect(transfers.at(-1)).toEqual([linked.port]);
    // The two ends are one channel.
    const heard = new Promise((resolve) => {
      linked.port.onmessage = (event: MessageEvent) => resolve(event.data);
    });
    keeper.postMessage({ type: 'hello' });
    expect(await heard).toEqual({ type: 'hello' });
    keeper.close();
    linked.port.close();
    main.fault({ type: 'fault', tid: 2, message: 'The engine stopped: thread 2 did not start in time' });
    expect(events.at(-1)).toEqual({ type: 'fault', generation: 3, tid: 2, message: 'The engine stopped: thread 2 did not start in time' });
  });
});
