// The Wasm EngineGateway's runtime generations (web/packages/losat/src/infra/engine-worker/gateway.ts) with a
// fake Engine worker: runs, phases, cancel, renewal, faults, errors, the use count of the
// threaded instance and the links of the thread workers and the memory keeper (link.ts). The real Engine worker runs
// in the browsers (tests/e2e/engine.spec.ts).
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ENGINE_ENDED, ENGINE_NOT_STARTED, WasmEngine, type EngineAssets, type LifeLocks, type WasmEngineOptions } from 'losat';
import { InputMismatchError, RunCancelledError, type EngineRunRequest, type RuntimeInfo } from 'losat/ports';
import { InstanceUse } from '../../../packages/losat/src/infra/engine-worker/instance-use';
import {
  SLOT_PREPARING,
  SLOT_RUNNING,
  type EngineCommand,
  type EngineEvent,
  type EngineInit,
  type LinkRole,
  type WorkerPrepare,
  type WorkerRun,
  type WorkerRunResult,
} from '../../../packages/losat/src/infra/engine-worker/protocol';

class FakeWorker {
  readonly posted: EngineCommand[] = [];
  terminated = false;
  onmessage: ((event: MessageEvent<EngineEvent>) => void) | null = null;
  onerror: ((event: ErrorEvent) => void) | null = null;
  onmessageerror: (() => void) | null = null;

  postMessage(command: EngineCommand): void {
    this.posted.push(command);
  }

  terminate(): void {
    this.terminated = true;
  }

  get init(): EngineInit {
    return this.posted.find((command): command is EngineInit => command.type === 'init')!;
  }

  get runs(): WorkerRun[] {
    return this.posted.filter((command): command is WorkerRun => command.type === 'run');
  }

  get prepares(): WorkerPrepare[] {
    return this.posted.filter((command): command is WorkerPrepare => command.type === 'prepare');
  }

  /** The use count of the threaded instance, as the Engine worker and its thread workers see it. */
  get use(): InstanceUse {
    return new InstanceUse(this.init.use!);
  }

  send(event: EngineEvent): void {
    this.onmessage?.({ data: event } as MessageEvent<EngineEvent>);
  }
}

/**
 * Web Locks as the Gateway uses them to learn that a worker ended (worker-life.ts): a request waits until `release`, as
 * the page's request for the lock of a worker waits until the worker ends, or until its signal aborts.
 */
class FakeLocks implements LifeLocks {
  readonly waiting = new Map<string, { readonly signal: AbortSignal | undefined; readonly grant: () => void }>();

  request(name: string, options: { readonly signal?: AbortSignal }, callback: () => unknown): Promise<unknown> {
    return new Promise((resolve, reject) => {
      this.waiting.set(name, { signal: options.signal, grant: () => resolve(callback()) });
      options.signal?.addEventListener('abort', () => {
        this.waiting.delete(name);
        reject(new DOMException('aborted', 'AbortError'));
      });
    });
  }

  /** The worker that held `name` ended: the page's request is granted. */
  release(name: string): void {
    const waiting = this.waiting.get(name);
    this.waiting.delete(name);
    waiting?.grant();
  }
}

const ASSETS: EngineAssets = {
  serial: { url: '/serial.wasm', sha256: 'a'.repeat(64), size: 1, artifactSha256: 'a'.repeat(64) },
  threads: { url: '/threads.wasm', sha256: 'b'.repeat(64), size: 1, artifactSha256: 'c'.repeat(64), transform: 'shared-memory guard' },
};

const ports: MessagePort[] = [];
afterEach(() => {
  for (const port of ports.splice(0)) port.close();
});

function setup(renewal?: { maxRuns?: number; highWaterBytes?: number }, extra: Partial<WasmEngineOptions> = {}) {
  const workers: FakeWorker[] = [];
  const engine = new WasmEngine({
    ...extra,
    assets: ASSETS,
    control: { describe: () => Promise.reject(new Error('unused')), start: () => Promise.reject(new Error('unused')), validate: () => Promise.reject(new Error('unused')) },
    hardwareConcurrency: 8,
    ...(renewal === undefined ? {} : { renewal }),
    createWorker: () => {
      const worker = new FakeWorker();
      workers.push(worker);
      return worker as unknown as Worker;
    },
    loadModules: async () => ({ serial: { asset: ASSETS.serial }, threads: { asset: ASSETS.threads } }),
  });
  return { engine, workers };
}

function request(runId: string, threads: number | 'auto' = 1): EngineRunRequest {
  const input = { bytes: new Uint8Array([62, 115, 10, 65, 10]), sha256: 'f'.repeat(64), revisionIds: ['r1'], records: [{ id: 's', length: 1 }] };
  return { runId, argv: ['blastn', '-query', 'q.fa', '-subject', 's.fa'], query: input, subject: input, requestedThreads: threads };
}

function output(): MessagePort {
  const channel = new MessageChannel();
  ports.push(channel.port1, channel.port2);
  return channel.port1;
}

/**
 * A link of `worker`'s generation (link.ts), as the Engine worker posts one for a thread worker or the keeper: returns
 * the end that the thread worker or the keeper holds. `answer` gets each message on it and returns what to post back,
 * if anything.
 */
function link(worker: FakeWorker, role: LinkRole, answer: (message: unknown) => unknown = () => undefined): MessagePort {
  const channel = new MessageChannel();
  ports.push(channel.port1, channel.port2);
  channel.port1.onmessage = (event: MessageEvent) => {
    const reply = answer(event.data);
    if (reply !== undefined) channel.port1.postMessage(reply);
  };
  worker.send({ type: 'link', generation: worker.init.generation, role, port: channel.port2 });
  return channel.port1;
}

/** As `link`; the main thread's end records `name` in `closed` when the Gateway closes it. */
function watchedLink(
  worker: FakeWorker,
  role: LinkRole,
  name: string,
  closed: string[],
  answer: (message: unknown) => unknown = () => undefined,
): MessagePort {
  const channel = new MessageChannel();
  ports.push(channel.port1, channel.port2);
  channel.port1.onmessage = (event: MessageEvent) => {
    const reply = answer(event.data);
    if (reply !== undefined) channel.port1.postMessage(reply);
  };
  const close = channel.port2.close.bind(channel.port2);
  channel.port2.close = () => {
    closed.push(name);
    close();
  };
  worker.send({ type: 'link', generation: worker.init.generation, role, port: channel.port2 });
  return channel.port1;
}

/** The keeper of `worker`'s generation on its link: `answer` gets each message and returns what to post back, if anything. */
function fakeKeeper(worker: FakeWorker, answer: (message: unknown) => unknown): MessagePort {
  return link(worker, 'keeper', answer);
}

/** Waits until `condition` holds. */
async function until(condition: () => boolean): Promise<void> {
  for (let i = 0; i < 1000 && !condition(); i++) await new Promise((resolve) => setTimeout(resolve, 1));
  expect(condition()).toBe(true);
}

const RESULT: WorkerRunResult = {
  path: 'serial',
  threads: 1,
  engineBuild: 'losat-web-serial.wasm sha256:aaaaaaaaaaaaaaaa',
  memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 },
  subjectRetained: false,
  messages: 'diagnostics',
};

/** Waits until `worker` has received `count` run commands. */
async function posted(worker: () => FakeWorker | undefined, count: number): Promise<FakeWorker> {
  for (let i = 0; i < 100; i++) {
    const current = worker();
    if (current !== undefined && current.runs.length >= count) return current;
    await new Promise((resolve) => setTimeout(resolve, 1));
  }
  throw new Error('the run was not posted');
}

/** Waits until `worker` has received `count` prepare commands. */
async function prepared(worker: () => FakeWorker | undefined, count: number): Promise<FakeWorker> {
  for (let i = 0; i < 100; i++) {
    const current = worker();
    if (current !== undefined && current.prepares.length >= count) return current;
    await new Promise((resolve) => setTimeout(resolve, 1));
  }
  throw new Error('the preparation was not posted');
}

describe('WasmEngine', () => {
  it('starts a generation, relays phases and resolves with the runtime of the search', async () => {
    const { engine, workers } = setup();
    const phases: string[] = [];
    const run = engine.run(request('r1', 4), output(), (phase) => phases.push(phase));
    const worker = await posted(() => workers[0], 1);
    expect(worker.init.generation).toBe(1);
    expect(worker.init.builds.threads).toBe('losat-web-threads.wasm sha256:cccccccccccccccc (shared-memory guard)');
    const command = worker.runs[0]!;
    expect(command).toMatchObject({ generation: 1, runId: 'r1', program: 'blastn', threads: 4 });
    expect(command.subject!.key).toBe(JSON.stringify(['blastn', 'f'.repeat(64), ['r1']]));
    worker.send({ type: 'phase', generation: 1, runId: 'r1', phase: 'running' });
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...RESULT, subjectRetained: true } });
    await expect(run).resolves.toEqual<RuntimeInfo>({
      path: 'serial',
      threads: 1,
      engineBuild: RESULT.engineBuild,
      runtimeGeneration: 1,
      memory: RESULT.memory,
      subjectRetained: true,
      messages: 'diagnostics',
    });
    expect(phases).toEqual(['running']);
  });

  it('runs the next search in the same generation', async () => {
    const { engine, workers } = setup();
    const first = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await first;
    const second = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[0], 2);
    worker.send({ type: 'done', generation: 1, runId: 'r2', ok: true, result: RESULT });
    expect((await second).runtimeGeneration).toBe(1);
    expect(workers).toHaveLength(1);
  });

  it('while the threaded instance is in use, ends the worker once the keeper of its memory grew it, and at once when not', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => {});
    const worker = await posted(() => workers[0], 1);
    const asked: Array<{ message: unknown; terminated: boolean }> = [];
    const keeper = fakeKeeper(worker, (message) => {
      asked.push({ message, terminated: worker.terminated });
      return { type: 'grown' };
    });
    // The Engine worker calls into the threaded instance; no message tells the Gateway.
    const use = worker.use;
    expect(use.enter()).toBe(true);
    engine.cancel('r1');
    // The search settles at once; the worker ends after the keeper's answer.
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    // A call that would start now, such as the next one of the search, does not.
    expect(use.enter()).toBe(false);
    expect(use.retiring()).toBe(true);
    await until(() => worker.terminated);
    expect(asked).toEqual([{ message: { type: 'grow' }, terminated: false }]);
    expect(ends).toEqual(['grown']);
    keeper.close();

    // After its search ended, a cancel of the next search before it calls into the threaded
    // instance ends the worker at once, and no call can start then.
    const next = engine.run(request('r2', 4), output(), () => {});
    const second = await posted(() => workers[1], 1);
    const other = fakeKeeper(second, () => {
      asked.push({ message: 'second keeper asked', terminated: second.terminated });
      return { type: 'grown' };
    });
    const secondUse = second.use;
    expect(secondUse.enter()).toBe(true);
    secondUse.leave();
    second.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    await next;
    const third = engine.run(request('r3', 4), output(), () => {});
    await posted(() => workers[1], 2);
    engine.cancel('r3');
    await expect(third).rejects.toBeInstanceOf(RunCancelledError);
    expect(second.terminated).toBe(true);
    expect(secondUse.enter()).toBe(false);
    expect(ends).toEqual(['grown', 'at-once']);
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(asked).toHaveLength(1);
    other.close();
  });

  it('a warm-up on the threaded instance is a use: a search queued behind it ends once the keeper grew the memory', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    engine.prepare({ requestedThreads: 4, inputBytes: 10, argv: ['tblastx', '-query', 'q.fa', '-subject', 's.fa'] });
    const worker = await prepared(() => workers[0], 1);
    let asked = 0;
    const keeper = fakeKeeper(worker, () => {
      asked++;
      return { type: 'grown' };
    });
    // The Engine worker warms the threaded instance up; the search waits behind it.
    const use = worker.use;
    expect(use.enter()).toBe(true);
    const run = engine.run(request('r1', 4), output(), () => {});
    await posted(() => workers[0], 1);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    await until(() => worker.terminated);
    expect(asked).toBe(1);
    expect(ends).toEqual(['grown']);
    keeper.close();
  });

  it('a thread of the pool that has not returned after its search keeps the instance in use', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => {});
    const worker = await posted(() => workers[0], 1);
    const keeper = fakeKeeper(worker, () => ({ type: 'grown' }));
    const use = worker.use;
    expect(use.enter()).toBe(true);
    // As ThreadHost.spawn counts a thread: its slot runs until the thread returns.
    const slot = new Int32Array(new SharedArrayBuffer(4));
    Atomics.store(slot, 0, SLOT_RUNNING);
    use.startThread();
    use.leave();
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await run;
    const next = engine.run(request('r2', 4), output(), () => {});
    await posted(() => workers[0], 2);
    engine.cancel('r2');
    await expect(next).rejects.toBeInstanceOf(RunCancelledError);
    await until(() => worker.terminated);
    expect(ends).toEqual(['grown']);
    // The thread returns while the runtime ends: it is counted off once.
    expect(use.endThread(slot, SLOT_PREPARING)).toBe(true);
    expect(use.uses()).toBe(0);
    keeper.close();
  });

  it('waits for the keeper however long its grow takes, also after the search returned', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => {});
    const worker = await posted(() => workers[0], 1);
    let asked = 0;
    const keeper = fakeKeeper(worker, () => {
      asked++;
      return undefined;
    });
    const use = worker.use;
    expect(use.enter()).toBe(true);
    vi.useFakeTimers({ toFake: ['setTimeout', 'clearTimeout'] });
    try {
      engine.cancel('r1');
      await expect(run).rejects.toBeInstanceOf(RunCancelledError);
      await vi.advanceTimersByTimeAsync(60_000);
    } finally {
      vi.useRealTimers();
    }
    await until(() => asked === 1);
    // The ended search returned meanwhile: its threads stopped, but the keeper may still be
    // growing the memory, and it is a thread on the memory too.
    use.leave();
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'Error', message: 'not stopped' } });
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'EngineStoppedError', message: 'trapped', stopped: true } });
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(worker.terminated).toBe(false);
    expect(ends).toEqual([]);
    keeper.postMessage({ type: 'grown' });
    await until(() => worker.terminated);
    expect(ends).toEqual(['grown']);
    keeper.close();
  });

  it('a thread fault during a threaded search ends the worker once the keeper grew the memory', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    let terminatedWhenAsked: boolean | undefined;
    const keeper = fakeKeeper(worker, () => {
      terminatedWhenAsked = worker.terminated;
      return { type: 'grown' };
    });
    expect(worker.use.enter()).toBe(true);
    link(worker, 'thread').postMessage({ type: 'fault', tid: 2, message: 'The engine stopped: thread 2 trapped' });
    await expect(run).rejects.toThrow('The engine stopped: thread 2 trapped');
    await until(() => worker.terminated);
    expect(terminatedWhenAsked).toBe(false);
    keeper.close();
  });

  it('ends the worker of an instance that stopped during a threaded search once the keeper grew the memory', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    let answer: () => void = () => undefined;
    const keeper = fakeKeeper(worker, () => {
      answer = () => keeper.postMessage({ type: 'grown' });
      return undefined;
    });
    // The call stopped, and a thread of its search still runs: the Engine worker kept it and the keeper.
    const use = worker.use;
    expect(use.enter()).toBe(true);
    use.startThread();
    use.leave();
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'EngineStoppedError', message: 'trapped', stopped: true } });
    await expect(run).rejects.toThrow('trapped');
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(worker.terminated).toBe(false);
    answer();
    await until(() => worker.terminated);
    expect(ends).toEqual(['grown']);
    keeper.close();
  });

  it('ends the worker of an instance that stopped at once when its threaded instance is not in use', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'EngineStoppedError', message: 'trapped', stopped: true } });
    await expect(run).rejects.toThrow('trapped');
    expect(worker.terminated).toBe(true);
    expect(ends).toEqual(['at-once']);
  });

  it('a page that cannot share memory gives the worker no use count, and its runtime ends at once', async () => {
    const ends: string[] = [];
    vi.stubGlobal('crossOriginIsolated', false);
    try {
      const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
      const run = engine.run(request('r1', 4), output(), () => {});
      const worker = await posted(() => workers[0], 1);
      expect(worker.init.use).toBeUndefined();
      engine.cancel('r1');
      await expect(run).rejects.toBeInstanceOf(RunCancelledError);
      expect(worker.terminated).toBe(true);
      expect(ends).toEqual(['at-once']);
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it('an observer that throws does not keep a cancel from settling', async () => {
    const { engine, workers } = setup(undefined, {
      onRuntimeEnd: () => {
        throw new Error('observer');
      },
    });
    const run = engine.run(request('r1'), output(), () => {});
    const worker = await posted(() => workers[0], 1);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    expect(worker.terminated).toBe(true);
  });

  it('cancels by ending the worker; the next search starts a new generation and late messages are dropped', async () => {
    const { engine, workers } = setup();
    const phases: string[] = [];
    const run = engine.run(request('r1'), output(), (phase) => phases.push(phase));
    const old = await posted(() => workers[0], 1);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    expect(old.terminated).toBe(true);
    expect(old.onmessage).toBeNull();

    const next = engine.run(request('r2'), output(), (phase) => phases.push(phase));
    const worker = await posted(() => workers[1], 1);
    expect(worker.init.generation).toBe(2);
    // A message of the old generation, as a late one from the ended worker would be.
    worker.send({ type: 'phase', generation: 1, runId: 'r2', phase: 'running' });
    worker.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
    expect(phases).toEqual([]);
  });

  it('ignores a cancel of another or a finished run', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    engine.cancel('other');
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await run;
    engine.cancel('r1');
    expect(worker.terminated).toBe(false);
  });

  it('renews the runtime before the next search when the instance ran too many searches', async () => {
    const { engine, workers } = setup({ maxRuns: 2 });
    for (const [i, runs] of [1, 2].entries()) {
      const run = engine.run(request(`r${i}`), output(), () => undefined);
      const worker = await posted(() => workers[0], i + 1);
      worker.send({ type: 'done', generation: 1, runId: `r${i}`, ok: true, result: { ...RESULT, memory: { ...RESULT.memory, instanceRuns: runs } } });
      await run;
    }
    expect(workers[0]!.terminated).toBe(false);
    const run = engine.run(request('r3'), output(), () => undefined);
    const worker = await posted(() => workers[1], 1);
    expect(workers[0]!.terminated).toBe(true);
    worker.send({ type: 'done', generation: 2, runId: 'r3', ok: true, result: RESULT });
    expect((await run).runtimeGeneration).toBe(2);
  });

  it('renews the runtime when the memory of the instance reached the high-water mark', async () => {
    const { engine, workers } = setup({ highWaterBytes: 100 });
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...RESULT, memory: { ...RESULT.memory, linearBytesAfter: 100 } } });
    await run;
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[1], 1);
    expect(workers[0]!.terminated).toBe(true);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    await next;
  });

  it('allows the high-water mark for the counts file that the instance holds, so that the counts are not read again after every search', async () => {
    const { engine, workers } = setup({ highWaterBytes: 100 });
    // 249 bytes of memory after the search: past the mark of 100, but under 100 + 2 x 75 for the counts.
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...RESULT, countsBytes: 75, memory: { ...RESULT.memory, linearBytesAfter: 249 } } });
    await run;
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[0], 2);
    expect(workers[0]!.terminated).toBe(false);
    // At that mark it is renewed: the counts go with the instance.
    workers[0]!.send({ type: 'done', generation: 1, runId: 'r2', ok: true, result: { ...RESULT, countsBytes: 75, memory: { ...RESULT.memory, linearBytesAfter: 250 } } });
    await next;
    const third = engine.run(request('r3'), output(), () => undefined);
    await posted(() => workers[1], 1);
    expect(workers[0]!.terminated).toBe(true);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r3', ok: true, result: RESULT });
    await third;
  });

  it('keeps the name, role and detail of an InputMismatchError', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({
      type: 'done',
      generation: 1,
      runId: 'r1',
      ok: false,
      error: { name: 'InputMismatchError', message: 'ignored', role: 'subject', detail: 'record 1 is s2 in the table' },
    });
    const error = await run.catch((e: unknown) => e);
    expect(error).toBeInstanceOf(InputMismatchError);
    expect(error).toMatchObject({
      name: 'InputMismatchError',
      role: 'subject',
      message: 'The subject records that the engine read differ from the record table: record 1 is s2 in the table',
    });
    expect(worker.terminated).toBe(false);
  });

  it('forwards assembled on WorkerInput and readRecords on RuntimeInfo (D17)', async () => {
    const { engine, workers } = setup();
    const base = request('r1');
    const run = engine.run(
      { ...base, query: { ...base.query!, assembled: true } },
      output(),
      () => undefined,
    );
    const worker = await posted(() => workers[0], 1);
    const postedRun = worker.posted.find((m) => m.type === 'run')!;
    expect(postedRun.type === 'run' && postedRun.query?.assembled).toBe(true);
    expect(postedRun.type === 'run' && postedRun.subject?.assembled).toBeUndefined();
    const readRecords = { query: [{ id: 'q1', length: 120 }] };
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...RESULT, readRecords } });
    expect((await run).readRecords).toEqual(readRecords);
  });

  it('ends the runtime after an error that stopped the instance', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'EngineStoppedError', message: 'The engine stopped: unreachable', stopped: true } });
    await expect(run).rejects.toMatchObject({ name: 'EngineStoppedError', message: 'The engine stopped: unreachable' });
    expect(worker.terminated).toBe(true);
  });

  // Line WQA's sweep, D2, its review, M1 (2026-11-19): the worker said that it started only once it held its lock, which
  // a search posted with `init` came before, so a worker that failed during its first search read as one that could not
  // start ("LOSAT could not start the part that runs searches. Reload the page."). It says so at once now.
  it('fails the search when the worker fails, also before it holds its lock, and starts a new generation for the next one', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'started', generation: 1 });
    worker.onerror?.({ message: 'out of memory', preventDefault: () => undefined } as ErrorEvent);
    await expect(run).rejects.toThrow('The engine worker stopped: out of memory');
    expect(worker.terminated).toBe(true);
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[1], 1);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
  });

  // Line WQA's sweep, D1 (2026-11-19): a worker script that could not be loaded gives an error event without a message,
  // and the search failed with "The engine worker stopped: an error without a message".
  it('fails the search with LOSAT\'s own failure when the worker could not start, and starts a new one for the next search', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.onerror?.({ message: '', preventDefault: () => undefined } as ErrorEvent);
    await expect(run).rejects.toMatchObject({ name: 'UnavailableError', message: ENGINE_NOT_STARTED });
    expect(worker.terminated).toBe(true);
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[1], 1);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
  });

  // Line WQA's sweep, D2 (2026-11-19): a worker that ends without an error event (Worker.terminate(), self.close(), a
  // browser that ends it) left its search "Searching" for good. The release of the lock that it held tells.
  it('fails the search when its worker ended without an error event, which the release of its lock tells', async () => {
    const locks = new FakeLocks();
    const { engine, workers } = setup(undefined, { locks });
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const name = worker.init.lifeLock;
    expect(name).toMatch(/^losat-web:engine-life:/);
    // The page waits for the lock only once the worker holds it.
    expect(locks.waiting.has(name)).toBe(false);
    worker.send({ type: 'started', generation: 1 });
    expect(locks.waiting.has(name)).toBe(false);
    worker.send({ type: 'life', generation: 1, held: true });
    expect(locks.waiting.has(name)).toBe(true);
    locks.release(name);
    await expect(run).rejects.toThrow(ENGINE_ENDED);
    expect(worker.terminated).toBe(true);
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[1], 1);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
  });

  // The second review of line WQA's D1 and D2, L1 (2026-11-19): the threaded instance was in use, so the Gateway waited
  // for the keeper to grow the memory; Firefox had ended the keeper with the worker, and the wait never ended.
  it('ends at once the generation of a worker that ended during a threaded search, without waiting for its keeper', async () => {
    const locks = new FakeLocks();
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { locks, onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const asked: unknown[] = [];
    // The keeper ended with the worker: it answers nothing.
    const keeper = fakeKeeper(worker, (message) => {
      asked.push(message);
      return undefined;
    });
    worker.send({ type: 'started', generation: 1 });
    worker.send({ type: 'life', generation: 1, held: true });
    expect(worker.use.enter()).toBe(true);
    locks.release(worker.init.lifeLock);
    await expect(run).rejects.toThrow(ENGINE_ENDED);
    expect(worker.terminated).toBe(true);
    expect(ends).toEqual(['at-once']);
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(asked).toEqual([]);
    keeper.close();
  });

  // The third review of line WQA's D1 and D2, L2 (2026-11-19): the Gateway stopped waiting for the worker's lock once the
  // generation began to end, so a generation already waiting for its keeper was never told that its worker ended: the
  // wait never ended, with two channels open. It waits for the lock until the generation has ended.
  describe('ends at once a generation that waits for its keeper when its worker ends', () => {
    /** A threaded search whose worker holds its lock (unless `life` is false) and calls into the threaded instance. */
    async function threadedSearch({ life = true } = {}) {
      const locks = new FakeLocks();
      const ends: string[] = [];
      const { engine, workers } = setup(undefined, { locks, onRuntimeEnd: (end) => ends.push(end) });
      const run = engine.run(request('r1', 4), output(), () => undefined);
      const worker = await posted(() => workers[0], 1);
      const asked: unknown[] = [];
      // The keeper has not answered yet; it ends with the worker.
      const keeper = fakeKeeper(worker, (message) => {
        asked.push(message);
        return undefined;
      });
      worker.send({ type: 'started', generation: 1 });
      if (life) worker.send({ type: 'life', generation: 1, held: true });
      expect(worker.use.enter()).toBe(true);
      return { engine, run, worker, locks, ends, asked, keeper };
    }

    /** The generation waits for the keeper; the worker's end ends it at once. */
    async function endsAtOnceWhenTheWorkerEnds({ worker, locks, ends, asked, keeper }: Awaited<ReturnType<typeof threadedSearch>>) {
      await until(() => asked.length === 1);
      expect(worker.terminated).toBe(false);
      expect(ends).toEqual([]);
      expect(locks.waiting.has(worker.init.lifeLock)).toBe(true);
      locks.release(worker.init.lifeLock);
      expect(worker.terminated).toBe(true);
      expect(ends).toEqual(['at-once']);
      keeper.close();
    }

    it('after a cancel', async () => {
      const search = await threadedSearch();
      search.engine.cancel('r1');
      await expect(search.run).rejects.toBeInstanceOf(RunCancelledError);
      await endsAtOnceWhenTheWorkerEnds(search);
    });

    it('after an error of the worker', async () => {
      const search = await threadedSearch();
      search.worker.onerror?.({ message: 'out of memory', preventDefault: () => undefined } as ErrorEvent);
      await expect(search.run).rejects.toThrow('The engine worker stopped: out of memory');
      await endsAtOnceWhenTheWorkerEnds(search);
    });

    it('after a cancel that came before the worker said that it holds its lock', async () => {
      const search = await threadedSearch({ life: false });
      search.engine.cancel('r1');
      await expect(search.run).rejects.toBeInstanceOf(RunCancelledError);
      search.worker.send({ type: 'life', generation: 1, held: true });
      await endsAtOnceWhenTheWorkerEnds(search);
    });

    it('and stops waiting for the lock once the keeper answered', async () => {
      const search = await threadedSearch();
      const keeper = fakeKeeper(search.worker, () => ({ type: 'grown' }));
      search.engine.cancel('r1');
      await expect(search.run).rejects.toBeInstanceOf(RunCancelledError);
      await until(() => search.worker.terminated);
      expect(search.ends).toEqual(['grown']);
      expect(search.locks.waiting.size).toBe(0);
      keeper.close();
      search.keeper.close();
    });
  });

  it('does not wait for the lock of a worker that it ends itself, nor of one that holds none', async () => {
    const locks = new FakeLocks();
    const { engine, workers } = setup(undefined, { locks });
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'started', generation: 1 });
    worker.send({ type: 'life', generation: 1, held: true });
    expect(locks.waiting.size).toBe(1);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    expect(locks.waiting.size).toBe(0);
    const next = engine.run(request('r2'), output(), () => undefined);
    const second = await posted(() => workers[1], 1);
    // Without Web Locks (or a refused request) the worker holds no lock: none is waited for.
    second.send({ type: 'started', generation: 2 });
    second.send({ type: 'life', generation: 2, held: false });
    expect(locks.waiting.size).toBe(0);
    second.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
  });

  it('fails the search when a thread worker reports a fault on its link', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    link(worker, 'thread').postMessage({ type: 'fault', tid: 2, message: 'The engine stopped: thread 2 trapped' });
    await expect(run).rejects.toThrow('The engine stopped: thread 2 trapped');
    expect(worker.terminated).toBe(true);
  });

  it('fails the search when the Engine worker reports a fault of its thread host', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'fault', generation: 1, tid: 3, message: 'The engine stopped: thread 3 did not start in time' });
    await expect(run).rejects.toThrow('The engine stopped: thread 3 did not start in time');
    expect(worker.terminated).toBe(true);
  });

  it("answers the hello on each link, which the thread worker or the keeper waits for before it is ready", async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    for (const role of ['thread', 'keeper'] as const) {
      const answers: unknown[] = [];
      link(worker, role, (message) => void answers.push(message)).postMessage({ type: 'hello' });
      await until(() => answers.length === 1);
      expect(answers).toEqual([{ type: 'hello' }]);
    }
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await run;
  });

  it('closes the links of a generation when it ends', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const closed: string[] = [];
    watchedLink(worker, 'keeper', 'keeper', closed);
    watchedLink(worker, 'thread', 'thread', closed);
    expect(closed).toEqual([]);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    expect(worker.terminated).toBe(true);
    expect(closed).toEqual(['keeper', 'thread']);
  });

  // The review of line WQA's D8, L1 (2026-11-20): the runtime makes a new ThreadHost after a search on the other module or
  // a failed preparation, and the links of the one before, whose workers had ended, stayed open until the generation
  // ended: with Auto threads alternating small and large searches, up to 32 more for each switch.
  it('closes the links of a thread host once the next host links its keeper, and asks only that keeper to grow', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const closed: string[] = [];
    const asked: string[] = [];
    const keeper = (name: string) =>
      watchedLink(worker, 'keeper', name, closed, () => {
        asked.push(name);
        return { type: 'grown' };
      });
    keeper('keeper 1');
    watchedLink(worker, 'thread', 'thread 1', closed);
    watchedLink(worker, 'thread', 'thread 2', closed);
    expect(closed).toEqual([]);
    keeper('keeper 2');
    expect(closed).toEqual(['keeper 1', 'thread 1', 'thread 2']);
    watchedLink(worker, 'thread', 'thread 3', closed);
    expect(worker.use.enter()).toBe(true);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    await until(() => worker.terminated);
    expect(asked).toEqual(['keeper 2']);
    expect(closed).toEqual(['keeper 1', 'thread 1', 'thread 2', 'keeper 2', 'thread 3']);
  });

  it('prepares the runtime of a search before it comes, once per thread count, and runs the search there', async () => {
    const { engine, workers } = setup();
    engine.prepare({ requestedThreads: 'auto', inputBytes: 2_000_000 });
    engine.prepare({ requestedThreads: 'auto', inputBytes: 2_000_000 });
    const worker = await prepared(() => workers[0], 1);
    expect(worker.init.generation).toBe(1);
    engine.prepare({ requestedThreads: 2, inputBytes: 10 });
    await prepared(() => workers[0], 2);
    // Auto on 8 hardware threads: 4 for an input of 1 MB or more.
    expect(worker.prepares).toEqual([
      { type: 'prepare', generation: 1, threads: 4 },
      { type: 'prepare', generation: 1, threads: 2 },
    ]);
    const run = engine.run(request('r1', 2), output(), () => undefined);
    await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    expect((await run).runtimeGeneration).toBe(1);
    expect(workers).toHaveLength(1);
    // A search may change what its runtime holds, so the same hint prepares again after it.
    engine.prepare({ requestedThreads: 2, inputBytes: 10 });
    await prepared(() => workers[0], 3);
  });

  it('passes the argv of the hint to the worker for the warm-up, and prepares again for another argv', async () => {
    const { engine, workers } = setup();
    const argv = ['tblastx', '-query', 'q.fa', '-subject', 's.fa'];
    engine.prepare({ requestedThreads: 2, inputBytes: 10, argv });
    engine.prepare({ requestedThreads: 2, inputBytes: 10, argv: [...argv] });
    const worker = await prepared(() => workers[0], 1);
    engine.prepare({ requestedThreads: 2, inputBytes: 10, argv: [...argv, '-seg', 'no'] });
    await prepared(() => workers[0], 2);
    expect(worker.prepares).toEqual([
      { type: 'prepare', generation: 1, threads: 2, argv },
      { type: 'prepare', generation: 1, threads: 2, argv: [...argv, '-seg', 'no'] },
    ]);
  });

  it('does not prepare while a search runs', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    engine.prepare({ requestedThreads: 2, inputBytes: 10 });
    await new Promise((resolve) => setTimeout(resolve, 5));
    expect(worker.prepares).toEqual([]);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await run;
  });

  it('prepares the renewal that the next search would make', async () => {
    const { engine, workers } = setup({ maxRuns: 1 });
    const run = engine.run(request('r1'), output(), () => undefined);
    const old = await posted(() => workers[0], 1);
    old.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await run;
    engine.prepare({ requestedThreads: 2, inputBytes: 10 });
    const worker = await prepared(() => workers[1], 1);
    expect(old.terminated).toBe(true);
    expect(worker.prepares).toEqual([{ type: 'prepare', generation: 2, threads: 2 }]);
  });

  it('runs one search at a time', async () => {
    const { engine, workers } = setup();
    const first = engine.run(request('r1'), output(), () => undefined);
    await expect(engine.run(request('r2'), output(), () => undefined)).rejects.toThrow('one search at a time');
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    await first;
  });
});

// Line WQA's sweep, D8 (2026-11-20): where the site's data is blocked (cookies blocked), Firefox refuses
// `new BroadcastChannel()` with SecurityError, on the page and in every worker. The Gateway opened one for the faults of
// the thread workers and one for the keeper with each runtime, and every search failed with "The operation is
// insecure." The thread workers and the keeper are linked to the Gateway by MessageChannels (link.ts).
describe("WasmEngine where the site's data is blocked", () => {
  beforeEach(() => {
    vi.stubGlobal(
      'BroadcastChannel',
      class {
        constructor() {
          throw new DOMException('The operation is insecure.', 'SecurityError');
        }
      },
    );
  });
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  it('runs a search', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    // The browser refuses the worker's lock too: none is waited for.
    worker.send({ type: 'started', generation: 1 });
    worker.send({ type: 'life', generation: 1, held: false });
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: RESULT });
    expect((await run).runtimeGeneration).toBe(1);
  });

  it('fails a threaded search whose thread stopped, which its thread worker reports on its link', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const keeper = fakeKeeper(worker, () => ({ type: 'grown' }));
    expect(worker.use.enter()).toBe(true);
    link(worker, 'thread').postMessage({ type: 'fault', tid: 2, message: 'The engine stopped: thread 2 trapped' });
    await expect(run).rejects.toThrow('The engine stopped: thread 2 trapped');
    await until(() => worker.terminated);
    keeper.close();
  });

  it('ends a cancelled threaded search in use once the keeper grew the memory, asked on its link', async () => {
    const ends: string[] = [];
    const { engine, workers } = setup(undefined, { onRuntimeEnd: (end) => ends.push(end) });
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const asked: Array<{ message: unknown; terminated: boolean }> = [];
    const keeper = fakeKeeper(worker, (message) => {
      asked.push({ message, terminated: worker.terminated });
      return { type: 'grown' };
    });
    expect(worker.use.enter()).toBe(true);
    engine.cancel('r1');
    await expect(run).rejects.toBeInstanceOf(RunCancelledError);
    await until(() => worker.terminated);
    expect(asked).toEqual([{ message: { type: 'grow' }, terminated: false }]);
    expect(ends).toEqual(['grown']);
    keeper.close();
  });
});
