// The Wasm EngineGateway's runtime generations (src/infra/engine-worker/gateway.ts) with a
// fake Engine worker: runs, phases, cancel, renewal, faults and errors. The real Engine
// worker runs in the browsers (tests/e2e/engine.spec.ts).
import { afterEach, describe, expect, it } from 'vitest';
import { InputMismatchError, RunCancelledError, type EngineRunRequest, type RuntimeInfo } from '../../src/ports/engine';
import { WasmEngine } from '../../src/infra/engine-worker/gateway';
import type { EngineCommand, EngineEvent, EngineInit, WorkerRun, WorkerRunResult } from '../../src/infra/engine-worker/protocol';
import type { EngineAssets } from '../../src/infra/reactor/assets';

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

  send(event: EngineEvent): void {
    this.onmessage?.({ data: event } as MessageEvent<EngineEvent>);
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

function setup(renewal?: { maxRuns?: number; highWaterBytes?: number }) {
  const workers: FakeWorker[] = [];
  const engine = new WasmEngine({
    assets: ASSETS,
    control: { describe: () => Promise.reject(new Error('unused')), validate: () => Promise.reject(new Error('unused')) },
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

const RESULT: WorkerRunResult = {
  path: 'serial',
  threads: 1,
  engineBuild: 'losat-web-serial.wasm sha256:aaaaaaaaaaaaaaaa',
  memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 },
  subjectRetained: false,
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
    expect(command.subject.key).toBe(JSON.stringify(['blastn', 'f'.repeat(64), ['r1']]));
    worker.send({ type: 'phase', generation: 1, runId: 'r1', phase: 'running' });
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...RESULT, subjectRetained: true } });
    await expect(run).resolves.toEqual<RuntimeInfo>({
      path: 'serial',
      threads: 1,
      engineBuild: RESULT.engineBuild,
      runtimeGeneration: 1,
      memory: RESULT.memory,
      subjectRetained: true,
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

  it('ends the runtime after an error that stopped the instance', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.send({ type: 'done', generation: 1, runId: 'r1', ok: false, error: { name: 'EngineStoppedError', message: 'The engine stopped: unreachable', stopped: true } });
    await expect(run).rejects.toMatchObject({ name: 'EngineStoppedError', message: 'The engine stopped: unreachable' });
    expect(worker.terminated).toBe(true);
  });

  it('fails the search when the worker fails, and starts a new generation for the next one', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1'), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    worker.onerror?.({ message: 'out of memory', preventDefault: () => undefined } as ErrorEvent);
    await expect(run).rejects.toThrow('The engine worker stopped: out of memory');
    expect(worker.terminated).toBe(true);
    const next = engine.run(request('r2'), output(), () => undefined);
    await posted(() => workers[1], 1);
    workers[1]!.send({ type: 'done', generation: 2, runId: 'r2', ok: true, result: RESULT });
    expect((await next).runtimeGeneration).toBe(2);
  });

  it('fails the search when a thread worker reports a fault', async () => {
    const { engine, workers } = setup();
    const run = engine.run(request('r1', 4), output(), () => undefined);
    const worker = await posted(() => workers[0], 1);
    const channel = new BroadcastChannel(worker.init.faultChannel);
    channel.postMessage({ type: 'fault', tid: 2, message: 'The engine stopped: thread 2 trapped' });
    await expect(run).rejects.toThrow('The engine stopped: thread 2 trapped');
    channel.close();
    expect(worker.terminated).toBe(true);
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
