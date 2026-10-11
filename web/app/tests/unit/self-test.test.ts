// The browser self-test (application/self-test.ts) against the FakeEngine and the Data service
// in memory: it compares each output's SHA-256 with the table, stays out of the queue (the
// Coordinator holds it while the self-test runs, and refuses the self-test while runs are
// active or queued), and leaves no source or run behind. The engine's real outputs are
// compared in self-test-engine.test.ts (with LOSAT_WEB_REACTORS) and in the browsers
// (tests/e2e/onboarding.spec.ts).
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { Coordinator, type SearchRequest } from '../../src/application/coordinator';
import { FAKE_ENGINE_REFUSAL, SelfTest, type SelfTestState } from '../../src/application/self-test';
import { exampleFileName, parseSelfTestTable, type SelfTestTable } from '../../src/domain/examples';
import type { ProgramId } from '../../src/domain/programs';
import type { RunStatus } from '../../src/domain/run';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { EngineGateway, EngineRunRequest } from '../../src/ports/engine';
import type { ExampleSource } from '../../src/ports/examples';
import TABLE from '../../src/examples/self-test.json';
import { memoryDownloader } from './support/memory-downloader';

const EXAMPLES = resolve(import.meta.dirname, '../../src/examples');
const RUNNABLE: readonly ProgramId[] = ['blastn', 'blastp', 'tblastn', 'tblastx'];

const examples: ExampleSource = {
  async files(program) {
    const file = (role: 'query' | 'subject') => new File([readFileSync(resolve(EXAMPLES, exampleFileName(program, role)))], exampleFileName(program, role));
    return { query: file('query'), subject: file('subject') };
  },
};

/** An engine whose runs wait, when `hold` is set, until the test releases them. */
class HeldEngine implements EngineGateway {
  readonly fake = new FakeEngine();
  readonly requests: EngineRunRequest[] = [];
  private releases: Array<() => void> = [];
  hold = false;

  describe = (program: ProgramId) => this.fake.describe(program);
  validate = (argv: readonly string[]) => this.fake.validate(argv);
  cancel = (runId: string) => this.fake.cancel(runId);

  async run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: 'preparing' | 'running' | 'finalizing') => void) {
    this.requests.push(request);
    if (this.hold) await new Promise<void>((release) => this.releases.push(release));
    return this.fake.run(request, output, onPhase);
  }

  releaseAll(): void {
    for (const release of this.releases.splice(0)) release();
  }

  get waiting(): number {
    return this.releases.length;
  }
}

function setup(options: { table?: SelfTestTable; source?: ExampleSource; fake?: boolean } = {}) {
  let token = 0;
  const data = new DataService({
    store: new MemoryBlockStore(),
    scanner: new FakeScanner(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const calls = { released: [] as string[], deleted: [] as string[], discarded: [] as string[] };
  const releaseSources = data.releaseSources.bind(data);
  data.releaseSources = async (ids) => {
    calls.released.push(...ids);
    return releaseSources(ids);
  };
  const deleteRun = data.deleteRun.bind(data);
  data.deleteRun = async (id) => {
    calls.deleted.push(id);
    return deleteRun(id);
  };
  const engine = new HeldEngine();
  let id = 0;
  const coordinator = new Coordinator({
    engine,
    data,
    downloader: memoryDownloader(() => undefined),
    now: () => 0,
    newRunId: () => `run-${++id}`,
  });
  let selfTestId = 0;
  const selfTest = new SelfTest({
    engine,
    data,
    examples: options.source ?? examples,
    table: options.table ?? parseSelfTestTable(TABLE),
    runs: coordinator.state,
    withEngine: (task) => coordinator.withEngine(task),
    digest: sha256Hex,
    browser: () => ({ crossOriginIsolated: true, logicalProcessors: 8 }),
    threads: [4, 1],
    newRunId: () => `self-test-${++selfTestId}`,
    now: () => 0,
    fake: options.fake ?? false,
  });
  return { data, engine, coordinator, selfTest, calls };
}

/** The table with the outputs that a first self-test got (the FakeEngine's), so that a second one passes. */
function tableOf(state: SelfTestState): SelfTestTable {
  const table = parseSelfTestTable(TABLE);
  return {
    ...table,
    cases: table.cases.map((c) => {
      const result = state.cases.find((r) => r.program === c.program);
      if (result?.outcome?.kind !== 'fail') return c;
      const outputs = { ...c.outputs };
      for (const mismatch of result.outcome.mismatches) outputs[mismatch.format] = { ...outputs[mismatch.format], sha256: mismatch.got };
      return { ...c, outputs };
    }),
  };
}

const request: SearchRequest = { program: 'blastn', query: { text: '>q\nACGT\n' }, subject: { text: '>s\nACGT\n' }, parameters: [], requestedThreads: 'auto' };
const statusOf = (coordinator: Coordinator, runId: string): RunStatus | undefined =>
  coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)?.status;

describe('the browser self-test', () => {
  it('runs each runnable example at each thread count and compares the SHA-256 of outfmt 0, 6 and 7', async () => {
    // The FakeEngine's outputs are not the engine's: every case fails, and says what it got.
    const first = setup();
    await first.selfTest.run();
    const failed = first.selfTest.state.get();
    expect(failed.status).toBe('finished');
    expect(failed.cases.map((c) => `${c.program}@${c.threads}:${c.outcome?.kind}`)).toEqual(
      RUNNABLE.flatMap((program) => [`${program}@4:fail`, `${program}@1:fail`]),
    );
    const mismatch = failed.cases[0]!.outcome!;
    expect(mismatch.kind === 'fail' && mismatch.mismatches.map((m) => m.format)).toEqual([0, 6, 7]);
    expect(mismatch.kind === 'fail' && mismatch.mismatches[0]!.expected).toBe(TABLE.cases[0]!.outputs['0'].sha256);
    expect(failed.summary).toBe('0 of 8 cases passed; 8 wrote other outputs than expected.');
    expect(failed.skipped).toEqual([{ program: 'blastx', reason: expect.stringContaining('BLASTX joins LOSAT Web') }]);
    expect(failed.facts).toEqual({ crossOriginIsolated: true, logicalProcessors: 8, storage: { backend: 'memory' } });
    // The engine was asked as the search form asks, at each thread count.
    expect(first.engine.requests.map((r) => [r.argv.join(' '), r.requestedThreads])).toEqual(
      RUNNABLE.flatMap((p) => [4, 1].map((threads) => [`${p} -query example-${p}-query.fasta -subject example-${p}-subject.fasta`, threads])),
    );

    // With the outputs it got as the expected ones, every case passes.
    const second = setup({ table: tableOf(failed) });
    await second.selfTest.run();
    const passed = second.selfTest.state.get();
    expect(passed.cases.every((c) => c.outcome?.kind === 'pass' && c.state === 'done')).toBe(true);
    // What the engine says of the run (the FakeEngine runs everything on one thread).
    expect(passed.cases[0]!.runtime).toEqual({ path: 'fake', threads: 1, engineBuild: 'fake-engine' });
    expect(passed.summary).toBe("All 8 cases passed: this browser's engine wrote the expected outputs.");
  });

  it('leaves nothing of it in the application or the Data worker', async () => {
    const { selfTest, coordinator, data, calls } = setup();
    const before = (await data.storageInfo()).sessionBytes;
    await selfTest.run();
    expect(coordinator.state.get().runs).toEqual([]);
    // Two sources a program, released; every run deleted.
    expect(calls.released).toHaveLength(2 * RUNNABLE.length);
    expect(calls.deleted).toEqual(Array.from({ length: 8 }, (_, i) => `self-test-${i + 1}`));
    expect((await data.storageInfo()).sessionBytes).toBe(before);
    await expect(data.readOutput('self-test-1', 6)).rejects.toThrow();
    await expect(data.buildRunInput(['revision-1'])).rejects.toThrow();
  });

  it('refuses to start while a run is active or queued, and says why', async () => {
    const { selfTest, coordinator, engine } = setup();
    engine.hold = true;
    const first = await coordinator.enqueue(request);
    await coordinator.enqueue(request);
    await expect.poll(() => engine.waiting).toBe(1);
    await selfTest.run();
    expect(selfTest.state.get()).toMatchObject({
      status: 'idle',
      cases: [],
      refusal: 'The self-test starts only when no search is running or waiting: Run 1 is running and 1 run is waiting. Try again when the queue is empty.',
    });
    expect(engine.requests.map((r) => r.runId)).toEqual([first.runId]);
    engine.hold = false;
    engine.releaseAll();
    await expect.poll(() => selfTest.refusal()).toBeUndefined();
  });

  it('holds the queue while it runs: a run added meanwhile waits and then runs', async () => {
    const { selfTest, coordinator, engine } = setup();
    engine.hold = true;
    const running = selfTest.run();
    await expect.poll(() => engine.waiting).toBe(1);
    expect(selfTest.state.get().status).toBe('running');
    expect(selfTest.state.get().cases[0]).toMatchObject({ program: 'blastn', threads: 4, state: 'running' });
    const added = await coordinator.enqueue(request);
    expect(added.ok).toBe(true);
    // A second request of the self-test is ignored while it runs.
    await selfTest.run();
    for (let i = 0; i < 8; i++) {
      await expect.poll(() => engine.waiting).toBe(1);
      expect(statusOf(coordinator, added.runId!)).toBe('queued');
      engine.releaseAll();
    }
    await running;
    expect(selfTest.state.get().status).toBe('finished');
    await expect.poll(() => engine.waiting).toBe(1);
    expect(engine.requests.at(-1)!.runId).toBe(added.runId);
    engine.hold = false;
    engine.releaseAll();
    await expect.poll(() => statusOf(coordinator, added.runId!)).toBe('completed');
  });

  it('reports a case that cannot run as an error, and goes on with the next', async () => {
    const source: ExampleSource = {
      async files(program) {
        if (program === 'blastp') throw new Error('The example file example-blastp-query.fasta could not be read from the site (HTTP 404).');
        const files = await examples.files(program);
        // A TBLASTN subject that is not the release's.
        return program === 'tblastn' ? { ...files, subject: new File(['>other\nACGT\n'], files.subject.name) } : files;
      },
    };
    const { selfTest, calls } = setup({ source });
    await selfTest.run();
    const outcomes = Object.fromEntries(selfTest.state.get().cases.map((c) => [`${c.program}@${c.threads}`, c.outcome]));
    expect(outcomes['blastp@4']).toEqual({ kind: 'error', message: 'The example file example-blastp-query.fasta could not be read from the site (HTTP 404).' });
    expect(outcomes['tblastn@1']).toEqual({
      kind: 'error',
      message: 'The example file example-tblastn-subject.fasta is not the one that the expected outputs were made from.',
    });
    expect(outcomes['tblastx@1']?.kind).toBe('fail');
    expect(selfTest.state.get().summary).toBe('0 of 8 cases passed; 4 wrote other outputs than expected; 4 could not run.');
    // The sources of the refused TBLASTN input were released too.
    expect(calls.released).toHaveLength(2 + 2 + 2);
  });

  it('a run that fails is an error of its case, and leaves no staged run', async () => {
    const { selfTest, engine, data } = setup();
    engine.run = async () => {
      throw new Error('The engine stopped: unreachable');
    };
    await selfTest.run();
    expect(selfTest.state.get().cases.every((c) => c.outcome?.kind === 'error')).toBe(true);
    expect(selfTest.state.get().cases[0]!.outcome).toEqual({ kind: 'error', message: 'The engine stopped: unreachable' });
    expect((await data.storageInfo()).sessionBytes).toBe(0);
  });

  it('says that a FakeEngine build cannot run it', async () => {
    const { selfTest, engine } = setup({ fake: true });
    expect(selfTest.unavailable()).toBe(FAKE_ENGINE_REFUSAL);
    await selfTest.run();
    expect(selfTest.state.get().refusal).toBe(FAKE_ENGINE_REFUSAL);
    expect(engine.requests).toEqual([]);
  });
});

describe('Coordinator.withEngine', () => {
  it('refuses while a run is active or queued, or another task holds the engine', async () => {
    const { coordinator, engine } = setup();
    engine.hold = true;
    await coordinator.enqueue(request);
    await expect.poll(() => engine.waiting).toBe(1);
    await expect(coordinator.withEngine(async () => 'ran')).rejects.toThrow('Runs are active or queued.');
    engine.hold = false;
    engine.releaseAll();
    await expect.poll(() => coordinator.state.get().runs[0]!.status).toBe('completed');
    let release!: () => void;
    const held = coordinator.withEngine(() => new Promise<string>((resolve) => (release = () => resolve('done'))));
    await expect(coordinator.withEngine(async () => 'ran')).rejects.toThrow('Another task is using the engine.');
    release();
    await expect(held).resolves.toBe('done');
    await expect(coordinator.withEngine(async () => 'ran')).resolves.toBe('ran');
  });

  it('starts the queue again after a task that failed', async () => {
    const { coordinator } = setup();
    let fail!: (error: Error) => void;
    const task = coordinator.withEngine(() => new Promise<never>((_, reject) => (fail = reject)));
    const queued = await coordinator.enqueue(request);
    expect(statusOf(coordinator, queued.runId!)).toBe('queued');
    fail(new Error('broken'));
    await expect(task).rejects.toThrow('broken');
    await expect.poll(() => statusOf(coordinator, queued.runId!)).toBe('completed');
  });
});
