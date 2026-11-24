import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import { RunOutputWriter, sha256Hex } from 'losat';
import { CUSTOM_OUTPUT, OUTPUT_FORMATS, type MessagePlace, type OutputFormat, type RegistryFile, type RunFile, type ShownIds } from 'losat/domain';
import {
  DIAGNOSTICS_STREAM,
  runOutputFailedOf,
  UnavailableError,
  withRunFiles,
  withRunOutputFailed,
  type EngineGateway,
  type EnginePhase,
  type EngineRunRequest,
  type RuntimeInfo,
  type StartInfo,
  type StartResult,
  type ValidationResult,
} from 'losat/ports';
import { Coordinator, isToolRun, type RunView, type SearchRequest } from '../../src/application/coordinator';
import type { RunStatus } from '../../src/domain/run';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine, FAKE_MARKER } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import type { DataGateway } from '../../src/ports/data';
import type { Downloader } from '../../src/ports/download';
import { NodeInflater } from './support/node-inflater';

const request: SearchRequest = {
  program: 'blastn',
  query: { text: '>q\nACGT\n' },
  subject: { text: '>s\nACGT\n' },
  parameters: [],
  requestedThreads: 'auto',
};

function setup(
  engine: EngineGateway = new FakeEngine(),
  wrap: (data: DataService) => DataGateway = (d) => d,
  dataWorkerEnded?: Promise<Error>,
) {
  const store = new MemoryBlockStore();
  let token = 0;
  const data = new DataService({
    store,
    scanner: new FakeScanner(),
    inflater: new NodeInflater(),
    checker: new FakeInputChecker(),
    digest: sha256Hex,
    newToken: () => `token-${++token}`,
    cleanup: Promise.resolve({ state: 'done', removedSessions: 0 }),
  });
  const saved: Array<{ fileName: string; text: string; type: string }> = [];
  const downloader: Downloader = {
    save: (fileName, bytes, type) => saved.push({ fileName, text: new TextDecoder().decode(bytes), type }),
  };
  let id = 0;
  const coordinator = new Coordinator({
    engine,
    data: wrap(data),
    downloader,
    now: () => 1000,
    newRunId: () => `run-${++id}`,
    digest: sha256Hex,
    ...(dataWorkerEnded === undefined ? {} : { dataWorkerEnded }),
  });
  return { coordinator, data, store, saved };
}

const viewOf = (coordinator: Coordinator, runId: string): RunView | undefined => {
  const run = coordinator.state.get().runs.find((entry) => entry.snapshot.runId === runId);
  return run !== undefined && !isToolRun(run) ? run : undefined;
};
const statusOf = (coordinator: Coordinator, runId: string): RunStatus | undefined => viewOf(coordinator, runId)?.status;

function waitFor(coordinator: Coordinator, runId: string, status: RunStatus): Promise<void> {
  return new Promise((resolve) => {
    const check = () => {
      if (statusOf(coordinator, runId) === status) {
        unsubscribe();
        resolve();
      }
    };
    const unsubscribe = coordinator.state.subscribe(check);
    check();
  });
}

/** An engine whose runs finish only when the test releases them. */
class ManualEngine implements EngineGateway {
  readonly started: string[] = [];
  readonly requests: EngineRunRequest[] = [];
  /** The registry files that `validate` was given, in order. */
  readonly registries: Array<RegistryFile | undefined> = [];
  /** What the engine's start reads of every argv, or its message. */
  startResult: StartResult = { ok: true, start: { conffile: null, logFile: null, dryRun: false, registerInputs: false, prints: null, database: false } };
  private readonly pending = new Map<
    string,
    { resolve: (runFiles?: readonly RunFile[]) => void; reject: (e: Error) => void; output: MessagePort }
  >();
  lateFinish?: () => void;
  /** The IDs that the outputs of the next runs show (a run of -parse_deflines). */
  shownIds?: ShownIds;
  /** Where the runs' messages go, as the engine says (`RuntimeInfo.messages`). */
  messages: MessagePlace = 'diagnostics';
  /** The text that the next runs print in place of a search (`RuntimeInfo.printed`). */
  printed?: Uint8Array;

  /** How many times the engine was asked to describe a program. */
  described = 0;
  async describe() {
    this.described += 1;
    return { program: 'blastn' as const, formats: [6] as const, parameters: [] };
  }
  async start(): Promise<StartResult> {
    return this.startResult;
  }
  async validate(_argv?: readonly string[], registry?: RegistryFile): Promise<ValidationResult> {
    this.registries.push(registry);
    return { ok: true, dataLoaders: true };
  }
  run(req: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    this.started.push(req.runId);
    this.requests.push(req);
    onPhase('running');
    return new Promise((resolve, reject) => {
      this.pending.set(req.runId, {
        resolve: (runFiles) => {
          const writer = new RunOutputWriter(output);
          writer.write(6, new TextEncoder().encode(`${req.runId}\n`));
          writer.write(DIAGNOSTICS_STREAM, new TextEncoder().encode(`Warning: ${req.runId}\n`));
          writer.end();
          resolve({
            path: 'fake',
            threads: 1,
            engineBuild: 'manual',
            messages: this.messages,
            ...(runFiles === undefined ? {} : { runFiles }),
            ...(this.shownIds === undefined ? {} : { shownIds: this.shownIds }),
            ...(this.printed === undefined ? {} : { printed: this.printed }),
          });
        },
        reject,
        output,
      });
      // Simulates a phase event that arrives after the run was cancelled.
      this.lateFinish = () => onPhase('finalizing');
    });
  }
  /** The runs that were cancelled, in order. */
  readonly cancelled: string[] = [];
  cancel(runId: string): void {
    this.cancelled.push(runId);
  }
  /** Ends a run, with the files that it made besides its outputs. */
  finish(runId: string, runFiles?: readonly RunFile[]) {
    this.pending.get(runId)?.resolve(runFiles);
  }
  /**
   * Fails a run, with the files that it made besides its outputs. With `diagnostics` the run fails as the
   * engine fails a run that began to write: it writes them to stream 3 and ends its output with `fail`,
   * then rejects (with `late`, it rejects first, as a slow port might deliver them after the rejection).
   */
  fail(runId: string, message: string, runFiles?: readonly RunFile[], diagnostics?: string, late = false, kept?: { readonly output: OutputFormat; readonly text: string }) {
    const pending = this.pending.get(runId);
    if (pending === undefined) return;
    const error = new Error(message);
    const failed = diagnostics === undefined ? error : withRunOutputFailed(error);
    const write = () => {
      const writer = new RunOutputWriter(pending.output);
      // The output that the CLI leaves before its error, and the others dropped (line WQA's sweep, D15, 2026-11-20).
      if (kept !== undefined) {
        writer.write(kept.output, new TextEncoder().encode(kept.text));
        for (const format of OUTPUT_FORMATS) if (format !== kept.output) writer.discard(format, 'run');
      }
      if (diagnostics !== '') writer.write(DIAGNOSTICS_STREAM, new TextEncoder().encode(diagnostics));
      writer.fail();
    };
    if (diagnostics !== undefined && !late) write();
    pending.reject(runFiles === undefined ? failed : withRunFiles(failed, runFiles));
    if (diagnostics !== undefined && late) setTimeout(write, 50);
  }
}

const startOf = (start: Partial<StartInfo>): StartResult => ({ ok: true, start: { conffile: null, logFile: null, dryRun: false, registerInputs: false, prints: null, database: false, ...start } });
const logFile = (text: string): RunFile => ({ name: 'run.log', bytes: new TextEncoder().encode(text) });

describe('Coordinator', () => {
  it('runs one job at a time in FIFO order and stores the outputs and diagnostics', async () => {
    const engine = new ManualEngine();
    const { coordinator, data } = setup(engine);
    const first = (await coordinator.enqueue(request)).runId!;
    const second = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'running');
    expect(engine.started).toEqual([first]);
    expect(statusOf(coordinator, second)).toBe('queued');

    engine.finish(first);
    await waitFor(coordinator, first, 'completed');
    await waitFor(coordinator, second, 'running');
    expect(engine.started).toEqual([first, second]);
    expect(new TextDecoder().decode(await data.readOutput(first, 6))).toBe(`${first}\n`);
    expect(await data.readDiagnostics(first)).toBe(`Warning: ${first}\n`);
  });

  // Line WQA's sweep, D2 (2026-11-19): a Data worker that did not answer the commit kept the run "Organizing results"
  // for good, with no Cancel, and the queue behind it waited.
  it('a cancel ends a run whose results are not committed, the next run starts, and results that come later are deleted', async () => {
    const engine = new ManualEngine();
    let commit: (value: unknown) => void = () => undefined;
    let made: Promise<unknown> = Promise.resolve();
    const deleted: string[] = [];
    const { coordinator } = setup(engine, (data) =>
      Object.assign(Object.create(data) as DataService, {
        // The first run's commit is made, and its answer comes only when the test lets it.
        commitRun: (runId: string) => {
          const commitment = data.commitRun(runId);
          if (viewOf(coordinator, runId)!.snapshot.number !== 1) return commitment;
          made = commitment;
          return new Promise((resolve) => (commit = resolve)).then(() => commitment);
        },
        deleteRun: async (runId: string) => {
          deleted.push(runId);
          await data.deleteRun(runId);
        },
      }),
    );
    const first = (await coordinator.enqueue(request)).runId!;
    const second = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'running');
    engine.finish(first);
    await waitFor(coordinator, first, 'finalizing');
    // The Data worker committed the results; its answer does not come.
    await made;
    coordinator.cancel(first);
    await waitFor(coordinator, first, 'cancelled');
    await waitFor(coordinator, second, 'running');
    expect(deleted).toEqual([]);
    commit(undefined);
    await new Promise((resolve) => setTimeout(resolve, 10));
    expect(deleted).toEqual([first]);
    expect(statusOf(coordinator, first)).toBe('cancelled');
  });

  // Line WQA's sweep, D1 and D2: a part of LOSAT that could not be loaded or stopped is said as the page says it.
  it('a run whose Data worker stopped fails with what stopped and the remedy', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine, (data) =>
      Object.assign(Object.create(data) as DataService, {
        commitRun: () => Promise.reject(new UnavailableError('The part of LOSAT that keeps the results stopped')),
      }),
    );
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'running');
    engine.finish(runId);
    await waitFor(coordinator, runId, 'failed');
    expect(viewOf(coordinator, runId)!.record.error).toBe('The part of LOSAT that keeps the results stopped. Reload the page.');
  });

  // Line WQA's sweep, D2, its review, L2 (2026-11-19): when the Data worker ended during a search, the page said so at
  // once, but the run went on "Searching" and failed only when the search ended, at the commit of its results.
  it('a run whose Data worker ended during its search fails at once, and its search ends', async () => {
    const engine = new ManualEngine();
    let end: (error: Error) => void = () => undefined;
    const ended = new Promise<Error>((resolve) => (end = resolve));
    const { coordinator } = setup(engine, undefined, ended);
    const first = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'running');
    end(new UnavailableError('The part of LOSAT that keeps the results stopped'));
    await waitFor(coordinator, first, 'failed');
    expect(viewOf(coordinator, first)!.record.error).toBe('The part of LOSAT that keeps the results stopped. Reload the page.');
    expect(engine.cancelled).toEqual([first]);
  });

  it('a run whose Data worker ended while its outputs were staged starts no search', async () => {
    const engine = new ManualEngine();
    let end: (error: Error) => void = () => undefined;
    const ended = new Promise<Error>((resolve) => (end = resolve));
    let stage: () => void = () => undefined;
    const staged = new Promise<void>((resolve) => (stage = resolve));
    const { coordinator } = setup(
      engine,
      (inner) => {
        const gateway = Object.create(inner) as DataGateway;
        // The outputs are staged, and the answer comes after the Data worker ended.
        gateway.openRun = async (runId, outputs) => {
          const port = await inner.openRun(runId, outputs);
          await staged;
          return port;
        };
        return gateway;
      },
      ended,
    );
    const first = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'preparing');
    end(new UnavailableError('The part of LOSAT that keeps the results stopped'));
    await ended;
    stage();
    await waitFor(coordinator, first, 'failed');
    expect(viewOf(coordinator, first)!.record.error).toBe('The part of LOSAT that keeps the results stopped. Reload the page.');
    expect(engine.started).toEqual([]);
  });

  it('refuses a request with the remedy, not as a fault of the options or the input, when a part of LOSAT is unavailable', async () => {
    const engine = new ManualEngine();
    const unavailable = new UnavailableError('LOSAT could not download its engine file /engine/serial.wasm (network error)');
    engine.start = () => Promise.reject(unavailable);
    const { coordinator } = setup(engine);
    expect(await coordinator.enqueue(request)).toEqual({
      ok: false,
      message: 'LOSAT could not download its engine file /engine/serial.wasm (network error). Reload the page.',
    });
    const reading = new ManualEngine();
    const { coordinator: other } = setup(reading, (data) =>
      Object.assign(Object.create(data) as DataService, { indexSource: () => Promise.reject(unavailable) }),
    );
    expect(await other.enqueue(request)).toEqual({
      ok: false,
      message: 'LOSAT could not download its engine file /engine/serial.wasm (network error). Reload the page.',
    });
  });

  it('freezes the snapshot so later edits cannot change a queued job', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    const params: Array<[string, string]> = [['-evalue', '10']];
    const runId = (await coordinator.enqueue({ ...request, parameters: params })).runId!;
    params[0]![1] = '1e-50';
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(Object.isFrozen(snapshot)).toBe(true);
    expect(snapshot.argv).toContain('10');
    expect(snapshot.argv).not.toContain('1e-50');
    // Pasted text is a file that the run names: its own prefix (number, program, time) and the role (DECISIONS 44).
    expect(snapshot.query.name).toBe('losat-run1-blastn-19700101T000001Z.query.fa');
    expect(snapshot.argv.slice(1, 5)).toEqual(['-query', snapshot.query.name, '-subject', snapshot.subject.name]);
    expect(snapshot.subject.sha256).toBe(createHash('sha256').update('>s\nACGT\n').digest('hex'));
    expect(snapshot.subject.records).toEqual([{ id: 's', length: 4 }]);
    expect(snapshot.subject.revisionIds).toHaveLength(1);
  });

  // Line WQA's D22 (2026-11-23): an input that is not a file of the user's as given is a file that the run's command names, which
  // the page offers: pasted text, and what the draft says is made of a file (InputSnapshot.from).
  it('keeps what a made input is made of, and saves its bytes under its name', async () => {
    const { coordinator, saved } = setup(new ManualEngine());
    const pasted = (await coordinator.enqueue(request)).runId!;
    expect(viewOf(coordinator, pasted)!.snapshot.query.from).toEqual([{ name: 'pasted text', records: 1, leftOut: [], leftOutCount: 0 }]);
    const file = new File(['>chr1\nACGTACGT\n'], 'genome.fa');
    const given = (await coordinator.enqueue({ ...request, subject: { file } })).runId!;
    expect(viewOf(coordinator, given)!.snapshot.subject.from).toBeUndefined();
    await expect(coordinator.exportInput(given, 'subject')).rejects.toThrow('a file of yours as you gave it');
    await coordinator.exportInput(pasted, 'query');
    expect(saved).toEqual([{ fileName: viewOf(coordinator, pasted)!.snapshot.query.name, text: '>q\nACGT\n', type: 'text/plain' }]);
    expect(saved[0]!.fileName).toMatch(/^losat-run1-blastn-\d{8}T\d{6}Z\.query\.fa$/);
  });

  it('reads a chosen file by reference and names the input after it', async () => {
    const { coordinator } = setup(new ManualEngine());
    const file = new File(['>chr1 genome\nACGTACGT\n>chr2\nGG\n'], 'genome.fa');
    const runId = (await coordinator.enqueue({ ...request, subject: { file } })).runId!;
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(snapshot.argv).toEqual(['blastn', '-query', snapshot.query.name, '-subject', 'genome.fa']);
    expect(snapshot.query.name).toMatch(/^losat-run1-blastn-\d{8}T\d{6}Z\.query\.fa$/);
    expect(snapshot.subject.name).toBe('genome.fa');
    expect(snapshot.subject.records).toEqual([
      { id: 'chr1', length: 8 },
      { id: 'chr2', length: 2 },
    ]);
  });

  it('does not queue an input that the index scan cannot read', async () => {
    const { coordinator } = setup(new ManualEngine());
    expect(await coordinator.enqueue({ ...request, query: { text: '>q\n>?2147483647\nA\n' } })).toEqual({
      ok: false,
      message: "Query FASTA: a query record longer than 2147483647 letters is not supported by LOSAT's BLASTN",
    });
    expect(coordinator.state.get().runs).toHaveLength(0);
  });

  it('cancels a queued job without starting it', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    const first = (await coordinator.enqueue(request)).runId!;
    const second = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'running');
    coordinator.cancel(second);
    expect(statusOf(coordinator, second)).toBe('cancelled');
    engine.finish(first);
    await waitFor(coordinator, first, 'completed');
    expect(engine.started).toEqual([first]);
  });

  it('discards a cancelled run and ignores its late events', async () => {
    const engine = new ManualEngine();
    const { coordinator, data, store } = setup(engine);
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'running');
    coordinator.cancel(runId);
    engine.lateFinish?.();
    expect(statusOf(coordinator, runId)).toBe('running');
    engine.finish(runId);
    await waitFor(coordinator, runId, 'cancelled');
    await expect(data.readOutput(runId, 6)).rejects.toThrow();
    expect(store.usage()).toBe(0);
  });

  it('cancels a run while its outputs are being staged, before the engine starts it (WA-RN-01)', async () => {
    // ManualEngine ignores a cancel of a run that it has not started, as the Wasm engine does.
    const engine = new ManualEngine();
    let release!: () => void;
    const staged = new Promise<void>((resolve) => (release = resolve));
    const { coordinator, data, store } = setup(engine, (inner) => {
      const gateway = Object.create(inner) as DataGateway;
      gateway.openRun = async (runId) => {
        const port = await inner.openRun(runId);
        await staged;
        return port;
      };
      return gateway;
    });
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'preparing');
    coordinator.cancel(runId);
    release();
    await waitFor(coordinator, runId, 'cancelled');
    expect(engine.started).toEqual([]);
    await expect(data.readOutput(runId, 6)).rejects.toThrow();
    expect(store.usage()).toBe(0);
    // The next run starts.
    const next = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, next, 'running');
    expect(engine.started).toEqual([next]);
  });

  it('cancels a FakeEngine run between phases', async () => {
    const { coordinator, data } = setup(new FakeEngine({ phaseDelayMs: 20 }));
    const runId = (await coordinator.enqueue(request)).runId!;
    coordinator.cancel(runId);
    await waitFor(coordinator, runId, 'cancelled');
    await expect(data.readOutput(runId, 0)).rejects.toThrow();
  });

  it('marks a failed run and keeps earlier results', async () => {
    const engine = new ManualEngine();
    const { coordinator, data } = setup(engine);
    const first = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'running');
    engine.finish(first);
    await waitFor(coordinator, first, 'completed');
    const second = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, second, 'running');
    engine.fail(second, 'engine error');
    await waitFor(coordinator, second, 'failed');
    expect(viewOf(coordinator, second)!.record.error).toBe('engine error');
    expect(new TextDecoder().decode(await data.readOutput(first, 6))).toBe(`${first}\n`);
  });

  it('fails a run whose records the engine reads differently from the record table', async () => {
    const { coordinator } = setup(new FakeEngine(), (data) => {
      const gateway = Object.create(data) as DataGateway;
      gateway.buildRunInput = async (ids) => {
        const input = await data.buildRunInput(ids);
        return { ...input, records: input.records.map((r) => ({ ...r, id: `${r.id}-stale` })) };
      };
      return gateway;
    });
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'failed');
    expect(viewOf(coordinator, runId)!.record.error).toBe(
      'The query records that the engine read differ from the record table: ' +
        'record 1 is "q-stale" (length 4) in the record table, but the engine read "q" (length 4)',
    );
  });

  it('fails a run with the reason when storage runs out, keeps earlier results and recovers', async () => {
    const { coordinator, store, saved } = setup();
    const first = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, first, 'completed');
    store.setCapacity(store.usage());
    const second = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, second, 'failed');
    expect(viewOf(coordinator, second)!.record.error).toMatch(/^Not enough temporary storage/);
    await coordinator.exportOutput(first, 6);
    expect(saved[0]!.text).toContain(FAKE_MARKER);
    store.setCapacity(Number.POSITIVE_INFINITY);
    const third = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, third, 'completed');
  });

  it('publishes the temporary storage in use', async () => {
    const { coordinator } = setup();
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    await coordinator.refreshStorage();
    const storage = coordinator.state.get().storage!;
    expect(storage.backend).toBe('memory');
    expect(storage.sessionBytes).toBeGreaterThan(0);
    expect(storage.cleanup).toEqual({ state: 'done', removedSessions: 0 });
  });

  it('runs BLASTX as the other programs (line SX)', async () => {
    const { coordinator } = setup();
    const runId = (await coordinator.enqueue({ ...request, program: 'blastx', subject: { text: '>s\nMKLV\n' } })).runId!;
    await waitFor(coordinator, runId, 'completed');
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(snapshot.argv).toEqual(['blastx', '-query', snapshot.query.name, '-subject', snapshot.subject.name]);
    expect(snapshot.subject.name).toMatch(/^losat-run1-blastx-blastx-\d{8}T\d{6}Z\.subject\.fa$/);
    expect(snapshot.subject.records).toEqual([{ id: 's', length: 4 }]);
  });

  it('does not queue a request that the engine rejects', async () => {
    const engine = new ManualEngine();
    engine.validate = async () => ({ ok: false, message: 'bad option' });
    const { coordinator } = setup(engine);
    expect(await coordinator.enqueue(request)).toEqual({ ok: false, message: 'bad option' });
    expect(coordinator.state.get().runs).toHaveLength(0);
  });

  // NCBI's application start (line RG E1w-3b): the engine's start of the argv names the registry
  // file of -conffile and tells a dry run. Enqueue reads the chosen file once (once for a group), and
  // validates and runs with those bytes under the argv's name; the snapshot keeps the file's names,
  // size and SHA-256, and the bytes go when the run ends.
  it("reads the engine's start at enqueue: the -conffile file read once, under the argv's name, and the dry run", async () => {
    const engine = new ManualEngine();
    engine.startResult = startOf({ conffile: 'settings.ini', dryRun: true });
    const { coordinator } = setup(engine);
    const text = '[BLAST]\n';
    const sha256 = createHash('sha256').update(text).digest('hex');
    let reads = 0;
    const chosen = new File([text], 'my settings.ini');
    const read = chosen.arrayBuffer.bind(chosen);
    chosen.arrayBuffer = () => (reads++, read());
    const held = (coordinator as unknown as { conffiles: Map<string, Blob> }).conffiles;
    const runId = (await coordinator.enqueue({ ...request, extra: ['-conffile', 'settings.ini', '-dryrun'], conffile: chosen })).runId!;
    expect(engine.registries.map((registry) => registry?.name)).toEqual(['settings.ini']);
    expect(await engine.registries[0]!.file.text()).toBe(text);
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(snapshot.conffile).toEqual({ name: 'settings.ini', fileName: 'my settings.ini', size: 8, sha256 });
    expect(Object.isFrozen(snapshot.conffile)).toBe(true);
    expect(snapshot.dryRun).toBe(true);
    await waitFor(coordinator, runId, 'running');
    // The run reads the bytes that enqueue read and validated.
    expect(engine.requests[0]!.conffile!.name).toBe('settings.ini');
    expect(engine.requests[0]!.conffile!.file).toBe(engine.registries[0]!.file);
    engine.finish(runId);
    await waitFor(coordinator, runId, 'completed');
    expect(held.size).toBe(0);

    // An argv without -conffile does not read a chosen file, and a run that is not a dry run says nothing of one.
    engine.startResult = startOf({});
    const plain = (await coordinator.enqueue({ ...request, conffile: chosen })).runId!;
    expect(engine.registries.at(-1)).toBeUndefined();
    expect('conffile' in viewOf(coordinator, plain)!.snapshot).toBe(false);
    expect('dryRun' in viewOf(coordinator, plain)!.snapshot).toBe(false);
    await waitFor(coordinator, plain, 'running');
    expect('conffile' in engine.requests[1]!).toBe(false);
    engine.finish(plain);
    await waitFor(coordinator, plain, 'completed');

    // `-conffile=` names no file that the engine can hold: the engine's message stands, without the file.
    engine.startResult = startOf({ conffile: '' });
    const empty = (await coordinator.enqueue({ ...request, extra: ['-conffile='], conffile: chosen })).runId!;
    expect(engine.registries.at(-1)).toBeUndefined();
    expect('conffile' in viewOf(coordinator, empty)!.snapshot).toBe(false);
    await waitFor(coordinator, empty, 'running');
    engine.finish(empty);
    await waitFor(coordinator, empty, 'completed');
    expect(reads).toBe(1);

    // A group reads the file once and shares its bytes; its SHA-256 is counted once. A cancelled run lets them go.
    engine.startResult = startOf({ conffile: 'settings.ini' });
    let digests = 0;
    (coordinator as unknown as { deps: { digest: (bytes: Uint8Array) => Promise<string> } }).deps.digest = async (b) => {
      digests++;
      return sha256Hex(b);
    };
    const group = await coordinator.enqueueAll([
      { ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen },
      { ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen },
      { ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen },
    ]);
    expect(reads).toBe(2);
    expect(digests).toBe(1);
    const ids = group.runIds!;
    const snapshots = ids.map((id) => viewOf(coordinator, id)!.snapshot);
    expect(snapshots.map((s) => s.conffile?.sha256)).toEqual([sha256, sha256, sha256]);
    expect(new Set(ids.map((id) => held.get(id))).size).toBe(1);
    coordinator.cancel(ids[2]!);
    expect(held.has(ids[2]!)).toBe(false);
    for (const id of ids.slice(0, 2)) {
      await waitFor(coordinator, id, 'running');
      engine.finish(id);
      await waitFor(coordinator, id, 'completed');
    }
    // The runs of the group read the one copy.
    expect(engine.requests.at(-2)!.conffile!.file).toBe(engine.requests.at(-1)!.conffile!.file);
    expect(held.size).toBe(0);
  });

  // Review 3 L4: the bytes of a run's registry file go when the run ends however it ends: completed (above),
  // failed, or cancelled while it runs.
  it("lets the run's registry file go when the run fails or is cancelled while it runs", async () => {
    const engine = new ManualEngine();
    engine.startResult = startOf({ conffile: 'settings.ini' });
    const { coordinator } = setup(engine);
    const held = (coordinator as unknown as { conffiles: Map<string, Blob> }).conffiles;
    const chosen = new File(['[BLAST]\n'], 'my settings.ini');
    const failed = (await coordinator.enqueue({ ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen })).runId!;
    await waitFor(coordinator, failed, 'running');
    expect(held.has(failed)).toBe(true);
    engine.fail(failed, 'LOSAT engine error: x');
    await waitFor(coordinator, failed, 'failed');
    expect(held.size).toBe(0);
    const cancelled = (await coordinator.enqueue({ ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen })).runId!;
    await waitFor(coordinator, cancelled, 'running');
    expect(held.has(cancelled)).toBe(true);
    coordinator.cancel(cancelled);
    engine.fail(cancelled, 'cancelled');
    await waitFor(coordinator, cancelled, 'cancelled');
    expect(held.size).toBe(0);
  });

  // Review 3 L2: a chosen file that changed or went after it was chosen cannot be read at "Run LOSAT".
  it('does not queue a request whose -conffile file cannot be read, says so, and holds nothing', async () => {
    const engine = new ManualEngine();
    engine.startResult = startOf({ conffile: 'settings.ini' });
    const { coordinator } = setup(engine);
    const held = (coordinator as unknown as { conffiles: Map<string, Blob> }).conffiles;
    const chosen = new File(['[BLAST]\n'], 'my settings.ini');
    chosen.arrayBuffer = () =>
      Promise.reject(new DOMException('The requested file could not be read, typically due to permission problems that have occurred after a reference to a file was acquired.', 'NotReadableError'));
    expect(await coordinator.enqueue({ ...request, extra: ['-conffile', 'settings.ini'], conffile: chosen })).toEqual({
      ok: false,
      message:
        'The file chosen for -conffile settings.ini (my settings.ini) could not be read: The requested file could not be read, typically due to permission problems that have occurred after a reference to a file was acquired. Choose it again.',
    });
    expect(engine.registries).toEqual([]);
    expect(coordinator.state.get().runs).toHaveLength(0);
    expect(held.size).toBe(0);
  });

  it('does not queue a request whose argv the engine cannot start, with its message', async () => {
    const engine = new ManualEngine();
    engine.startResult = { ok: false, message: 'the engine cannot start this argv' };
    const { coordinator } = setup(engine);
    expect(await coordinator.enqueue({ ...request, extra: ['-unstartable'] })).toEqual({ ok: false, message: 'the engine cannot start this argv' });
    expect(engine.registries).toEqual([]);
    expect(coordinator.state.get().runs).toHaveLength(0);
  });

  // ABI v2 stream 8 (line RG E1w-3b): the file of -logfile comes with the engine's answer, also when
  // the run fails, and the run keeps it; a cancelled run keeps none.
  // line WQA's sweep, D7, 2026-11-19: the results name the records of a -parse_deflines run as its outputs do.
  it('keeps the IDs that the outputs of a run show, where the engine read them from the deflines', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    const plain = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, plain, 'running');
    engine.finish(plain);
    await waitFor(coordinator, plain, 'completed');
    expect('shownIds' in viewOf(coordinator, plain)!.record).toBe(false);
    engine.shownIds = { query: ['AB000001.1'], subject: ['NP_000006.1'] };
    const parsed = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, parsed, 'running');
    engine.finish(parsed);
    await waitFor(coordinator, parsed, 'completed');
    expect(viewOf(coordinator, parsed)!.record.shownIds).toEqual({ query: ['AB000001.1'], subject: ['NP_000006.1'] });
  });

  // line WQA's sweep, D7, 2026-11-20: where the run's messages went is the engine's word, kept with the completed run.
  it("keeps where the run's messages went, as the engine said", async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    for (const messages of ['diagnostics', 'logfile', 'nowhere'] as const) {
      engine.messages = messages;
      const runId = (await coordinator.enqueue(request)).runId!;
      await waitFor(coordinator, runId, 'running');
      engine.finish(runId);
      await waitFor(coordinator, runId, 'completed');
      expect(viewOf(coordinator, runId)!.record.messages).toBe(messages);
    }
  });

  // Line WQA's D14, 2026-11-20: an argv with a help or version word prints the command line's text and searches nothing,
  // as the CLI does. The run reads no input (none is chosen), keeps the text, and saves it under the run's names.
  it("runs an argv that prints the command line's help or version without inputs, and keeps and saves its text", async () => {
    const engine = new ManualEngine();
    engine.startResult = startOf({ prints: 'version' });
    engine.printed = new TextEncoder().encode('blastn: LOSAT 0.2.0\n');
    const { coordinator, saved } = setup(engine);
    const runId = (await coordinator.enqueue({ ...request, query: { text: '' }, subject: { text: '' }, extra: ['-version'] })).runId!;
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(snapshot.prints).toBe('version');
    expect(snapshot.query).toMatchObject({ bytes: new Uint8Array(0), records: [], revisionIds: [] });
    expect(snapshot.subject).toMatchObject({ bytes: new Uint8Array(0), records: [], revisionIds: [] });
    await waitFor(coordinator, runId, 'running');
    engine.finish(runId);
    await waitFor(coordinator, runId, 'completed');
    expect(new TextDecoder().decode(viewOf(coordinator, runId)!.record.printed)).toBe('blastn: LOSAT 0.2.0\n');
    await coordinator.exportPrinted(runId);
    expect(saved).toEqual([{ fileName: 'losat-run1-blastn-19700101T000001Z.version.txt', text: 'blastn: LOSAT 0.2.0\n', type: 'text/plain' }]);
    // NCBI's XML and JSON are saved as what they are (the bytes unchanged).
    for (const [word, kind, text, name, type] of [
      ['-xmlhelp', 'help', '<?xml version="1.0"?>\n<ncbi_application/>\n', 'help.xml', 'application/xml'],
      ['-version-full-json', 'version', '{\n  "appname": "blastn"\n}\n', 'version.json', 'application/json'],
    ] as const) {
      engine.printed = new TextEncoder().encode(text);
      engine.startResult = startOf({ prints: kind });
      const other = (await coordinator.enqueue({ ...request, query: { text: '' }, subject: { text: '' }, extra: [word] })).runId!;
      await waitFor(coordinator, other, 'running');
      engine.finish(other);
      await waitFor(coordinator, other, 'completed');
      saved.length = 0;
      await coordinator.exportPrinted(other);
      expect(saved).toEqual([{ fileName: expect.stringMatching(new RegExp(`\\.${name.replace('.', '\\.')}$`)), text, type }]);
    }

    // A search has no printed text, and says nothing of one.
    engine.startResult = startOf({});
    delete engine.printed;
    const search = (await coordinator.enqueue(request)).runId!;
    expect('prints' in viewOf(coordinator, search)!.snapshot).toBe(false);
    await waitFor(coordinator, search, 'running');
    engine.finish(search);
    await waitFor(coordinator, search, 'completed');
    expect('printed' in viewOf(coordinator, search)!.record).toBe(false);
    await expect(coordinator.exportPrinted(search)).rejects.toThrow('the run printed no text');
  });

  it('keeps the files that a run made besides its outputs, of a completed and of a failed run, and none of a cancelled one', async () => {
    const engine = new ManualEngine();
    const { coordinator, saved } = setup(engine);
    const completed = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, completed, 'running');
    engine.finish(completed, [logFile('Warning: [tblastx] Examining 5 or more matches is recommended\n')]);
    await waitFor(coordinator, completed, 'completed');
    expect(viewOf(coordinator, completed)!.record.runFiles).toEqual([logFile('Warning: [tblastx] Examining 5 or more matches is recommended\n')]);

    const failed = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, failed, 'running');
    engine.fail(failed, 'LOSAT engine error: The average subject length is too short', [logFile('LOSAT engine error: The average subject length is too short\n')]);
    await waitFor(coordinator, failed, 'failed');
    expect(viewOf(coordinator, failed)!.record).toMatchObject({
      error: 'LOSAT engine error: The average subject length is too short',
      runFiles: [logFile('LOSAT engine error: The average subject length is too short\n')],
    });

    const plain = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, plain, 'running');
    engine.fail(plain, 'no files');
    await waitFor(coordinator, plain, 'failed');
    expect('runFiles' in viewOf(coordinator, plain)!.record).toBe(false);

    const cancelled = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, cancelled, 'running');
    coordinator.cancel(cancelled);
    engine.finish(cancelled, [logFile('late\n')]);
    await waitFor(coordinator, cancelled, 'cancelled');
    expect('runFiles' in viewOf(coordinator, cancelled)!.record).toBe(false);

    // A file is saved byte for byte under the run's names (a program without a default task here), a failed run's too.
    await coordinator.exportRunFile(completed, 0);
    await coordinator.exportRunFile(failed, 0);
    expect(saved).toEqual([
      { fileName: 'losat-run1-blastn-19700101T000001Z.logfile.run.log', text: 'Warning: [tblastx] Examining 5 or more matches is recommended\n', type: 'text/plain' },
      { fileName: 'losat-run2-blastn-19700101T000001Z.logfile.run.log', text: 'LOSAT engine error: The average subject length is too short\n', type: 'text/plain' },
    ]);
    await expect(coordinator.exportRunFile(plain, 0)).rejects.toThrow('the run has no such file');
  });

  // Line WQA's sweep, D3 (2026-11-19): the warnings that the command line writes before its error
  // (stream 3) arrive at the data layer before the engine's `fail`, and the failed run keeps them, also
  // where the rejection comes first; a cancelled run keeps none, and a failed run without warnings has
  // none.
  it('keeps the warnings that a failed run wrote before its error, and none of a cancelled run', async () => {
    const engine = new ManualEngine();
    const { coordinator, data } = setup(engine);
    const warning = 'Warning: [blastn] Subject_1 empty record: Subject sequence contains no data\n';
    const failed = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, failed, 'running');
    engine.fail(failed, 'LOSAT engine error: The average subject length is too short', undefined, warning);
    await waitFor(coordinator, failed, 'failed');
    expect(viewOf(coordinator, failed)!.record).toMatchObject({ error: 'LOSAT engine error: The average subject length is too short', diagnostics: warning });
    expect('runFiles' in viewOf(coordinator, failed)!.record).toBe(false);
    // The run kept no output (none has bytes), so it is dropped as before version 7: only its record holds the warnings.
    await expect(data.readDiagnostics(failed)).rejects.toThrow();
    await expect(data.keepFailedRun(failed)).resolves.toBeUndefined();
    expect(viewOf(coordinator, failed)!.result).toBeUndefined();

    // The rejection before the chunks and the `fail`: the data layer waits for them.
    const late = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, late, 'running');
    engine.fail(late, 'LOSAT engine error: The average subject length is too short', undefined, warning, true);
    await waitFor(coordinator, late, 'failed');
    expect(viewOf(coordinator, late)!.record).toMatchObject({ diagnostics: warning });

    // A failed run whose output ended with `fail` but holds no warning keeps none.
    const silent = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, silent, 'running');
    engine.fail(silent, 'LOSAT engine error: Empty CBlastQueryVector', undefined, '');
    await waitFor(coordinator, silent, 'failed');
    expect('diagnostics' in viewOf(coordinator, silent)!.record).toBe(false);

    const plain = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, plain, 'running');
    engine.fail(plain, 'LOSAT engine error: Empty CBlastQueryVector');
    await waitFor(coordinator, plain, 'failed');
    expect('diagnostics' in viewOf(coordinator, plain)!.record).toBe(false);

    const cancelled = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, cancelled, 'running');
    coordinator.cancel(cancelled);
    engine.fail(cancelled, 'cancelled', undefined, warning);
    await waitFor(coordinator, cancelled, 'cancelled');
    expect('diagnostics' in viewOf(coordinator, cancelled)!.record).toBe(false);
  });

  // Line WQA's sweep, D15 (2026-11-20): a failed run keeps the output that the CLI leaves before its error, the bytes
  // written to the output of the argv, and the engine drops the other outputs (stream 9, place `run`); the data layer
  // commits the run as a completed run's, so that the output is read and exported as a completed run's is. A
  // cancelled run keeps nothing.
  it('keeps the output that a failed run wrote before its error, reads and exports it', async () => {
    const engine = new ManualEngine();
    const { coordinator, data, saved } = setup(engine);
    const failed = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, failed, 'running');
    engine.fail(failed, 'LOSAT engine error: The average subject length is too short', undefined, 'Warning: one\n', false, { output: 0, text: 'BLASTN 2.17.0+\n\nDatabase: x\n' });
    await waitFor(coordinator, failed, 'failed');
    const view = viewOf(coordinator, failed)!;
    expect(view.result).toMatchObject({ byteLengths: { 0: 28, 6: 0, 7: 0 }, hitCount: 0 });
    expect(view.result!.failed).toEqual([
      { output: 6, error: 'run' },
      { output: 7, error: 'run' },
    ]);
    expect(await coordinator.readOutput(failed, 0)).toBe('BLASTN 2.17.0+\n\nDatabase: x\n');
    expect(new TextDecoder().decode(await data.readOutput(failed, 6))).toBe('');
    await coordinator.exportOutput(failed, 0);
    expect(saved.at(-1)!.text).toBe('BLASTN 2.17.0+\n\nDatabase: x\n');
    await expect(coordinator.readHits(failed)).rejects.toThrow('run is not completed');

    const cancelled = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, cancelled, 'running');
    coordinator.cancel(cancelled);
    engine.fail(cancelled, 'cancelled', undefined, 'Warning: one\n', false, { output: 0, text: 'kept?' });
    await waitFor(coordinator, cancelled, 'cancelled');
    expect(viewOf(coordinator, cancelled)!.result).toBeUndefined();
    await expect(coordinator.readOutput(cancelled, 0)).rejects.toThrow('run has no outputs');
  });

  it('gives the FakeEngine\'s failed run its marked warning', async () => {
    const { coordinator } = setup();
    const runId = (await coordinator.enqueue({ ...request, subject: { text: '' } })).runId!;
    await waitFor(coordinator, runId, 'failed');
    expect(viewOf(coordinator, runId)!.record).toMatchObject({ error: expect.stringContaining('No sequence input was provided'), diagnostics: `Warning: ${FAKE_MARKER}\n` });
    expect(runOutputFailedOf(new Error('x'))).toBe(false);
  });

  it('exports a completed output byte for byte under a name with the task and the time of the run (WA-RS-05)', async () => {
    const { coordinator, data, saved } = setup();
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    await coordinator.exportOutput(runId, 7);
    expect(saved).toHaveLength(1);
    // The engine's default task, and the UTC time of "Run LOSAT" (now() is 1000 ms).
    expect(saved[0]!.fileName).toBe('losat-run1-blastn-megablast-19700101T000001Z.outfmt7.txt');
    expect(saved[0]!.text).toBe(new TextDecoder().decode(await data.readOutput(runId, 7)));
    expect(saved[0]!.text).toContain(FAKE_MARKER);
    const short = (await coordinator.enqueue({ ...request, parameters: [['-task', 'blastn-short']] })).runId!;
    await waitFor(coordinator, short, 'completed');
    await coordinator.exportOutput(short, 0);
    expect(saved[1]!.fileName).toBe('losat-run2-blastn-blastn-short-19700101T000001Z.outfmt0.txt');
    // A program without tasks.
    const tblastx = (await coordinator.enqueue({ ...request, program: 'tblastx' })).runId!;
    await waitFor(coordinator, tblastx, 'completed');
    await coordinator.exportOutput(tblastx, 6);
    expect(saved[2]!.fileName).toBe('losat-run3-tblastx-19700101T000001Z.outfmt6.txt');
  });

  // The output of the argv's -outfmt with fields (WA2-OP-03): run, stored and exported besides the
  // formats of every run; the engine refuses the formats that a run does not write.
  it("runs, stores and exports the output of the argv's -outfmt with fields", async () => {
    const { coordinator, data, saved } = setup();
    const runId = (await coordinator.enqueue({ ...request, extra: ['-outfmt', '7 std qcovs'] })).runId!;
    await waitFor(coordinator, runId, 'completed');
    const view = coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)!;
    expect(Object.keys(view.result!.byteLengths)).toEqual(['0', '4', '6', '7']);
    const custom = await coordinator.readOutput(runId, CUSTOM_OUTPUT);
    expect(custom).toContain('# -outfmt 7 std qcovs');
    await coordinator.exportOutput(runId, CUSTOM_OUTPUT);
    expect(saved[0]!.fileName).toBe('losat-run1-blastn-megablast-19700101T000001Z.outfmt7-custom.txt');
    expect(saved[0]!.text).toBe(custom);
    expect(saved[0]!.type).toBe('text/plain');
    // Its formats are those of a run without -outfmt.
    const plain = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, plain, 'completed');
    for (const format of [0, 6, 7] as const) {
      expect(await coordinator.readOutput(runId, format)).toBe(await coordinator.readOutput(plain, format));
    }
    await expect(data.readOutput(plain, CUSTOM_OUTPUT)).rejects.toThrow('does not write output 4');
  });

  // Every other format that the CLI writes is the output of the argv's -outfmt too (W11), saved
  // with its format's extension and media type.
  it("runs, stores and exports the output of an -outfmt besides 0, 6 and 7", async () => {
    const { coordinator, saved } = setup();
    const runId = (await coordinator.enqueue({ ...request, extra: ['-outfmt', '5'] })).runId!;
    await waitFor(coordinator, runId, 'completed');
    const view = coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)!;
    expect(Object.keys(view.result!.byteLengths)).toEqual(['0', '4', '6', '7']);
    expect(await coordinator.readOutput(runId, CUSTOM_OUTPUT)).toContain('# -outfmt 5');
    await coordinator.exportOutput(runId, CUSTOM_OUTPUT);
    expect(saved[0]!.fileName).toBe('losat-run1-blastn-megablast-19700101T000001Z.outfmt5.xml');
    expect(saved[0]!.type).toBe('application/xml');
    await coordinator.exportOutput(runId, 0);
    expect(saved[1]!.type).toBe('text/plain');
    // An output is opened without being read (a binary one, which the Outputs do not show).
    await expect(coordinator.checkOutput(runId, CUSTOM_OUTPUT)).resolves.toBeUndefined();
    // The files of 13 and 14 are one ZIP archive.
    const files = (await coordinator.enqueue({ ...request, extra: ['-outfmt', '13'] })).runId!;
    await waitFor(coordinator, files, 'completed');
    await coordinator.exportOutput(files, CUSTOM_OUTPUT);
    expect(saved[2]!.fileName).toBe('losat-run2-blastn-megablast-19700101T000001Z.outfmt13.zip');
    expect(saved[2]!.type).toBe('application/zip');
  });

  it('queues a group in order, numbered without gaps, and cancels the group together', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    const result = await coordinator.enqueueAll([request, { ...request, query: { text: '>q2\nACGT\n' } }, request]);
    expect(result.ok).toBe(true);
    const ids = result.runIds!;
    const runs = coordinator.state.get().runs;
    expect(runs.map((run) => run.snapshot.number)).toEqual([1, 2, 3]);
    expect(runs.map((run) => run.snapshot.group?.position)).toEqual([1, 2, 3]);
    expect(new Set(runs.map((run) => run.snapshot.group?.groupId)).size).toBe(1);
    expect(runs.every((run) => run.snapshot.group?.size === 3)).toBe(true);
    await waitFor(coordinator, ids[0]!, 'running');
    coordinator.cancelGroup(runs[0]!.snapshot.group!.groupId);
    expect(statusOf(coordinator, ids[1]!)).toBe('cancelled');
    expect(statusOf(coordinator, ids[2]!)).toBe('cancelled');
    engine.finish(ids[0]!);
    await waitFor(coordinator, ids[0]!, 'cancelled');
    // The cancelled runs of the group never started.
    expect(engine.started).toEqual([ids[0]]);
  });

  it('keeps the Job Title trimmed in each snapshot of a group, outside the argv, and none when it is empty', async () => {
    const { coordinator } = setup(new ManualEngine());
    const titled = { ...request, title: '  Plasmid screen  ' };
    const group = await coordinator.enqueueAll([titled, titled]);
    const snapshots = group.runIds!.map((id) => viewOf(coordinator, id)!.snapshot);
    expect(snapshots.map((snapshot) => snapshot.title)).toEqual(['Plasmid screen', 'Plasmid screen']);
    // Each run of the group has its own number, so its own names (the same time).
    expect(snapshots.map((snapshot) => snapshot.argv.slice(1, 5))).toEqual(
      snapshots.map((snapshot) => ['-query', snapshot.query.name, '-subject', snapshot.subject.name]),
    );
    expect(snapshots.map((snapshot) => snapshot.query.name)).toEqual([
      expect.stringMatching(/^losat-run1-blastn-.*\.query\.fa$/),
      expect.stringMatching(/^losat-run2-blastn-.*\.query\.fa$/),
    ]);
    for (const title of [undefined, '', '   ']) {
      const runId = (await coordinator.enqueue(title === undefined ? request : { ...request, title })).runId!;
      expect('title' in viewOf(coordinator, runId)!.snapshot).toBe(false);
    }
  });

  it('queues nothing of a group when one request is invalid, and keeps the numbering', async () => {
    const engine = new ManualEngine();
    engine.validate = async (argv?: readonly string[]): Promise<ValidationResult> =>
      argv?.includes('-word_size') === true ? { ok: false, message: 'bad word size' } : { ok: true, dataLoaders: true };
    const { coordinator } = setup(engine);
    const bad = { ...request, parameters: [['-word_size', '1']] as const };
    expect(await coordinator.enqueueAll([request, bad])).toEqual({ ok: false, message: 'bad word size' });
    expect(coordinator.state.get().runs).toHaveLength(0);
    const runId = (await coordinator.enqueue(request)).runId!;
    expect(viewOf(coordinator, runId)!.snapshot.number).toBe(1);
    expect(viewOf(coordinator, runId)!.snapshot.group).toBeUndefined();
  });

  // DECISIONS 44 (line WQA's D22): every input that the page makes is named with the run's own prefix and its role, in the argv that
  // runs and in the file that the page saves; the user's file as given keeps the user's name (D21).
  it('names a made input after its run: the prefix of the run\'s outputs and the role, each run of a group its own', async () => {
    const { coordinator, data } = setup(new ManualEngine());
    const made = async (text: string) => ({ dataset: { revisionIds: [(await data.indexSource((await data.addSource(new File([text], 'x.fa'))).sourceId, { program: 'blastn', role: 'query', dataLoaders: true })).revisionId], from: [{ name: 'x.fa', records: 1, leftOut: [], leftOutCount: 1 }] } });
    const given = { file: new File(['>g\nACGT\n'], 'My file.fa') };
    const group = await coordinator.enqueueAll([
      { ...request, query: await made('>a\nACGT\n'), subject: given, parameters: [['-task', 'blastn']] },
      { ...request, query: await made('>b\nACGT\n'), subject: given, parameters: [['-task', 'blastn']] },
    ]);
    const snapshots = group.runIds!.map((id) => viewOf(coordinator, id)!.snapshot);
    expect(snapshots.map((snapshot) => snapshot.query.name)).toEqual([
      'losat-run1-blastn-blastn-19700101T000001Z.query.fa',
      'losat-run2-blastn-blastn-19700101T000001Z.query.fa',
    ]);
    // The user's own name is kept as it is, and the argv names what the run holds.
    expect(snapshots.map((snapshot) => snapshot.subject.name)).toEqual(['My file.fa', 'My file.fa']);
    for (const snapshot of snapshots) expect(snapshot.argv.slice(1, 5)).toEqual(['-query', snapshot.query.name, '-subject', 'My file.fa']);
  });

  it('asks the engine to describe a program once for the requests of a press', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    await coordinator.enqueueAll([request, request, request, request]);
    expect(engine.described).toBe(1);
    // A press of its own asks again (a description is not kept between presses: the engine may change).
    await coordinator.enqueueAll([request]);
    expect(engine.described).toBe(2);
  });

  it('numbers the runs of two enqueues made at once one after the other, and names them by those numbers', async () => {
    const { coordinator } = setup(new ManualEngine());
    const [first, second] = await Promise.all([coordinator.enqueueAll([request, request]), coordinator.enqueueAll([request])]);
    const names = [...first.runIds!, ...second.runIds!].map((id) => viewOf(coordinator, id)!.snapshot).map((snapshot) => [snapshot.number, snapshot.query.name.slice(0, 15)]);
    expect(names).toEqual([
      [1, 'losat-run1-blas'],
      [2, 'losat-run2-blas'],
      [3, 'losat-run3-blas'],
    ]);
  });

  it('drops the sources that it makes of text and files once the snapshots have their inputs (SURVEY #15)', async () => {
    const { coordinator, data } = setup();
    const inner = data as unknown as { sources: Map<string, File>; revisions: Map<string, unknown> };
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    expect([inner.sources.size, inner.revisions.size]).toEqual([0, 0]);
    expect(new TextDecoder().decode(viewOf(coordinator, runId)!.snapshot.query.bytes)).toBe('>q\nACGT\n');
    expect(await coordinator.readOutput(runId, 6)).toContain('q\ts\t');
  });

  it('keeps nothing of a request whose input cannot be read: no source and no cached run input', async () => {
    const { coordinator, data } = setup();
    const inner = data as unknown as { sources: Map<string, File>; revisions: Map<string, unknown> };
    const runInputs = (coordinator as unknown as { runInputs: Map<string, unknown> }).runInputs;
    const source = await data.addSource(new File(['>q\nACGT\n'], 'q.fa'));
    const query = (await data.indexSource(source.sourceId, { program: 'blastn', role: 'query', dataLoaders: true })).revisionId;
    const result = await coordinator.enqueueAll([
      // A record over 2147483647 letters, which LOSAT's reader refuses (a reader's error ends the records
      // instead: line WQA's D10, 2026-11-19).
      { ...request, query: { dataset: { name: 'q.fa', revisionIds: [query] } }, subject: { text: '>s\nAC\n>?2147483647\nAC\n' } },
    ]);
    expect(result.ok).toBe(false);
    expect(runInputs.size).toBe(0);
    // Only the caller's source is left: the subject's is dropped.
    expect([...inner.sources.keys()]).toEqual([source.sourceId]);
  });

  it('runs dataset inputs as the revisions say, and shares the bytes of the same revisions', async () => {
    const { coordinator, data } = setup();
    const add = async (name: string, text: string) => {
      const source = await data.addSource(new File([text], name));
      return (await data.indexSource(source.sourceId, { program: 'blastn', role: 'query', dataLoaders: true })).revisionId;
    };
    const subject = await add('s.fa', '>s1\nACGT\n>s2\nGGCC\n');
    const query = await add('q.fa', '>q\nACGT\n');
    const revised = (await data.reviseDataset(subject, [0])).revisionId;
    const result = await coordinator.enqueueAll([
      { ...request, query: { dataset: { name: 'q.fa', revisionIds: [query] } }, subject: { dataset: { name: 's.fa', revisionIds: [subject] } } },
      { ...request, query: { dataset: { name: 'q.fa', revisionIds: [query] } }, subject: { dataset: { name: 's.fa', revisionIds: [subject] } } },
      { ...request, query: { dataset: { name: 'q.fa', revisionIds: [query] } }, subject: { dataset: { name: 's2.fa', revisionIds: [revised] } } },
    ]);
    const [a, b, c] = result.runIds!.map((id) => viewOf(coordinator, id)!.snapshot);
    expect(a!.argv).toEqual(['blastn', '-query', 'q.fa', '-subject', 's.fa']);
    expect(a!.subject.revisionIds).toEqual([subject]);
    expect(a!.subject.bytes).toBe(b!.subject.bytes);
    expect(new TextDecoder().decode(c!.subject.bytes)).toBe('>s2\nGGCC\n');
    expect(c!.subject.records).toEqual([{ id: 's2', length: 4 }]);
  });
});
