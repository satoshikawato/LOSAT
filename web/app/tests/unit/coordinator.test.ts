import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import { Coordinator, type SearchRequest } from '../../src/application/coordinator';
import { CUSTOM_OUTPUT } from '../../src/domain/output-format';
import type { RunStatus } from '../../src/domain/run';
import { sha256Hex } from '../../src/infra/browser/platform';
import { DataService } from '../../src/infra/data/data-service';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { FakeEngine, FAKE_MARKER } from '../../src/infra/fake/fake-engine';
import { FakeInputChecker, FakeScanner } from '../../src/infra/fake/fake-fasta';
import { RunOutputWriter } from '../../src/infra/run-output/writer';
import type { DataGateway } from '../../src/ports/data';
import type { Downloader } from '../../src/ports/download';
import type {
  EngineGateway,
  EnginePhase,
  EngineRunRequest,
  RuntimeInfo,
  ValidationResult,
} from '../../src/ports/engine';
import { DIAGNOSTICS_STREAM } from '../../src/ports/run-output';

const request: SearchRequest = {
  program: 'blastn',
  query: { text: '>q\nACGT\n' },
  subject: { text: '>s\nACGT\n' },
  parameters: [],
  requestedThreads: 'auto',
};

function setup(engine: EngineGateway = new FakeEngine(), wrap: (data: DataService) => DataGateway = (d) => d) {
  const store = new MemoryBlockStore();
  let token = 0;
  const data = new DataService({
    store,
    scanner: new FakeScanner(),
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
  });
  return { coordinator, data, store, saved };
}

const viewOf = (coordinator: Coordinator, runId: string) =>
  coordinator.state.get().runs.find((run) => run.snapshot.runId === runId);
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
  private readonly pending = new Map<string, { resolve: () => void; reject: (e: Error) => void }>();
  lateFinish?: () => void;

  async describe() {
    return { program: 'blastn' as const, formats: [6] as const, parameters: [] };
  }
  async validate(): Promise<ValidationResult> {
    return { ok: true };
  }
  run(req: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    this.started.push(req.runId);
    onPhase('running');
    return new Promise((resolve, reject) => {
      this.pending.set(req.runId, {
        resolve: () => {
          const writer = new RunOutputWriter(output);
          writer.write(6, new TextEncoder().encode(`${req.runId}\n`));
          writer.write(DIAGNOSTICS_STREAM, new TextEncoder().encode(`Warning: ${req.runId}\n`));
          writer.end();
          resolve({ path: 'fake', threads: 1, engineBuild: 'manual' });
        },
        reject,
      });
      // Simulates a phase event that arrives after the run was cancelled.
      this.lateFinish = () => onPhase('finalizing');
    });
  }
  cancel(): void {}
  finish(runId: string) {
    this.pending.get(runId)?.resolve();
  }
  fail(runId: string, message: string) {
    this.pending.get(runId)?.reject(new Error(message));
  }
}

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
    expect(snapshot.query.name).toBe('query.fa');
    expect(snapshot.subject.sha256).toBe(createHash('sha256').update('>s\nACGT\n').digest('hex'));
    expect(snapshot.subject.records).toEqual([{ id: 's', length: 4 }]);
    expect(snapshot.subject.revisionIds).toHaveLength(1);
  });

  it('reads a chosen file by reference and names the input after it', async () => {
    const { coordinator } = setup(new ManualEngine());
    const file = new File(['>chr1 genome\nACGTACGT\n>chr2\nGG\n'], 'genome.fa');
    const runId = (await coordinator.enqueue({ ...request, subject: { file } })).runId!;
    const snapshot = viewOf(coordinator, runId)!.snapshot;
    expect(snapshot.argv).toEqual(['blastn', '-query', 'query.fa', '-subject', 'genome.fa']);
    expect(snapshot.subject.name).toBe('genome.fa');
    expect(snapshot.subject.records).toEqual([
      { id: 'chr1', length: 8 },
      { id: 'chr2', length: 2 },
    ]);
  });

  it('does not queue an input that the index scan cannot read', async () => {
    const { coordinator } = setup(new ManualEngine());
    expect(await coordinator.enqueue({ ...request, query: { text: '>q\nACGT\n>?10\n' } })).toEqual({
      ok: false,
      message:
        "Query FASTA: line 3 is a gap line ('>?'), whose residues have no bytes in the input for LOSAT's index to locate; gap lines are not supported in the browser",
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
    expect(snapshot.argv).toEqual(['blastx', '-query', 'query.fa', '-subject', 'subject.fa']);
    expect(snapshot.subject.records).toEqual([{ id: 's', length: 4 }]);
  });

  it('does not queue a request that the engine rejects', async () => {
    const engine = new ManualEngine();
    engine.validate = async () => ({ ok: false, message: 'bad option' });
    const { coordinator } = setup(engine);
    expect(await coordinator.enqueue(request)).toEqual({ ok: false, message: 'bad option' });
    expect(coordinator.state.get().runs).toHaveLength(0);
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
    expect(snapshots[0]!.argv).toEqual(['blastn', '-query', 'query.fa', '-subject', 'subject.fa']);
    for (const title of [undefined, '', '   ']) {
      const runId = (await coordinator.enqueue(title === undefined ? request : { ...request, title })).runId!;
      expect('title' in viewOf(coordinator, runId)!.snapshot).toBe(false);
    }
  });

  it('queues nothing of a group when one request is invalid, and keeps the numbering', async () => {
    const engine = new ManualEngine();
    engine.validate = async (argv?: readonly string[]): Promise<ValidationResult> =>
      argv?.includes('-word_size') === true ? { ok: false, message: 'bad word size' } : { ok: true };
    const { coordinator } = setup(engine);
    const bad = { ...request, parameters: [['-word_size', '1']] as const };
    expect(await coordinator.enqueueAll([request, bad])).toEqual({ ok: false, message: 'bad word size' });
    expect(coordinator.state.get().runs).toHaveLength(0);
    const runId = (await coordinator.enqueue(request)).runId!;
    expect(viewOf(coordinator, runId)!.snapshot.number).toBe(1);
    expect(viewOf(coordinator, runId)!.snapshot.group).toBeUndefined();
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
    const query = (await data.indexSource(source.sourceId, 1)).revisionId;
    const result = await coordinator.enqueueAll([
      { ...request, query: { dataset: { name: 'q.fa', revisionIds: [query] } }, subject: { text: '$$$$\n' } },
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
      return (await data.indexSource(source.sourceId, 1)).revisionId;
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
