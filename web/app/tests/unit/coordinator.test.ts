import { createHash } from 'node:crypto';
import { describe, expect, it } from 'vitest';
import { Coordinator, type SearchRequest } from '../../src/application/coordinator';
import type { FastaParserKind } from '../../src/domain/dataset';
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
import { memoryDownloader } from './support/memory-downloader';

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
  options: { readonly exportRangeBytes?: number } = {},
) {
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
  const saved: Array<{ fileName: string; text: string; bytes: Uint8Array; blocks: number }> = [];
  const downloader: Downloader = memoryDownloader((file) =>
    saved.push({ fileName: file.name, text: new TextDecoder().decode(file.bytes), bytes: file.bytes, blocks: file.blocks }),
  );
  let id = 0;
  const coordinator = new Coordinator({
    engine,
    data: wrap(data),
    downloader,
    now: () => 1000,
    newRunId: () => `run-${++id}`,
    ...options,
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
    const message = "CFastaReader: Near line 2, there's a line that doesn't look like plausible data, but it's not marked as defline or comment.";
    const { coordinator } = setup(new ManualEngine(), (data) =>
      Object.assign(data, {
        indexSource: async () => {
          throw new Error(message);
        },
      }),
    );
    expect(await coordinator.enqueue({ ...request, query: { text: '>q\n1234567890\n' } })).toEqual({
      ok: false,
      message: `Query FASTA: ${message}`,
    });
    expect(coordinator.state.get().runs).toHaveLength(0);
  });

  it("indexes each text or file input with the reader of the program's kind for its role", async () => {
    const parsers: string[] = [];
    const { coordinator } = setup(new ManualEngine(), (data) => {
      const indexSource = data.indexSource.bind(data);
      return Object.assign(data, {
        indexSource: async (sourceId: string, parser: FastaParserKind) => {
          parsers.push(String(parser));
          return indexSource(sourceId, parser);
        },
      });
    });
    await coordinator.enqueue({ ...request, program: 'tblastn', query: { text: '>q\nMKV*LL\n' }, subject: { text: '>s\nACGU\n' } });
    await coordinator.enqueue({ ...request, program: 'blastp' });
    expect(parsers).toEqual(['2', '1', '2', '2']);
    const [tblastn] = coordinator.state.get().runs;
    // The protein query keeps L and *; the nucleotide subject stores U as T.
    expect(tblastn!.snapshot.query.records).toEqual([{ id: 'q', length: 6 }]);
    expect(tblastn!.snapshot.subject.records).toEqual([{ id: 's', length: 4 }]);
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

  it('refuses BLASTX as ABI v2 does until session SX', async () => {
    const { coordinator } = setup();
    expect(await coordinator.enqueue({ ...request, program: 'blastx' })).toEqual({
      ok: false,
      message: 'blastx is not available in LOSAT Web ABI v2 yet',
    });
  });

  it('does not queue a request that the engine rejects', async () => {
    const engine = new ManualEngine();
    engine.validate = async () => ({ ok: false, message: 'bad option' });
    const { coordinator } = setup(engine);
    expect(await coordinator.enqueue(request)).toEqual({ ok: false, message: 'bad option' });
    expect(coordinator.state.get().runs).toHaveLength(0);
  });

  it('exports a completed output byte for byte under a descriptive name', async () => {
    const { coordinator, saved } = setup();
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    await coordinator.exportOutput(runId, 7);
    expect(saved).toHaveLength(1);
    expect(saved[0]!.fileName).toBe('losat-run1-blastn.outfmt7.txt');
    expect(saved[0]!.text).toContain(FAKE_MARKER);
  });

  it('exports the whole stored output byte for byte, read in bounded ranges and written in order', async () => {
    const { coordinator, data, saved } = setup(new FakeEngine(), (d) => d, { exportRangeBytes: 16 });
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    for (const format of [0, 6, 7] as const) {
      saved.length = 0;
      await coordinator.exportOutput(runId, format);
      const stored = await data.readOutput(runId, format);
      expect(stored.length).toBeGreaterThan(32);
      expect(saved.map((file) => file.fileName)).toEqual([`losat-run1-blastn.outfmt${format}.txt`]);
      expect(saved[0]!.bytes).toEqual(stored);
      expect(saved[0]!.blocks).toBe(Math.ceil(stored.length / 16));
    }
  });

  it('saves nothing when a read of the stored output fails', async () => {
    let reads = 0;
    const { coordinator, saved } = setup(
      new FakeEngine(),
      (data) => {
        const gateway = Object.create(data) as DataGateway;
        gateway.readOutputRange = async (...args) => {
          if (++reads === 2) throw new Error('the storage is gone');
          return data.readOutputRange(...args);
        };
        return gateway;
      },
      { exportRangeBytes: 16 },
    );
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'completed');
    await expect(coordinator.exportOutput(runId, 6)).rejects.toThrow('the storage is gone');
    expect(saved).toEqual([]);
    await coordinator.exportOutput(runId, 6);
    expect(saved[0]!.text).toContain(FAKE_MARKER);
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
