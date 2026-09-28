import { describe, expect, it } from 'vitest';
import { Coordinator, type SearchRequest } from '../../src/application/coordinator';
import type { RunStatus } from '../../src/domain/run';
import { FakeEngine, FAKE_MARKER } from '../../src/infra/fake/fake-engine';
import { MemoryDataGateway } from '../../src/infra/memory/memory-data-gateway';
import type { Downloader } from '../../src/ports/download';
import type {
  EngineGateway,
  EnginePhase,
  EngineRunRequest,
  RunSink,
  RuntimeInfo,
  ValidationResult,
} from '../../src/ports/engine';

const request: SearchRequest = {
  program: 'blastn',
  query: { text: '>q\nACGT\n' },
  subject: { text: '>s\nACGT\n' },
  parameters: [],
  requestedThreads: 'auto',
};

function setup(engine: EngineGateway = new FakeEngine()) {
  const data = new MemoryDataGateway();
  const saved: Array<{ fileName: string; text: string }> = [];
  const downloader: Downloader = {
    save: (fileName, bytes) => saved.push({ fileName, text: new TextDecoder().decode(bytes) }),
  };
  let id = 0;
  const coordinator = new Coordinator({
    engine,
    data,
    downloader,
    digest: async (bytes) => `sha-${bytes.length}`,
    now: () => 1000,
    newRunId: () => `run-${++id}`,
  });
  return { coordinator, data, saved };
}

const statusOf = (coordinator: Coordinator, runId: string): RunStatus | undefined =>
  coordinator.state.get().runs.find((run) => run.snapshot.runId === runId)?.status;

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
  run(req: EngineRunRequest, sink: RunSink, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    this.started.push(req.runId);
    onPhase('running');
    return new Promise((resolve, reject) => {
      this.pending.set(req.runId, {
        resolve: () => {
          sink.write(6, new TextEncoder().encode(`${req.runId}\n`));
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
  it('runs one job at a time in FIFO order and stores all three outputs', async () => {
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
  });

  it('freezes the snapshot so later edits cannot change a queued job', async () => {
    const engine = new ManualEngine();
    const { coordinator } = setup(engine);
    const params: Array<[string, string]> = [['-evalue', '10']];
    const runId = (await coordinator.enqueue({ ...request, parameters: params })).runId!;
    params[0]![1] = '1e-50';
    const snapshot = coordinator.state.get().runs.find((r) => r.snapshot.runId === runId)!.snapshot;
    expect(Object.isFrozen(snapshot)).toBe(true);
    expect(snapshot.argv).toContain('10');
    expect(snapshot.argv).not.toContain('1e-50');
    expect(snapshot.query.name).toBe('query.fa');
    expect(snapshot.subject.sha256).toBe('sha-8');
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
    const { coordinator, data } = setup(engine);
    const runId = (await coordinator.enqueue(request)).runId!;
    await waitFor(coordinator, runId, 'running');
    coordinator.cancel(runId);
    engine.lateFinish?.();
    expect(statusOf(coordinator, runId)).toBe('running');
    engine.finish(runId);
    await waitFor(coordinator, runId, 'cancelled');
    await expect(data.readOutput(runId, 6)).rejects.toThrow();
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
    const failed = coordinator.state.get().runs.find((r) => r.snapshot.runId === second)!;
    expect(failed.record.error).toBe('engine error');
    expect(new TextDecoder().decode(await data.readOutput(first, 6))).toBe(`${first}\n`);
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
});
