// Application coordinator (plan §5.5): one active run, a FIFO queue, and the run state
// machine  queued -> preparing -> running -> finalizing -> completed | cancelled | failed.
//
// Ordering rule for cancel versus commit: a run can be cancelled until it enters
// `finalizing`; after that, cancel requests are ignored and the run completes.
import { buildArgv } from '../domain/argv';
import type { OutputFormat } from '../domain/output-format';
import type { ProgramId } from '../domain/programs';
import { isTerminal, type RunRecord, type RunSnapshot, type RunStatus } from '../domain/run';
import type { DataGateway, ResultSetRef, RunStaging } from '../ports/data';
import type { Downloader } from '../ports/download';
import {
  RunCancelledError,
  type EngineGateway,
  type EnginePhase,
  type ValidationResult,
} from '../ports/engine';
import { Store } from './store';

export interface SequenceInput {
  /** File name, or undefined for pasted text. */
  readonly fileName?: string;
  readonly text: string;
}

export interface SearchRequest {
  readonly program: ProgramId;
  readonly query: SequenceInput;
  readonly subject: SequenceInput;
  readonly parameters: ReadonlyArray<readonly [string, string | true]>;
  readonly requestedThreads: number | 'auto';
}

export interface RunView {
  readonly snapshot: RunSnapshot;
  readonly status: RunStatus;
  readonly record: RunRecord;
  readonly result?: ResultSetRef;
}

export interface AppState {
  readonly runs: readonly RunView[];
}

export interface CoordinatorDeps {
  readonly engine: EngineGateway;
  readonly data: DataGateway;
  readonly downloader: Downloader;
  readonly digest: (bytes: Uint8Array) => Promise<string>;
  readonly now: () => number;
  readonly newRunId: () => string;
}

export type EnqueueResult = ValidationResult & { readonly runId?: string };

const PASTED_NAMES = { query: 'query.fa', subject: 'subject.fa' } as const;
const PHASE_STATUS: Readonly<Record<EnginePhase, RunStatus>> = {
  preparing: 'preparing',
  running: 'running',
  finalizing: 'finalizing',
};

export class Coordinator {
  readonly state = new Store<AppState>({ runs: [] });
  private readonly queue: string[] = [];
  private readonly cancelRequested = new Set<string>();
  private active: string | undefined;
  private nextNumber = 1;

  constructor(private readonly deps: CoordinatorDeps) {}

  /** Validates the request with the engine and, if valid, freezes it as a queued run. */
  async enqueue(request: SearchRequest): Promise<EnqueueResult> {
    const queryName = request.query.fileName ?? PASTED_NAMES.query;
    const subjectName = request.subject.fileName ?? PASTED_NAMES.subject;
    let argv: readonly string[];
    try {
      argv = buildArgv({ program: request.program, queryName, subjectName, parameters: request.parameters });
    } catch (error) {
      return { ok: false, message: errorMessage(error) };
    }
    const validation = await this.deps.engine.validate(argv);
    if (!validation.ok) return validation;

    const encoder = new TextEncoder();
    const queryBytes = encoder.encode(request.query.text);
    const subjectBytes = encoder.encode(request.subject.text);
    const runId = this.deps.newRunId();
    const snapshot: RunSnapshot = Object.freeze({
      runId,
      number: this.nextNumber++,
      program: request.program,
      argv,
      query: Object.freeze({ name: queryName, bytes: queryBytes, sha256: await this.deps.digest(queryBytes) }),
      subject: Object.freeze({
        name: subjectName,
        bytes: subjectBytes,
        sha256: await this.deps.digest(subjectBytes),
      }),
      requestedThreads: request.requestedThreads,
      queuedAt: this.deps.now(),
    });
    this.setRuns([...this.state.get().runs, { snapshot, status: 'queued', record: {} }]);
    this.queue.push(runId);
    this.pump();
    return { ok: true, runId };
  }

  /** Cancels a queued or running run. Runs that are finalizing or finished are left alone. */
  cancel(runId: string): void {
    const view = this.find(runId);
    if (view === undefined || isTerminal(view.status) || view.status === 'finalizing') return;
    const queuedIndex = this.queue.indexOf(runId);
    if (queuedIndex >= 0) {
      this.queue.splice(queuedIndex, 1);
      this.update(runId, { status: 'cancelled', record: { ...view.record, endedAt: this.deps.now() } });
      return;
    }
    this.cancelRequested.add(runId);
    this.deps.engine.cancel(runId);
  }

  /** Saves one compatibility output of a completed run, byte for byte. */
  async exportOutput(runId: string, format: OutputFormat): Promise<void> {
    const view = this.find(runId);
    if (view?.status !== 'completed') throw new Error('only completed runs can be exported');
    const bytes = await this.deps.data.readOutput(runId, format);
    const { number, program } = view.snapshot;
    this.deps.downloader.save(`losat-run${number}-${program}.outfmt${format}.txt`, bytes, 'text/plain');
  }

  async readOutput(runId: string, format: OutputFormat): Promise<string> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    return new TextDecoder().decode(await this.deps.data.readOutput(runId, format));
  }

  private pump(): void {
    if (this.active !== undefined) return;
    const next = this.queue.shift();
    if (next === undefined) return;
    this.active = next;
    void this.execute(next);
  }

  private async execute(runId: string): Promise<void> {
    const view = this.find(runId);
    if (view === undefined) throw new Error(`missing run ${runId}`);
    const { snapshot } = view;
    let staging: RunStaging | undefined;
    const startedAt = this.deps.now();
    try {
      this.update(runId, { status: 'preparing', record: { startedAt } });
      staging = await this.deps.data.openRun(runId);
      const runtime = await this.deps.engine.run(
        {
          runId,
          argv: snapshot.argv,
          query: snapshot.query.bytes,
          subject: snapshot.subject.bytes,
          requestedThreads: snapshot.requestedThreads,
        },
        staging,
        (phase) => this.onPhase(runId, phase),
      );
      if (this.cancelRequested.has(runId)) throw new RunCancelledError(runId);
      this.update(runId, { status: 'finalizing' });
      const result = await staging.commit();
      this.update(runId, {
        status: 'completed',
        result,
        record: {
          startedAt,
          endedAt: this.deps.now(),
          runtimePath: runtime.path,
          threads: runtime.threads,
          engineBuild: runtime.engineBuild,
          ...(runtime.fallbackReason === undefined ? {} : { fallbackReason: runtime.fallbackReason }),
        },
      });
    } catch (error) {
      await staging?.discard();
      const cancelled = error instanceof RunCancelledError || this.cancelRequested.has(runId);
      this.update(runId, {
        status: cancelled ? 'cancelled' : 'failed',
        record: {
          startedAt,
          endedAt: this.deps.now(),
          ...(cancelled ? {} : { error: errorMessage(error) }),
        },
      });
    } finally {
      this.cancelRequested.delete(runId);
      this.active = undefined;
      this.pump();
    }
  }

  /** Phase events are accepted only for the active, non-cancelled run. */
  private onPhase(runId: string, phase: EnginePhase): void {
    if (this.active !== runId || this.cancelRequested.has(runId)) return;
    const view = this.find(runId);
    if (view === undefined || isTerminal(view.status) || view.status === 'finalizing') return;
    this.update(runId, { status: PHASE_STATUS[phase] });
  }

  private find(runId: string): RunView | undefined {
    return this.state.get().runs.find((run) => run.snapshot.runId === runId);
  }

  private update(runId: string, change: Partial<Omit<RunView, 'snapshot'>>): void {
    this.setRuns(
      this.state.get().runs.map((run) =>
        run.snapshot.runId === runId ? { ...run, ...change, record: { ...run.record, ...change.record } } : run,
      ),
    );
  }

  private setRuns(runs: readonly RunView[]): void {
    this.state.set({ runs: Object.freeze([...runs]) });
  }
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
