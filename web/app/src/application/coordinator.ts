// Application coordinator (plan §5.5): one active run, a FIFO queue, and the run state
// machine  queued -> preparing -> running -> finalizing -> completed | cancelled | failed.
//
// Ordering rule for cancel versus commit: a run can be cancelled until it enters
// `finalizing`; after that, cancel requests are ignored and the run completes.
//
// Inputs go through the data layer (plan §5.4): each input becomes a source, the source
// gets a record table (a dataset revision), and the run snapshot holds the bytes of the
// included records with their IDs and lengths, which the engine checks at `register`.
import { buildArgv, optionValue, PASTED_NAMES } from '../domain/argv';
import type { FastaParserKind } from '../domain/dataset';
import { outputFileName, outputMediaType, runOutputs, type RunOutput } from '../domain/output-format';
import { fastaParser, programById, type ProgramId } from '../domain/programs';
import { isTerminal, type InputSnapshot, type RunRecord, type RunSnapshot, type RunStatus } from '../domain/run';
import type { DataGateway, ResultSetRef, RunInput, StorageInfo } from '../ports/data';
import type { Downloader } from '../ports/download';
import {
  RunCancelledError,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type HspRecord,
  type ValidationResult,
} from '../ports/engine';
import { taskOf } from './results';
import { Store } from './store';

/**
 * One input of a search: dataset revisions of the data layer (the search form's sources,
 * already indexed, plan §5.4) with the name passed as -query / -subject, whose sources the
 * caller keeps until enqueueAll returns; or pasted text or a file, which the coordinator
 * indexes itself. A run keeps its input in its snapshot (bytes and records), not its sources.
 */
export type SequenceInput =
  | { readonly dataset: { readonly name: string; readonly revisionIds: readonly string[] } }
  | { readonly text: string }
  | { readonly file: File };

export interface SearchRequest {
  readonly program: ProgramId;
  readonly query: SequenceInput;
  readonly subject: SequenceInput;
  readonly parameters: ReadonlyArray<readonly [string, string | true]>;
  /** Words of other options, after the parameters ("Additional options", ArgvInput.extra). */
  readonly extra?: readonly string[];
  readonly requestedThreads: number | 'auto';
  /** The run's name (NCBI's "Job Title"); the snapshot keeps it trimmed, or not at all when empty. */
  readonly title?: string;
}

export interface RunView {
  readonly snapshot: RunSnapshot;
  readonly status: RunStatus;
  readonly record: RunRecord;
  readonly result?: ResultSetRef;
}

export interface AppState {
  readonly runs: readonly RunView[];
  /** Temporary storage in use, shown without blocking any action (plan §5.6). */
  readonly storage?: StorageInfo;
}

export interface CoordinatorDeps {
  readonly engine: EngineGateway;
  readonly data: DataGateway;
  readonly downloader: Downloader;
  readonly now: () => number;
  readonly newRunId: () => string;
}

export type EnqueueResult = ValidationResult & { readonly runId?: string };
export type EnqueueAllResult = ValidationResult & { readonly runIds?: readonly string[] };

const PHASE_STATUS: Readonly<Record<EnginePhase, RunStatus>> = {
  preparing: 'preparing',
  running: 'running',
  finalizing: 'finalizing',
};
/** How often the storage status is read again while the start-up cleanup runs. */
const CLEANUP_POLL_MS = 250;

export class Coordinator {
  readonly state = new Store<AppState>({ runs: [] });
  private readonly queue: string[] = [];
  private readonly cancelRequested = new Set<string>();
  /**
   * Run inputs of dataset revisions, by the revisions (snapshotInput). An entry outlives the
   * sources of its revisions and shares its bytes with the snapshots that use it; an operation
   * that removes runs must remove the entries that no run uses.
   */
  private readonly runInputs = new Map<string, Promise<RunInput>>();
  private active: string | undefined;
  private nextNumber = 1;

  constructor(private readonly deps: CoordinatorDeps) {
    void this.refreshStorage();
  }

  /** Validates the request with the engine and, if valid, freezes it as a queued run. */
  async enqueue(request: SearchRequest): Promise<EnqueueResult> {
    const result = await this.enqueueAll([request]);
    return result.ok ? { ok: true, runId: result.runIds![0]! } : result;
  }

  /**
   * Validates every request and, if all are valid, freezes them as queued runs in order.
   * Two or more requests are one group of separate searches (plan §5.2): they share a
   * group ID, and `cancelGroup` cancels them together.
   */
  async enqueueAll(requests: readonly SearchRequest[]): Promise<EnqueueAllResult> {
    if (requests.length === 0) return { ok: false, message: 'There is nothing to search.' };
    const prepared: Prepared[] = [];
    for (const request of requests) {
      const queryName = inputName(request.query, PASTED_NAMES.query);
      const subjectName = inputName(request.subject, PASTED_NAMES.subject);
      let argv: readonly string[];
      try {
        argv = buildArgv({ program: request.program, queryName, subjectName, parameters: request.parameters, extra: request.extra ?? [] });
      } catch (error) {
        return { ok: false, message: errorMessage(error) };
      }
      let validation: ValidationResult;
      try {
        validation = await this.deps.engine.validate(argv);
      } catch (error) {
        return { ok: false, message: `The options could not be checked: ${errorMessage(error)}` };
      }
      if (!validation.ok) return validation;
      prepared.push({ request, argv, queryName, subjectName });
    }

    // The snapshots hold what the runs need: the sources made here from text and files are
    // dropped once they are read (SURVEY #15), and a request that fails caches nothing.
    const made: Made = { sources: [], runInputs: [] };
    let result: EnqueueAllResult | undefined;
    try {
      result = await this.freeze(prepared, made);
      return result;
    } finally {
      for (const sourceId of made.sources) void this.deps.data.dropSource(sourceId).catch(() => undefined);
      if (result?.ok !== true) for (const key of made.runInputs) this.runInputs.delete(key);
    }
  }

  private async freeze(prepared: readonly Prepared[], made: Made): Promise<EnqueueAllResult> {
    const snapshots: RunSnapshot[] = [];
    const groupId = prepared.length > 1 ? this.deps.newRunId() : undefined;
    for (const [index, { request, argv, queryName, subjectName }] of prepared.entries()) {
      const program = programById(request.program);
      let query: InputSnapshot;
      let subject: InputSnapshot;
      try {
        query = await this.snapshotInput('Query', request.query, queryName, fastaParser(program, 'query'), made);
        subject = await this.snapshotInput('Subject', request.subject, subjectName, fastaParser(program, 'subject'), made);
      } catch (error) {
        return { ok: false, message: errorMessage(error) };
      }
      const title = request.title?.trim();
      snapshots.push(
        Object.freeze({
          runId: this.deps.newRunId(),
          number: 0,
          program: request.program,
          ...(title ? { title } : {}),
          argv,
          query,
          subject,
          requestedThreads: request.requestedThreads,
          queuedAt: this.deps.now(),
          ...(groupId === undefined ? {} : { group: Object.freeze({ groupId, position: index + 1, size: prepared.length }) }),
        }),
      );
    }
    // Numbers are given only now, so that a request that fails leaves no gap.
    const numbered = snapshots.map((snapshot) => Object.freeze({ ...snapshot, number: this.nextNumber++ }));
    this.setRuns([
      ...this.state.get().runs,
      ...numbered.map((snapshot): RunView => ({ snapshot, status: 'queued', record: {} })),
    ]);
    for (const snapshot of numbered) this.queue.push(snapshot.runId);
    this.pump();
    return { ok: true, runIds: numbered.map((snapshot) => snapshot.runId) };
  }

  /**
   * Cancels the runs of a group: the running one, if it can still be cancelled, and the
   * queued ones (plan §5.2). Finished runs and a run that is finalizing are left alone;
   * the operation is not atomic.
   */
  cancelGroup(groupId: string): void {
    const members = this.state.get().runs.filter((run) => run.snapshot.group?.groupId === groupId);
    // Queued runs first, so that cancelling the active run does not start the next one of the group.
    for (const run of members) if (run.status === 'queued') this.cancel(run.snapshot.runId);
    for (const run of members) this.cancel(run.snapshot.runId);
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

  /** Saves one compatibility output of a completed run, byte for byte, named by outputFileName, with its format's media type. */
  async exportOutput(runId: string, output: RunOutput): Promise<void> {
    const view = this.find(runId);
    if (view?.status !== 'completed') throw new Error('only completed runs can be exported');
    const { number, program, argv, queuedAt } = view.snapshot;
    const [bytes, description] = await Promise.all([
      this.deps.data.readOutput(runId, output),
      // Without the engine's default task the name has none; the output is still saved.
      this.deps.engine.describe(program).catch(() => undefined),
    ]);
    const { task } = description === undefined ? { task: optionValue(argv, '-task') } : taskOf(argv, description);
    this.deps.downloader.save(outputFileName({ number, program, task, queuedAt, argv }, output), bytes, outputMediaType(argv, output));
  }

  async readOutput(runId: string, output: RunOutput): Promise<string> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    return new TextDecoder().decode(await this.deps.data.readOutput(runId, output));
  }

  /** The HSP records of a completed run (docs/web/abi_v2.md §8). */
  async readHits(runId: string): Promise<readonly HspRecord[]> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    return this.deps.data.readHits(runId);
  }

  /** The warnings of a completed run, as the CLI writes them to standard error. */
  async readDiagnostics(runId: string): Promise<string> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    return this.deps.data.readDiagnostics(runId);
  }

  /** Reads the storage status again; it keeps polling while the start-up cleanup runs. */
  async refreshStorage(): Promise<void> {
    let storage: StorageInfo;
    try {
      storage = await this.deps.data.storageInfo();
    } catch {
      return; // The status is informative; a failure to read it must not affect runs.
    }
    this.state.set({ ...this.state.get(), storage });
    if (storage.cleanup.state === 'pending') setTimeout(() => void this.refreshStorage(), CLEANUP_POLL_MS);
  }

  /**
   * Freezes the run input of an input. Dataset inputs use the revisions as they are; text
   * and files first become a source with a record table, which `made` lists to be dropped. The
   * run input of the same dataset revisions is built once and shared by the snapshots that
   * use it (a subject kept for many searches is held once); the revisions of text and files
   * are never used again.
   */
  private async snapshotInput(
    role: string,
    input: SequenceInput,
    name: string,
    parser: FastaParserKind,
    made: Made,
  ): Promise<InputSnapshot> {
    try {
      let runInput: Promise<RunInput>;
      let revisionIds: readonly string[];
      if ('dataset' in input) {
        revisionIds = Object.freeze([...input.dataset.revisionIds]);
        const key = JSON.stringify(revisionIds);
        const cached = this.runInputs.get(key);
        if (cached === undefined) {
          runInput = this.deps.data.buildRunInput(revisionIds);
          this.runInputs.set(key, runInput);
          made.runInputs.push(key);
          runInput.catch(() => this.runInputs.delete(key));
        } else {
          runInput = cached;
        }
      } else {
        const file = 'file' in input ? input.file : new File([input.text], name, { type: 'text/plain' });
        const source = await this.deps.data.addSource(file);
        made.sources.push(source.sourceId);
        const revision = await this.deps.data.indexSource(source.sourceId, parser);
        revisionIds = Object.freeze([revision.revisionId]);
        runInput = this.deps.data.buildRunInput(revisionIds);
      }
      const { bytes, sha256, records } = await runInput;
      return Object.freeze({ name, bytes, sha256, revisionIds, records });
    } catch (error) {
      throw new Error(`${role} FASTA: ${errorMessage(error)}`, { cause: error });
    }
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
    let staged = false;
    const startedAt = this.deps.now();
    try {
      this.update(runId, { status: 'preparing', record: { startedAt } });
      const output = await this.deps.data.openRun(runId, runOutputs(snapshot.argv));
      staged = true;
      // A cancel while the outputs were staged found no search to end: the engine has not started it.
      if (this.cancelRequested.has(runId)) {
        output.close();
        throw new RunCancelledError(runId);
      }
      const runtime = await this.deps.engine.run(
        {
          runId,
          argv: snapshot.argv,
          query: engineInput(snapshot.query),
          subject: engineInput(snapshot.subject),
          requestedThreads: snapshot.requestedThreads,
        },
        output,
        (phase) => this.onPhase(runId, phase),
      );
      if (this.cancelRequested.has(runId)) throw new RunCancelledError(runId);
      this.update(runId, { status: 'finalizing' });
      const result = await this.deps.data.commitRun(runId);
      this.update(runId, {
        status: 'completed',
        result,
        record: {
          endedAt: this.deps.now(),
          runtimePath: runtime.path,
          threads: runtime.threads,
          engineBuild: runtime.engineBuild,
          ...(runtime.fallbackReason === undefined ? {} : { fallbackReason: runtime.fallbackReason }),
          ...(runtime.runtimeGeneration === undefined ? {} : { runtimeGeneration: runtime.runtimeGeneration }),
          ...(runtime.memory === undefined ? {} : { memory: runtime.memory }),
          ...(runtime.subjectRetained === undefined ? {} : { subjectRetained: runtime.subjectRetained }),
        },
      });
    } catch (error) {
      if (staged) await this.deps.data.discardRun(runId).catch(() => undefined);
      const cancelled = error instanceof RunCancelledError || this.cancelRequested.has(runId);
      this.update(runId, {
        status: cancelled ? 'cancelled' : 'failed',
        record: {
          endedAt: this.deps.now(),
          ...(cancelled ? {} : { error: errorMessage(error) }),
        },
      });
    } finally {
      this.cancelRequested.delete(runId);
      this.active = undefined;
      void this.refreshStorage();
      this.pump();
    }
  }

  /** Phase events are accepted only for the active, non-cancelled run. */
  private onPhase(runId: string, phase: EnginePhase): void {
    if (this.active !== runId || this.cancelRequested.has(runId)) return;
    const view = this.find(runId);
    if (view === undefined || isTerminal(view.status) || view.status === 'finalizing') return;
    const phaseTimes = { ...view.record.phaseTimes, [phase]: this.deps.now() };
    this.update(runId, { status: PHASE_STATUS[phase], record: { phaseTimes } });
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
    this.state.set({ ...this.state.get(), runs: Object.freeze([...runs]) });
  }
}

/** A request that the engine validated, with its argv and the names of its inputs. */
interface Prepared {
  readonly request: SearchRequest;
  readonly argv: readonly string[];
  readonly queryName: string;
  readonly subjectName: string;
}

/** What one enqueueAll made: sources of text and files, and the cached run inputs of dataset revisions. */
interface Made {
  readonly sources: string[];
  readonly runInputs: string[];
}

function inputName(input: SequenceInput, pastedName: string): string {
  if ('dataset' in input) return input.dataset.name;
  return 'file' in input ? input.file.name : pastedName;
}

function engineInput(input: InputSnapshot): EngineInput {
  return { bytes: input.bytes, sha256: input.sha256, revisionIds: input.revisionIds, records: input.records };
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
