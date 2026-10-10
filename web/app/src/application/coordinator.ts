// Application coordinator (plan §5.5): one active run, a FIFO queue, and the run state
// machine  queued -> preparing -> running -> finalizing -> completed | cancelled | failed.
//
// Ordering rule for cancel versus commit: a run can be cancelled until it enters
// `finalizing`; after that, cancel requests are ignored and the run completes.
//
// Inputs go through the data layer (plan §5.4): each input becomes a source, the source
// gets a record table (a dataset revision), and the run snapshot holds the bytes of the
// included records with their IDs and lengths, which the engine checks at `register`.
import { buildArgv, PASTED_NAMES } from '../domain/argv';
import type { FastaParserKind } from '../domain/dataset';
import type { OutputFormat } from '../domain/output-format';
import { indexParser, type InputRole, type ProgramId } from '../domain/programs';
import { isTerminal, type InputSnapshot, type RunRecord, type RunSnapshot, type RunStatus } from '../domain/run';
import type { SessionInput } from '../domain/session-file';
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
import { writeFile } from './export-writer';
import { Store } from './store';

/**
 * One input of a search: dataset revisions of the data layer (the search form's sources,
 * already indexed, plan §5.4) with the name passed as -query / -subject; or pasted text or
 * a file, which the coordinator indexes itself.
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
  readonly requestedThreads: number | 'auto';
  /** The run's name (NCBI's "Job Title"); the snapshot keeps it trimmed, or not at all when empty. */
  readonly title?: string;
}

export interface RunView {
  readonly snapshot: RunSnapshot;
  readonly status: RunStatus;
  readonly record: RunRecord;
  readonly result?: ResultSetRef;
  /** Set for a run loaded from a session file: where it comes from. Such a run is completed and never searched here. */
  readonly fromSession?: SessionOrigin;
  /**
   * The original FASTA re-attached to a run loaded from a session file, per role (REQ-23). It is
   * kept apart from the snapshot, which stays as the file recorded it.
   */
  readonly attached?: Readonly<Partial<Record<InputRole, Attachment>>>;
}

/** Where a run loaded from a session file comes from (design §12.2). */
export interface SessionOrigin {
  /** The session file's name, as the chosen file has it. */
  readonly fileName: string;
  /** The run's number in the working session that saved the file. */
  readonly number: number;
  /** When the file was saved (ms since the epoch). */
  readonly savedAt: number;
  /** The identity of the run's inputs that the file recorded: what an original FASTA must match to be attached. */
  readonly inputs: Readonly<Record<InputRole, SessionInput>>;
}

/** An original FASTA re-attached to a loaded run: the revisions whose run input matched the session file's record of it. */
export interface Attachment {
  readonly revisionIds: readonly string[];
  /** The chosen files' names, in order. */
  readonly fileNames: readonly string[];
}

/** A run loaded from a session file, before the coordinator numbers it (`addSessionRuns`). */
export interface SessionRunInit {
  readonly snapshot: Omit<RunSnapshot, 'number'>;
  readonly record: RunRecord;
  readonly result: ResultSetRef;
  readonly fromSession: SessionOrigin;
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
  /** Largest range of a stored output read at once by `exportOutput` (`EXPORT_RANGE_BYTES`); tests lower it. */
  readonly exportRangeBytes?: number;
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
/** The ranges in which `exportOutput` reads a stored output: the Data worker's read size. */
export const EXPORT_RANGE_BYTES = 8 * 1024 * 1024;

export class Coordinator {
  readonly state = new Store<AppState>({ runs: [] });
  private readonly queue: string[] = [];
  private readonly cancelRequested = new Set<string>();
  /** Run inputs by their revisions (snapshotInput). */
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
    const prepared: Array<{ request: SearchRequest; argv: readonly string[]; queryName: string; subjectName: string }> = [];
    for (const request of requests) {
      const queryName = inputName(request.query, PASTED_NAMES.query);
      const subjectName = inputName(request.subject, PASTED_NAMES.subject);
      let argv: readonly string[];
      try {
        argv = buildArgv({ program: request.program, queryName, subjectName, parameters: request.parameters });
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

    const snapshots: RunSnapshot[] = [];
    const groupId = prepared.length > 1 ? this.deps.newRunId() : undefined;
    for (const [index, { request, argv, queryName, subjectName }] of prepared.entries()) {
      let query: InputSnapshot;
      let subject: InputSnapshot;
      try {
        query = await this.snapshotInput('Query', request.query, queryName, indexParser(request.program, 'query'));
        subject = await this.snapshotInput('Subject', request.subject, subjectName, indexParser(request.program, 'subject'));
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

  /**
   * Adds runs loaded from a session file (application/session.ts) as completed runs, numbered
   * after the runs of this working session, in order. They are never queued, validated or
   * searched: their outputs are already in the Data worker.
   */
  addSessionRuns(runs: readonly SessionRunInit[]): readonly RunView[] {
    const views = runs.map(
      (run): RunView => ({
        snapshot: Object.freeze({ ...run.snapshot, number: this.nextNumber++ }),
        status: 'completed',
        record: run.record,
        result: run.result,
        fromSession: run.fromSession,
      }),
    );
    this.setRuns([...this.state.get().runs, ...views]);
    void this.refreshStorage();
    return views;
  }

  /** Records the original FASTA of a role of a loaded run, once the session found that it matches (REQ-23). */
  attach(runId: string, role: InputRole, attachment: Attachment): void {
    const view = this.find(runId);
    if (view?.fromSession === undefined) throw new Error('only a run loaded from a session file takes its original FASTA again');
    this.update(runId, { attached: { ...view.attached, [role]: attachment } });
  }

  /**
   * Saves one compatibility output of a completed run, the whole of it byte for byte as stored:
   * read in ranges of `exportRangeBytes` and written in order (design §12.1), so the output is
   * never held whole. A failed read saves nothing.
   */
  async exportOutput(runId: string, format: OutputFormat): Promise<void> {
    const view = this.find(runId);
    if (view?.status !== 'completed') throw new Error('only completed runs can be exported');
    if (view.result === undefined) throw new Error('the run has no stored result');
    const length = view.result.byteLengths[format];
    const range = this.deps.exportRangeBytes ?? EXPORT_RANGE_BYTES;
    if (!Number.isSafeInteger(range) || range < 1) throw new RangeError(`a range of ${range} bytes`);
    const { number, program } = view.snapshot;
    await writeFile(this.deps.downloader, `losat-run${number}-${program}.outfmt${format}.txt`, 'text/plain', async (writer) => {
      for (let start = 0; start < length; start += range) {
        await writer.bytes(await this.deps.data.readOutputRange(runId, format, start, Math.min(length, start + range)));
      }
    });
  }

  async readOutput(runId: string, format: OutputFormat): Promise<string> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    return new TextDecoder().decode(await this.deps.data.readOutput(runId, format));
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
   * and files first become a source with a record table, read with `parser`, the index
   * scan's reader for the program's sequence kind of the role. The run input of the same
   * revisions is built once and shared by the snapshots that use it (a subject kept for
   * many searches is held once).
   */
  private async snapshotInput(
    role: string,
    input: SequenceInput,
    name: string,
    parser: FastaParserKind,
  ): Promise<InputSnapshot> {
    try {
      let revisionIds: readonly string[];
      if ('dataset' in input) {
        revisionIds = Object.freeze([...input.dataset.revisionIds]);
      } else {
        const file = 'file' in input ? input.file : new File([input.text], name, { type: 'text/plain' });
        const source = await this.deps.data.addSource(file);
        const revision = await this.deps.data.indexSource(source.sourceId, parser);
        revisionIds = Object.freeze([revision.revisionId]);
      }
      const key = JSON.stringify(revisionIds);
      let runInput = this.runInputs.get(key);
      if (runInput === undefined) {
        runInput = this.deps.data.buildRunInput(revisionIds);
        this.runInputs.set(key, runInput);
        runInput.catch(() => this.runInputs.delete(key));
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
      const output = await this.deps.data.openRun(runId);
      staged = true;
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

function inputName(input: SequenceInput, pastedName: string): string {
  if ('dataset' in input) return input.dataset.name;
  return 'file' in input ? input.file.name : pastedName;
}

function engineInput(input: InputSnapshot): EngineInput {
  // Only queued runs are searched, and a run loaded from a session file is never queued.
  if (input.bytes === undefined) throw new Error('a run loaded from a session file is never searched again');
  return { bytes: input.bytes, sha256: input.sha256, revisionIds: input.revisionIds, records: input.records };
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
