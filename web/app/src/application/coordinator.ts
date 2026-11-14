// Application coordinator (plan §5.5): one active run, a FIFO queue, and the run state
// machine  queued -> preparing -> running -> finalizing -> completed | cancelled | failed.
//
// Ordering rule for cancel versus commit: a run can be cancelled until it enters
// `finalizing`; after that, cancel requests are ignored and the run completes.
//
// Inputs go through the data layer (plan §5.4): each input becomes a source, the source
// gets a record table (a dataset revision), and the run snapshot holds the bytes of the
// included records with their IDs and lengths, which the engine checks at `register`.
import { WINDOW_MASKER_FLAG, WINDOW_MASKER_TAXID_FLAG, type FastaParserKind, type ProgramId, type RunFile, type RunOutput } from 'losat/domain';
import {
  RunCancelledError,
  runFilesOf,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type HspRecord,
  type StartResult,
  type ValidationResult,
} from 'losat/ports';
import { buildArgv, optionValue, PASTED_NAMES } from '../domain/argv';
import { countsName, parseTaxid, type CountsOrigin } from '../domain/counts';
import { outputFileName, outputMediaType, runFileName, runOutputs } from '../domain/output-format';
import { fastaParser, programById } from '../domain/programs';
import { isTerminal, type ConffileSnapshot, type InputSnapshot, type RunRecord, type RunSnapshot, type RunStatus } from '../domain/run';
import type { DataGateway, ResultSetRef, RunInput, StorageInfo } from '../ports/data';
import type { Downloader } from '../ports/download';
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
  /**
   * The counts file of the search (BLASTN): the file that the parameters' `-window_masker_db` names,
   * which the run reads under that name, or the file kept for the `-window_masker_taxid` (line WW2),
   * which it reads under the name that the engine looks for (domain/counts.ts).
   */
  readonly windowMasker?: File;
  /** Where the counts of a taxonomy ID came from; the run records it. */
  readonly windowMaskerOrigin?: CountsOrigin;
  /**
   * The file chosen for the argv's `-conffile` (Additional options; line RG E1w-3b): read once when
   * the runs are queued (once for the requests of a group that share it), the run reads its bytes under
   * the name that the engine's start reads of the argv. Not used by an argv without `-conffile`.
   */
  readonly conffile?: File;
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
  /** The lower-case hex SHA-256 of bytes (the package's `sha256Hex`), for the registry file of a snapshot. */
  readonly digest: (bytes: Uint8Array) => Promise<string>;
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
  /**
   * The bytes of the registry file of each queued or running run whose argv has `-conffile` (read at
   * enqueue, shared by a group), by run: the run's input until it ends, when they go; its snapshot
   * keeps the name, the size and the SHA-256.
   */
  private readonly conffiles = new Map<string, Blob>();
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
    // Each registry file read once, with its SHA-256, for the requests of a group that share it.
    const reads = new Map<File, Promise<{ readonly file: Blob; readonly sha256: string }>>();
    for (const request of requests) {
      const queryName = inputName(request.query, PASTED_NAMES.query);
      const subjectName = inputName(request.subject, PASTED_NAMES.subject);
      let argv: readonly string[];
      try {
        argv = buildArgv({ program: request.program, queryName, subjectName, parameters: request.parameters, extra: request.extra ?? [] });
      } catch (error) {
        return { ok: false, message: errorMessage(error) };
      }
      // NCBI's application start first: the name of the registry file that the argv reads, and
      // whether the run is a dry run (line RG E1w-3b).
      let start: StartResult;
      try {
        start = await this.deps.engine.start(argv);
      } catch (error) {
        return { ok: false, message: `The options could not be checked: ${errorMessage(error)}` };
      }
      if (!start.ok) return start;
      let conffile: Prepared['conffile'];
      const { conffile: conffileName, dryRun } = start.start;
      // An empty name (`-conffile=`) names no file that the engine can hold: its message stands.
      if (conffileName !== null && conffileName !== '' && request.conffile !== undefined) {
        const chosen = request.conffile;
        let read = reads.get(chosen);
        if (read === undefined) reads.set(chosen, (read = this.readConffile(chosen)));
        try {
          const { file, sha256 } = await read;
          conffile = { snapshot: Object.freeze({ name: conffileName, fileName: chosen.name, size: file.size, sha256 }), file };
        } catch (error) {
          return { ok: false, message: `The file of -conffile ${conffileName} (${chosen.name}) could not be read: ${errorMessage(error)}` };
        }
      }
      let validation: ValidationResult;
      try {
        validation = await this.deps.engine.validate(argv, conffile === undefined ? undefined : { name: conffile.snapshot.name, file: conffile.file });
      } catch (error) {
        return { ok: false, message: `The options could not be checked: ${errorMessage(error)}` };
      }
      if (!validation.ok) return validation;
      prepared.push({ request, argv, queryName, subjectName, conffile, dryRun });
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
    /** The bytes of the registry file of each snapshot, in their order. */
    const conffiles: Array<Blob | undefined> = [];
    const groupId = prepared.length > 1 ? this.deps.newRunId() : undefined;
    for (const [index, { request, argv, queryName, subjectName, conffile, dryRun }] of prepared.entries()) {
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
      // The run reads the file under the name that its argv gives the option.
      const taxid = parseTaxid(optionValue(argv, WINDOW_MASKER_TAXID_FLAG) ?? '');
      const name = optionValue(argv, WINDOW_MASKER_FLAG) ?? (taxid === undefined ? undefined : countsName(taxid));
      const windowMasker =
        request.windowMasker === undefined || name === undefined
          ? undefined
          : Object.freeze({
              name,
              file: request.windowMasker,
              size: request.windowMasker.size,
              lastModified: request.windowMasker.lastModified,
              ...(taxid === undefined || request.windowMaskerOrigin?.taxid !== taxid ? {} : { origin: request.windowMaskerOrigin }),
            });
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
          ...(windowMasker === undefined ? {} : { windowMasker }),
          ...(conffile === undefined ? {} : { conffile: conffile.snapshot }),
          ...(dryRun ? { dryRun: true as const } : {}),
          queuedAt: this.deps.now(),
          ...(groupId === undefined ? {} : { group: Object.freeze({ groupId, position: index + 1, size: prepared.length }) }),
        }),
      );
      conffiles.push(conffile?.file);
    }
    // Numbers are given only now, so that a request that fails leaves no gap.
    const numbered = snapshots.map((snapshot) => Object.freeze({ ...snapshot, number: this.nextNumber++ }));
    this.setRuns([
      ...this.state.get().runs,
      ...numbered.map((snapshot): RunView => ({ snapshot, status: 'queued', record: {} })),
    ]);
    for (const [index, snapshot] of numbered.entries()) {
      const file = conffiles[index];
      if (file !== undefined) this.conffiles.set(snapshot.runId, file);
      this.queue.push(snapshot.runId);
    }
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
      this.conffiles.delete(runId);
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
    const [bytes, task] = await Promise.all([this.deps.data.readOutput(runId, output), this.taskOf(view.snapshot)]);
    this.deps.downloader.save(outputFileName({ number, program, task, queuedAt, argv }, output), bytes, outputMediaType(argv, output));
  }

  /**
   * Saves a file that a run made besides its outputs (RunRecord.runFiles: the file of `-logfile`), its
   * bytes as the engine wrote them, named by runFileName. A failed run keeps its files too.
   */
  async exportRunFile(runId: string, index: number): Promise<void> {
    const view = this.find(runId);
    const file: RunFile | undefined = view?.record.runFiles?.[index];
    if (view === undefined || file === undefined) throw new Error('the run has no such file');
    const { number, program, argv, queuedAt } = view.snapshot;
    const task = await this.taskOf(view.snapshot);
    this.deps.downloader.save(runFileName({ number, program, task, queuedAt, argv }, file.name), file.bytes, 'text/plain');
  }

  /**
   * Reads the file chosen for `-conffile` once: its bytes as a Blob (a copy that the file's later
   * changes do not reach), the run's input, and their SHA-256.
   */
  private async readConffile(chosen: File): Promise<{ readonly file: Blob; readonly sha256: string }> {
    const bytes = new Uint8Array(await chosen.arrayBuffer());
    return { file: new Blob([bytes]), sha256: await this.deps.digest(bytes) };
  }

  /** The task of a run for its files' names: the argv's, or the engine's default; none without the engine's description. */
  private async taskOf(snapshot: RunSnapshot): Promise<string | undefined> {
    // Without the engine's default task the name has none; the file is still saved.
    const description = await this.deps.engine.describe(snapshot.program).catch(() => undefined);
    return description === undefined ? optionValue(snapshot.argv, '-task') : taskOf(snapshot.argv, description).task;
  }

  /** Opens one output of a completed run without reading it (a binary one, which is not shown): rejects as `readOutput` when it cannot be read. */
  async checkOutput(runId: string, output: RunOutput): Promise<void> {
    if (this.find(runId)?.status !== 'completed') throw new Error('run is not completed');
    await this.deps.data.readOutputRange(runId, output, 0, 0);
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
    const conffile = this.conffiles.get(runId);
    let staged = false;
    const startedAt = this.deps.now();
    // The files that the run made besides its outputs, which a completed and a failed run keep.
    let runFiles: readonly RunFile[] | undefined;
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
          ...(snapshot.windowMasker === undefined ? {} : { windowMasker: snapshot.windowMasker }),
          ...(snapshot.conffile === undefined || conffile === undefined ? {} : { conffile: { name: snapshot.conffile.name, file: conffile } }),
        },
        output,
        (phase) => this.onPhase(runId, phase),
      ).catch((error: unknown) => {
        runFiles = runFilesOf(error);
        throw error;
      });
      runFiles = runtime.runFiles;
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
          ...(runtime.windowMasker === undefined
            ? {}
            : { windowMasker: { ...runtime.windowMasker, ...(snapshot.windowMasker?.origin === undefined ? {} : { origin: snapshot.windowMasker.origin }) } }),
          ...(runFiles === undefined ? {} : { runFiles }),
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
          ...(cancelled || runFiles === undefined ? {} : { runFiles }),
        },
      });
    } finally {
      // The run's registry file goes with its run; the snapshot keeps its name, size and SHA-256.
      this.conffiles.delete(runId);
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

/** A request that the engine validated, with its argv, the names of its inputs, and what the engine's start read of it. */
interface Prepared {
  readonly request: SearchRequest;
  readonly argv: readonly string[];
  readonly queryName: string;
  readonly subjectName: string;
  /** The registry file of the argv's `-conffile`, read from the chosen file: what the snapshot keeps of it, and its bytes. */
  readonly conffile: { readonly snapshot: ConffileSnapshot; readonly file: Blob } | undefined;
  readonly dryRun: boolean;
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
