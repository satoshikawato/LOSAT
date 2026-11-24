// Application coordinator (plan §5.5): one active run, a FIFO queue, and the run state
// machine  queued -> preparing -> running -> finalizing -> completed | cancelled | failed.
//
// Ordering rule for cancel versus commit: a run can be cancelled until it ends. A cancel while it
// is `finalizing` ends the wait for its results, which are deleted if they come later: a Data
// worker that stopped answering kept such a run "Organizing results" for good, with no Cancel
// (line WQA's sweep, D2, 2026-11-19). When the Data worker ends, the active run ends at once,
// with its search: it could not keep its results (the review of D2, L2).
//
// Inputs go through the data layer (plan §5.4): each input becomes a source, the source
// gets a record table (a dataset revision), and the run snapshot holds the bytes of the
// included records with their IDs and lengths, which the engine checks at `register`.
import {
  WINDOW_MASKER_FLAG,
  WINDOW_MASKER_TAXID_FLAG,
  type DirectoryFile,
  type FastaReader,
  type PrintedKind,
  type ProgramId,
  runFileOption,
  type RunFile,
  type RunOutput,
} from 'losat/domain';
import {
  isUnavailable,
  RunCancelledError,
  runFilesOf,
  runOutputFailedOf,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type HspRecord,
  type ProgramDescription,
  type StartResult,
  type ToolEngine,
  type ToolValidation,
  type ValidationResult,
} from 'losat/ports';
import { buildArgv, optionValue, runTask } from '../domain/argv';
import { countsName, parseTaxid, type CountsOrigin } from '../domain/counts';
import { inputFileName, outputFileName, outputMediaType, printedFileName, printedFileType, runFileName, runOutputs } from '../domain/output-format';
import {
  isTerminal,
  isToolSnapshot,
  type ConffileSnapshot,
  type InputPart,
  type InputSnapshot,
  type RunRecord,
  type RunSnapshot,
  type RunStatus,
  type ToolFileOrigin,
  type ToolOutputEntry,
  type ToolRunRecord,
  type ToolRunSnapshot,
} from '../domain/run';
import { isCountsOutputFile, toolFileOptionLabel } from '../domain/tool-argv';
import { toolOutputFileName, type ToolId, type WindowMaskerStage } from '../domain/tools';
import type { CountsLocks } from '../ports/counts';
import type { DataGateway, KeptRun, ResultSetRef, RunInput, StorageInfo } from '../ports/data';
import type { Downloader } from '../ports/download';
import { failureText } from './availability';
import { registryName, unreadableConffile } from './conffile';
import type { SessionFiles } from './session-files';
import { Store } from './store';
import { ToolStorageSink } from './tool-storage-sink';

/**
 * One input of a search: dataset revisions of the data layer (the search form's sources,
 * already indexed, plan §5.4) with the name passed as -query / -subject, whose sources the
 * caller keeps until enqueueAll returns; or pasted text or a file, which the coordinator
 * indexes itself. A run keeps its input in its snapshot (bytes and records), not its sources.
 */
export type SequenceInput =
  | {
      readonly dataset: {
        /** The user's file name, for the user's file as given (no `from`); a made input has none: the coordinator names it. */
        readonly name?: string;
        readonly revisionIds: readonly string[];
        /** What a file that LOSAT makes (pasted text, several files, a file with records left out) is made of. */
        readonly from?: readonly InputPart[];
      };
    }
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

export interface ToolDirectoryFile {
  readonly name: string;
  readonly file: Blob;
  readonly origin: ToolFileOrigin;
}

export interface ToolRequest {
  readonly tool: ToolId;
  readonly stage: WindowMaskerStage;
  readonly argv: readonly string[];
  readonly directory: readonly ToolDirectoryFile[];
  readonly prints?: PrintedKind;
  /** OPFS-backed inflated session file tokens read by this run (directly or inside combined_input.fa). */
  readonly sessionTokens?: readonly string[];
}

export interface SearchRunView {
  readonly snapshot: RunSnapshot;
  readonly status: RunStatus;
  readonly record: RunRecord;
  readonly result?: ResultSetRef;
}

export interface ToolRunView {
  readonly snapshot: ToolRunSnapshot;
  readonly status: RunStatus;
  readonly record: ToolRunRecord;
  readonly result?: undefined;
}

export type AnyRunView = SearchRunView | ToolRunView;
export type RunView = SearchRunView;

export function isToolRun(run: AnyRunView): run is ToolRunView {
  return isToolSnapshot(run.snapshot);
}

/** A counts file produced in this session by WindowMasker (`-mk_counts` or `-convert`), offered to Mask and Convert. */
export interface SessionCountsEntry {
  readonly runId: string;
  readonly runNumber: number;
  readonly key: string;
  readonly name: string;
  readonly size: number;
}

export interface AppState {
  readonly runs: readonly AnyRunView[];
  /** Temporary storage in use, shown without blocking any action (plan §5.6). */
  readonly storage?: StorageInfo;
}

export interface CoordinatorDeps {
  readonly engine: EngineGateway & Partial<ToolEngine>;
  readonly data: DataGateway;
  readonly downloader: Downloader;
  readonly now: () => number;
  readonly newRunId: () => string;
  /** The lower-case hex SHA-256 of bytes (the package's `sha256Hex`), for the registry file of a snapshot. */
  readonly digest: (bytes: Uint8Array) => Promise<string>;
  /** Settles when the Data worker could not start or ended, with the error of its calls from then on. */
  readonly dataWorkerEnded?: Promise<Error>;
  /** Runtime path and build identifier recorded on completed/failed tool runs. */
  readonly toolRuntime?: { readonly path: 'serial' | 'fake'; readonly engineBuild: string };
  /** Shared holder tracker for inflated session files used by WindowMasker tool runs. */
  readonly sessionFiles?: SessionFiles;
  /**
   * The taxonomy IDs' locks on saved counts: a run that reads a saved copy holds its ID's lock in shared mode from
   * enqueue to its end, so that Replace, Delete and a download of that ID refuse while it is queued or running.
   * Without it (a browser without Web Locks), runs hold nothing, as before.
   */
  readonly countsLocks?: CountsLocks;
}

/** Whether a request was queued, or why not (a refusal of its argv is the validation's). */
type Queued = { readonly ok: true } | { readonly ok: false; readonly message: string };
export type EnqueueResult = Queued & { readonly runId?: string };
export type EnqueueAllResult = Queued & { readonly runIds?: readonly string[] };

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
  /** Working directory files held for queued/running tool runs; cleared when the run finishes or is cancelled. */
  private readonly toolDirectories = new Map<string, readonly DirectoryFile[]>();
  /** The releases of the saved-counts locks that each queued or running run holds (holdCounts), by run. */
  private readonly countsHolds = new Map<string, () => void>();
  private active: string | undefined;
  /**
   * Ends the active run's wait for its search or for the commit of its results with an error: a cancel while its results
   * are organized (cancel), or the end of the Data worker (dataWorkerEnded).
   */
  private stopActive: ((error: Error) => void) | undefined;
  private toolAbort: AbortController | undefined;
  private nextNumber = 1;
  /** The enqueue before this one, which this one waits for (enqueueAll). */
  private enqueueTurn: Promise<void> = Promise.resolve();

  constructor(private readonly deps: CoordinatorDeps) {
    void this.refreshStorage();
    void deps.dataWorkerEnded?.then((error) => this.dataWorkerEnded(error));
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
    // One at a time: a made input is named by its run's number and time, which are known before the argv is built (`inputName`) and
    // must stay the run's own until its snapshot is made.
    const turn = this.enqueueTurn.then(() => this.enqueueInTurn(requests));
    this.enqueueTurn = turn.then(() => undefined, () => undefined);
    return turn;
  }

  /** Validates and queues one tool run in the shared queue. */
  async enqueueTool(request: ToolRequest): Promise<EnqueueResult> {
    const result = await this.enqueueToolAll([request]);
    return result.ok ? { ok: true, runId: result.runIds![0]! } : result;
  }

  /** Validates every tool request and queues them together in the single queue shared with searches. */
  async enqueueToolAll(requests: readonly ToolRequest[]): Promise<EnqueueAllResult> {
    if (requests.length === 0) return { ok: false, message: 'There is nothing to run.' };
    const turn = this.enqueueTurn.then(() => this.enqueueToolInTurn(requests));
    this.enqueueTurn = turn.then(() => undefined, () => undefined);
    return turn;
  }

  /** Whether a search or tool run is currently executing on the Engine worker. */
  isBusy(): boolean {
    return this.active !== undefined;
  }

  private async enqueueToolInTurn(requests: readonly ToolRequest[]): Promise<EnqueueAllResult> {
    if (this.deps.engine.validateTool === undefined || this.deps.engine.runTool === undefined) {
      return { ok: false, message: 'This engine does not support running tools.' };
    }
    const directories: Array<readonly DirectoryFile[]> = [];
    // `WasmEngine.validateTool` executes on the Engine worker behind any active run (`inTurn`).
    // Validate eagerly when the engine is idle; when another run is active, queue immediately so
    // the new run enters the shared queue (`Waiting for run N`) instead of blocking on the worker.
    const validateNow = this.active === undefined;
    for (const request of requests) {
      const directory: DirectoryFile[] = request.directory.map((entry) => ({ name: entry.name, file: entry.file }));
      if (validateNow) {
        let verdict: ToolValidation;
        try {
          verdict = await this.deps.engine.validateTool(request.argv, { directory });
        } catch (error) {
          return { ok: false, message: uncheckedOptions(error) };
        }
        if (!verdict.ok) return verdict;
      }
      directories.push(directory);
    }
    const queuedAt = this.deps.now();
    const groupId = requests.length > 1 ? this.deps.newRunId() : undefined;
    const snapshots: ToolRunSnapshot[] = requests.map((request, index) =>
      Object.freeze({
        kind: 'tool' as const,
        runId: this.deps.newRunId(),
        number: this.nextNumber + index,
        tool: request.tool,
        stage: request.stage,
        argv: Object.freeze([...request.argv]),
        directory: Object.freeze(
          request.directory.map((entry) =>
            Object.freeze({ name: entry.name, size: entry.file.size, origin: entry.origin }),
          ),
        ),
        ...(request.prints === undefined ? {} : { prints: request.prints }),
        queuedAt,
        ...(groupId === undefined ? {} : { group: Object.freeze({ groupId, position: index + 1, size: requests.length }) }),
      }),
    );
    const refused = await this.holdCounts(
      snapshots.flatMap((snapshot, index) =>
        requests[index]!.directory.flatMap((entry) =>
          entry.origin.kind === 'saved' ? [{ runId: snapshot.runId, taxid: entry.origin.taxid, file: entry.file }] : [],
        ),
      ),
    );
    if (refused !== undefined) return { ok: false, message: refused };
    this.nextNumber += snapshots.length;
    this.setRuns([
      ...this.state.get().runs,
      ...snapshots.map((snapshot): ToolRunView => ({ snapshot, status: 'queued', record: {} })),
    ]);
    for (const [index, snapshot] of snapshots.entries()) {
      const tokens = requests[index]!.sessionTokens;
      if (tokens !== undefined && tokens.length > 0) this.deps.sessionFiles?.set(`run:${snapshot.runId}`, tokens);
      this.toolDirectories.set(snapshot.runId, directories[index]!);
      this.queue.push(snapshot.runId);
    }
    this.pump();
    return { ok: true, runIds: snapshots.map((snapshot) => snapshot.runId) };
  }

  private async enqueueInTurn(requests: readonly SearchRequest[]): Promise<EnqueueAllResult> {
    const prepared: Prepared[] = [];
    const queuedAt = this.deps.now();
    // The engine's description of a program, asked once for the requests of this press (a group has one program).
    const descriptions = new Map<ProgramId, Promise<ProgramDescription | undefined>>();
    // Each registry file read once, with its SHA-256, for the requests of a group that share it.
    const reads = new Map<File, Promise<{ readonly file: Blob; readonly sha256: string }>>();
    for (const [index, request] of requests.entries()) {
      // The number that the run will have (numbers are given on success only, so that a request that fails leaves no gap) and its time.
      const number = this.nextNumber + index;
      let argv: readonly string[];
      let queryName: string;
      let subjectName: string;
      try {
        // The task of the run's names is the argv's, whatever the names of the inputs: a made input is named after the task is known.
        const probe = buildArgv({ program: request.program, queryName: 'query.fa', subjectName: 'subject.fa', parameters: request.parameters, extra: request.extra ?? [] });
        const task = await this.taskOfArgv(request.program, probe, descriptions);
        const run = { number, program: request.program, task, queuedAt };
        queryName = inputName(request.query, 'query', run);
        subjectName = inputName(request.subject, 'subject', run);
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
        return { ok: false, message: uncheckedOptions(error) };
      }
      if (!start.ok) return start;
      let conffile: Prepared['conffile'];
      const { dryRun, prints } = start.start;
      // The file chosen for the argv's -conffile, under its name (`registryName`: none for `-conffile=`).
      const conffileName = registryName(start.start);
      if (conffileName !== undefined && request.conffile !== undefined) {
        const chosen = request.conffile;
        let read = reads.get(chosen);
        if (read === undefined) reads.set(chosen, (read = this.readConffile(chosen)));
        try {
          const { file, sha256 } = await read;
          conffile = { snapshot: Object.freeze({ name: conffileName, fileName: chosen.name, size: file.size, sha256 }), file };
        } catch (error) {
          return { ok: false, message: unreadableConffile(conffileName, chosen.name, error) };
        }
      }
      let validation: ValidationResult;
      try {
        validation = await this.deps.engine.validate(argv, conffile === undefined ? undefined : { name: conffile.snapshot.name, file: conffile.file });
      } catch (error) {
        return { ok: false, message: uncheckedOptions(error) };
      }
      if (!validation.ok) return validation;
      prepared.push({ request, number, queuedAt, argv, queryName, subjectName, conffile, dryRun, prints, dataLoaders: validation.dataLoaders });
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
    for (const [index, { request, number, queuedAt, argv, queryName, subjectName, conffile, dryRun, prints, dataLoaders }] of prepared.entries()) {
      let query: InputSnapshot;
      let subject: InputSnapshot;
      // The engine's readers of the run's inputs, with the data loaders of its registry (line WQA's D16, 2026-11-20).
      const reader = (role: 'query' | 'subject'): FastaReader => ({ program: request.program, role, dataLoaders });
      try {
        // A run that prints the command line's help or version reads no input (line WQA's D14, 2026-11-20).
        query = prints !== null ? await this.noInput(queryName) : await this.snapshotInput('Query', request.query, queryName, reader('query'), made);
        subject = prints !== null ? await this.noInput(subjectName) : await this.snapshotInput('Subject', request.subject, subjectName, reader('subject'), made);
      } catch (error) {
        return { ok: false, message: failureText(error) };
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
          number,
          program: request.program,
          ...(title ? { title } : {}),
          argv,
          query,
          subject,
          requestedThreads: request.requestedThreads,
          ...(windowMasker === undefined ? {} : { windowMasker }),
          ...(conffile === undefined ? {} : { conffile: conffile.snapshot }),
          ...(dryRun ? { dryRun: true as const } : {}),
          ...(prints === null ? {} : { prints }),
          queuedAt,
          ...(groupId === undefined ? {} : { group: Object.freeze({ groupId, position: index + 1, size: prepared.length }) }),
        }),
      );
      conffiles.push(conffile?.file);
    }
    const refused = await this.holdCounts(
      snapshots.flatMap((snapshot, index) => {
        const taxid = snapshot.windowMasker?.origin?.taxid;
        const file = prepared[index]!.request.windowMasker;
        return taxid === undefined || file === undefined ? [] : [{ runId: snapshot.runId, taxid, file }];
      }),
    );
    if (refused !== undefined) return { ok: false, message: refused };
    // The numbers that the names were made with are taken only now, so that a request that fails leaves no gap.
    this.nextNumber += snapshots.length;
    const numbered = snapshots;
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
   * Cancels the runs of a group: the active one and the queued ones (plan §5.2). Finished
   * runs are left alone; the operation is not atomic.
   */
  cancelGroup(groupId: string): void {
    const members = this.state.get().runs.filter((run) => run.snapshot.group?.groupId === groupId);
    // Queued runs first, so that cancelling the active run does not start the next one of the group.
    for (const run of members) if (run.status === 'queued') this.cancel(run.snapshot.runId);
    for (const run of members) this.cancel(run.snapshot.runId);
  }

  /** Cancels a queued or active run (one that is finalizing too). Finished runs are left alone. */
  cancel(runId: string): void {
    const view = this.find(runId);
    if (view === undefined || isTerminal(view.status)) return;
    const queuedIndex = this.queue.indexOf(runId);
    if (queuedIndex >= 0) {
      this.queue.splice(queuedIndex, 1);
      this.conffiles.delete(runId);
      this.toolDirectories.delete(runId);
      this.deps.sessionFiles?.release(`run:${runId}`);
      this.releaseCounts(runId);
      if (isToolRun(view)) {
        this.updateTool(runId, { status: 'cancelled', record: { ...view.record, endedAt: this.deps.now() } });
      } else {
        this.update(runId, { status: 'cancelled', record: { ...view.record, endedAt: this.deps.now() } });
      }
      return;
    }
    this.cancelRequested.add(runId);
    if (this.active === runId) this.toolAbort?.abort();
    this.deps.engine.cancel(runId);
    if (this.active === runId) this.stopActive?.(new RunCancelledError(runId));
  }

  /**
   * The Data worker ended: the active run cannot keep its results, so it fails now with the error of the Data worker's
   * calls, and its search ends. Before, the search went on, and the run failed only at its commit, when the search
   * ended, while the page already said that the Data worker had stopped (line WQA's sweep, D2; its review, L2, 2026-11-19).
   * The runs queued behind it fail when they start, as each call of the Data worker fails at once.
   */
  private dataWorkerEnded(error: Error): void {
    const runId = this.active;
    if (runId === undefined) return;
    // First the run's wait, so that it ends with this error rather than with the cancel of its search.
    this.stopActive?.(error);
    this.toolAbort?.abort();
    this.deps.engine.cancel(runId);
  }

  /**
   * Saves one compatibility output of a run that kept its outputs (a completed run, or a failed run: the bytes that the
   * CLI leaves before its error), byte for byte, named by outputFileName, with its format's media type.
   */
  async exportOutput(runId: string, output: RunOutput): Promise<void> {
    const view = this.findSearch(runId);
    if (view?.result === undefined) throw new Error('only runs that kept their outputs can be exported');
    const { number, program, argv, queuedAt } = view.snapshot;
    const [bytes, task] = await Promise.all([this.deps.data.readOutput(runId, output), this.taskOf(view.snapshot)]);
    this.deps.downloader.save(outputFileName({ number, program, task, queuedAt, argv }, output), bytes, outputMediaType(argv, output));
  }

  /**
   * Saves a file that a run made besides its outputs (RunRecord.runFiles: the file of `-logfile`, the
   * search strategy file of `-export_search_strategy`), its bytes as the engine wrote them, named by
   * runFileName. A failed run keeps its files too.
   */
  async exportRunFile(runId: string, index: number): Promise<void> {
    const view = this.findSearch(runId);
    const file: RunFile | undefined = view?.record.runFiles?.[index];
    if (view === undefined || file === undefined) throw new Error('the run has no such file');
    const { number, program, argv, queuedAt } = view.snapshot;
    const task = await this.taskOf(view.snapshot);
    this.deps.downloader.save(runFileName({ number, program, task, queuedAt, argv }, file.name, runFileOption(file)), file.bytes, 'text/plain');
  }

  /**
   * Saves the file that a run's query or subject is, where it is not a file of the user's as given (InputSnapshot.from:
   * pasted text, several inputs together, a file with records left out): its bytes, as the run searched them, under the
   * name that the run's command gives them, so that the command runs on it (line WQA's D22, 2026-11-23).
   */
  async exportInput(runId: string, role: 'query' | 'subject'): Promise<void> {
    const input = this.findSearch(runId)?.snapshot[role];
    if (input?.from === undefined) throw new Error(`the run's ${role} is a file of yours as you gave it`);
    this.deps.downloader.save(input.name, input.bytes, 'text/plain');
  }

  /** The input of a run that reads none (one that prints the command line's help or version): no bytes, no records. */
  private async noInput(name: string): Promise<InputSnapshot> {
    const bytes = new Uint8Array(0);
    return Object.freeze({ name, bytes, sha256: await this.deps.digest(bytes), revisionIds: Object.freeze([]), records: Object.freeze([]) });
  }

  /**
   * Saves the text that a run printed in place of a search (RunRecord.printed: the command line's help or version; line
   * WQA's D14, 2026-11-20), its bytes as the engine wrote them, named by printedFileName.
   */
  async exportPrinted(runId: string): Promise<void> {
    const view = this.findSearch(runId);
    const bytes = view?.record.printed;
    const kind = view?.snapshot.prints;
    if (view?.status !== 'completed' || bytes === undefined || kind === undefined) throw new Error('the run printed no text');
    const { number, program, argv, queuedAt } = view.snapshot;
    const task = await this.taskOf(view.snapshot);
    this.deps.downloader.save(printedFileName({ number, program, task, queuedAt, argv }, kind, bytes), bytes, printedFileType(bytes).type);
  }

  /** Reads up to `maxBytes` of a committed tool run's output (`'stdout'`, `'stderr'`, or `'file:N'`). */
  async readToolOutputPreview(
    runId: string,
    key: string,
    maxBytes: number,
  ): Promise<{ readonly bytes: Uint8Array; readonly totalSize: number; readonly truncated: boolean }> {
    const { size } = this.toolOutputEntry(runId, key);
    const end = Math.min(size, Math.max(0, maxBytes));
    const bytes = await this.deps.data.readToolOutputRange(runId, key, 0, end);
    return { bytes, totalSize: size, truncated: size > end };
  }

  /** Returns an OPFS/storage-backed `File` for a committed tool run's output. */
  async getToolOutputFile(runId: string, key: string, fileName?: string): Promise<File> {
    const name = fileName ?? this.toolOutputEntry(runId, key).name;
    return this.deps.data.getToolOutputFile(runId, key, name);
  }

  /** Streams a committed tool run output to a browser download without buffering the full file in page JS memory. */
  async exportToolOutput(runId: string, key: string): Promise<void> {
    const view = this.findTool(runId);
    if (view === undefined || view.record.exitStatus === undefined) throw new Error('the run has no outputs');
    const { name } = this.toolOutputEntry(runId, key);
    const baseName = key === 'stdout' ? 'stdout.txt' : key === 'stderr' ? 'stderr.txt' : name;
    const downloadName = toolOutputFileName(view.snapshot, baseName);
    const file = await this.deps.data.getToolOutputFile(runId, key, downloadName);
    if (this.deps.downloader.saveBlob !== undefined) {
      this.deps.downloader.saveBlob(downloadName, file);
    } else {
      this.deps.downloader.save(downloadName, new Uint8Array(await file.arrayBuffer()), 'application/octet-stream');
    }
  }

  /** Counts files created by completed `-mk_counts` and `-convert` runs in this session, newest first. */
  sessionMadeCounts(): readonly SessionCountsEntry[] {
    const entries: SessionCountsEntry[] = [];
    for (const run of this.state.get().runs) {
      if (!isToolRun(run) || run.status !== 'completed') continue;
      for (const file of run.record.files ?? []) {
        if (file.isCounts === true) {
          entries.push({
            runId: run.snapshot.runId,
            runNumber: run.snapshot.number,
            key: file.key,
            name: file.name,
            size: file.size,
          });
        }
      }
    }
    return entries.reverse();
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
  private taskOf(snapshot: RunSnapshot): Promise<string | undefined> {
    return this.taskOfArgv(snapshot.program, snapshot.argv);
  }

  /** `descriptions` keeps what the engine said of each program, for the requests of one press. */
  private async taskOfArgv(
    program: ProgramId,
    argv: readonly string[],
    descriptions: Map<ProgramId, Promise<ProgramDescription | undefined>> = new Map(),
  ): Promise<string | undefined> {
    let asked = descriptions.get(program);
    if (asked === undefined) descriptions.set(program, (asked = this.deps.engine.describe(program).catch(() => undefined)));
    // Without the engine's default task the name has none; the file is still saved.
    const description = await asked;
    return runTask(argv, description?.parameters ?? []);
  }

  /** Opens one output of a run that kept its outputs without reading it (a binary one, which is not shown): rejects as `readOutput` when it cannot be read. */
  async checkOutput(runId: string, output: RunOutput): Promise<void> {
    if (this.findSearch(runId)?.result === undefined) throw new Error('run has no outputs');
    await this.deps.data.readOutputRange(runId, output, 0, 0);
  }

  /** One output of a run that kept its outputs (a completed run, or a failed run: what the CLI leaves before its error). */
  async readOutput(runId: string, output: RunOutput): Promise<string> {
    if (this.findSearch(runId)?.result === undefined) throw new Error('run has no outputs');
    return new TextDecoder().decode(await this.deps.data.readOutput(runId, output));
  }

  /** The HSP records of a completed run (docs/web/abi_v2.md §8). */
  async readHits(runId: string): Promise<readonly HspRecord[]> {
    if (this.findSearch(runId)?.status !== 'completed') throw new Error('run is not completed');
    return this.deps.data.readHits(runId);
  }

  /** The warnings of a completed run, as the CLI writes them to standard error. */
  async readDiagnostics(runId: string): Promise<string> {
    if (this.findSearch(runId)?.status !== 'completed') throw new Error('run is not completed');
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
    reader: FastaReader,
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
        const revision = await this.deps.data.indexSource(source.sourceId, reader);
        revisionIds = Object.freeze([revision.revisionId]);
        runInput = this.deps.data.buildRunInput(revisionIds);
      }
      const { bytes, sha256, records } = await runInput;
      // Pasted text is a file that LOSAT made, named by the run (`inputName`): the page offers it with the run (line WQA's D22, DECISIONS 44).
      const from = 'dataset' in input ? input.dataset.from : 'text' in input ? [{ name: 'pasted text', records: records.length, leftOut: [], leftOutCount: 0 }] : undefined;
      return Object.freeze({ name, bytes, sha256, revisionIds, records, ...(from === undefined ? {} : { from: Object.freeze(from) }) });
    } catch (error) {
      // A part of LOSAT that could not be loaded or stopped is no fault of the input (line WQA's sweep, D1).
      if (isUnavailable(error)) throw error;
      throw new Error(`${role} FASTA: ${errorMessage(error)}`, { cause: error });
    }
  }

  /**
   * Takes, for runs about to be queued, the shared lock of each taxonomy ID whose saved counts they read, and checks
   * that each such file can still be read (a copy replaced or removed after it was chosen cannot). From then until a
   * run ends, Replace, Delete and a download of that ID refuse (the counts store takes the lock exclusive with
   * `ifAvailable`), in this tab and the others. Returns why the runs cannot be queued, having let go what it took, or
   * undefined when each run holds its locks. Without `countsLocks`, runs hold nothing.
   */
  private async holdCounts(uses: ReadonlyArray<{ readonly runId: string; readonly taxid: number; readonly file: Blob }>): Promise<string | undefined> {
    const locks = this.deps.countsLocks;
    if (locks === undefined || uses.length === 0) return undefined;
    const taken: Array<{ readonly runId: string; readonly release: () => void }> = [];
    const refuse = (message: string): string => {
      for (const { release } of taken) release();
      return message;
    };
    for (const { runId, taxid, file } of uses) {
      const release = await locks.holdShared(taxid);
      if (release === undefined) return refuse(`The saved counts for ${taxid} are being written or removed. Run again when that ends.`);
      taken.push({ runId, release });
      if (!(await readable(file))) return refuse(`The saved counts for ${taxid} were replaced or removed after they were chosen. Run again.`);
    }
    for (const { runId, release } of taken) {
      const before = this.countsHolds.get(runId);
      this.countsHolds.set(
        runId,
        before === undefined
          ? release
          : () => {
              before();
              release();
            },
      );
    }
    return undefined;
  }

  /** Lets go the saved-counts locks of a run that ended or left the queue. */
  private releaseCounts(runId: string): void {
    const release = this.countsHolds.get(runId);
    this.countsHolds.delete(runId);
    release?.();
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
    if (isToolRun(view)) {
      await this.executeTool(runId, view);
      return;
    }
    const { snapshot } = view;
    const conffile = this.conffiles.get(runId);
    let staged = false;
    const startedAt = this.deps.now();
    // The files that the run made besides its outputs, which a completed and a failed run keep, and
    // whether the engine ended the output of a failed run with `fail`, after which the data layer has
    // the warnings that it wrote before its error (a completed run's are in its results).
    let runFiles: readonly RunFile[] | undefined;
    let outputFailed = false;
    // Why the run was stopped (stopActive), if it was.
    let stop: Error | undefined;
    const stopped = new Promise<never>((_, reject) => {
      this.stopActive = (error) => {
        stop ??= error;
        reject(error);
      };
    });
    stopped.catch(() => undefined);
    try {
      this.update(runId, { status: 'preparing', record: { startedAt } });
      const output = await this.deps.data.openRun(runId, runOutputs(snapshot.argv));
      staged = true;
      // A stop while the outputs were staged found no search to end: the engine has not started it, and does not.
      if (stop !== undefined) {
        output.close();
        throw stop;
      }
      const search = this.deps.engine.run(
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
        outputFailed = runOutputFailedOf(error);
        throw error;
      });
      const runtime = await Promise.race([search, stopped]);
      runFiles = runtime.runFiles;
      if (stop !== undefined) throw stop;
      this.update(runId, { status: 'finalizing' });
      const commit = this.deps.data.commitRun(runId);
      const result = await Promise.race([commit, stopped]).catch((error: unknown) => {
        // Results that the commit makes after the run ended are deleted.
        void commit.then(() => this.deps.data.deleteRun(runId)).catch(() => undefined);
        throw error;
      });
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
          messages: runtime.messages,
          ...(runtime.shownIds === undefined ? {} : { shownIds: runtime.shownIds }),
          ...(runtime.windowMasker === undefined
            ? {}
            : { windowMasker: { ...runtime.windowMasker, ...(snapshot.windowMasker?.origin === undefined ? {} : { origin: snapshot.windowMasker.origin }) } }),
          ...(runFiles === undefined ? {} : { runFiles }),
          ...(runtime.printed === undefined ? {} : { printed: runtime.printed }),
        },
      });
    } catch (error) {
      const cancelled = error instanceof RunCancelledError || this.cancelRequested.has(runId);
      // A failed run keeps the outputs that the CLI leaves before its error, and the warnings that it wrote before it
      // (stream 3), which the data layer commits as it commits a completed run's (docs/web/abi_v2.md §5; line WQA's
      // sweep, D3, 2026-11-19, and D15, 2026-11-20). A stop (Cancel, the end of the Data worker) ends that wait, and the
      // other discards are not waited for: a Data worker that does not answer must not keep the run from ending (line
      // WQA's sweep, D2). What was kept after the wait ended is deleted.
      let kept: KeptRun | undefined;
      if (staged && outputFailed && !cancelled && stop === undefined) {
        const keeping = this.deps.data.keepFailedRun(runId);
        kept = await Promise.race([keeping, stopped]).catch(() => {
          void keeping.then((late) => (late === undefined ? undefined : this.deps.data.deleteRun(runId))).catch(() => undefined);
          return undefined;
        });
      } else if (staged) void this.deps.data.discardRun(runId).catch(() => undefined);
      this.update(runId, {
        status: cancelled ? 'cancelled' : 'failed',
        ...(kept?.result === undefined ? {} : { result: kept.result }),
        record: {
          endedAt: this.deps.now(),
          ...(cancelled ? {} : { error: failureText(error) }),
          ...(cancelled || runFiles === undefined ? {} : { runFiles }),
          ...(cancelled || !kept?.diagnostics ? {} : { diagnostics: kept.diagnostics }),
        },
      });
    } finally {
      this.stopActive = undefined;
      // The run's registry file goes with its run; the snapshot keeps its name, size and SHA-256.
      this.conffiles.delete(runId);
      this.releaseCounts(runId);
      this.cancelRequested.delete(runId);
      this.active = undefined;
      void this.refreshStorage();
      this.pump();
    }
  }

  private async executeTool(runId: string, view: ToolRunView): Promise<void> {
    const { snapshot } = view;
    const directory = this.toolDirectories.get(runId) ?? [];
    const startedAt = this.deps.now();
    let stop: Error | undefined;
    const stopped = new Promise<never>((_, reject) => {
      this.stopActive = (error) => {
        stop ??= error;
        reject(error);
      };
    });
    stopped.catch(() => undefined);
    const controller = new AbortController();
    this.toolAbort = controller;
    let staged = false;
    let committing = false;
    let sink: ToolStorageSink | undefined;
    try {
      this.updateTool(runId, { status: 'preparing', record: { startedAt, phaseTimes: { preparing: startedAt } } });
      // Same order as search's openRun: wait for openToolRun instead of racing it so a cancel in
      // preparing always closes the port and deletes the staged tool run on the Data side.
      const port = await this.deps.data.openToolRun(runId);
      staged = true;
      if (stop !== undefined) {
        port.close();
        throw stop;
      }
      sink = new ToolStorageSink(port);
      const runningAt = this.deps.now();
      this.updateTool(runId, { status: 'running', record: { phaseTimes: { preparing: startedAt, running: runningAt } } });
      const runTool = this.deps.engine.runTool;
      if (runTool === undefined) throw new Error('This engine does not support running tools.');
      const runPromise = runTool.call(this.deps.engine, snapshot.argv, {
        directory,
        sinks: sink.sinks,
        signal: controller.signal,
      });
      const result = await Promise.race([runPromise, stopped]);
      if (stop !== undefined) throw stop;
      const finalizingAt = this.deps.now();
      this.updateTool(runId, {
        status: 'finalizing',
        record: { phaseTimes: { preparing: startedAt, running: runningAt, finalizing: finalizingAt } },
      });
      committing = true;
      const commit = sink.commit();
      const summary = await Promise.race([commit, stopped]).catch((error: unknown) => {
        // Outputs sealed after the run ended are deleted (same order as search's commitRun).
        void commit.then(() => this.deps.data.deleteRun(runId)).catch(() => undefined);
        throw error;
      });
      const stdout: ToolOutputEntry = Object.freeze({ key: 'stdout', name: 'stdout', size: summary.stdoutSize });
      const stderr: ToolOutputEntry = Object.freeze({ key: 'stderr', name: 'stderr', size: summary.stderrSize });
      const files: ToolOutputEntry[] = summary.files.map((file) => {
        const option = toolFileOptionLabel(snapshot.argv, file.name);
        const counts = isCountsOutputFile(snapshot.stage, snapshot.argv, file.name);
        return Object.freeze({
          key: file.key,
          name: file.name,
          size: file.size,
          ...(option === undefined ? {} : { option }),
          ...(counts ? { isCounts: true as const } : {}),
        });
      });
      this.updateTool(runId, {
        status: result.status === 0 ? 'completed' : 'failed',
        record: {
          endedAt: this.deps.now(),
          exitStatus: result.status,
          stdout,
          stderr,
          files: Object.freeze(files),
          ...(this.deps.toolRuntime === undefined
            ? {}
            : { runtimePath: this.deps.toolRuntime.path, engineBuild: this.deps.toolRuntime.engineBuild }),
          ...(result.status === 0 ? {} : { error: `Exit status ${result.status}` }),
        },
      });
    } catch (error) {
      const cancelled = error instanceof RunCancelledError || this.cancelRequested.has(runId);
      if (!committing) sink?.fail(error);
      if (staged) void this.deps.data.deleteRun(runId).catch(() => undefined);
      this.updateTool(runId, {
        status: cancelled ? 'cancelled' : 'failed',
        record: {
          endedAt: this.deps.now(),
          ...(cancelled ? {} : { error: failureText(error) }),
        },
      });
    } finally {
      this.stopActive = undefined;
      this.toolAbort = undefined;
      this.toolDirectories.delete(runId);
      this.deps.sessionFiles?.release(`run:${runId}`);
      this.releaseCounts(runId);
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

  private find(runId: string): AnyRunView | undefined {
    return this.state.get().runs.find((run) => run.snapshot.runId === runId);
  }

  private findSearch(runId: string): SearchRunView | undefined {
    const view = this.find(runId);
    return view !== undefined && !isToolRun(view) ? view : undefined;
  }

  private findTool(runId: string): ToolRunView | undefined {
    const view = this.find(runId);
    return view !== undefined && isToolRun(view) ? view : undefined;
  }

  private toolOutputEntry(runId: string, key: string): ToolOutputEntry {
    const view = this.findTool(runId);
    if (view === undefined || view.record.exitStatus === undefined) throw new Error('the run has no outputs');
    if (key === 'stdout' && view.record.stdout !== undefined) return view.record.stdout;
    if (key === 'stderr' && view.record.stderr !== undefined) return view.record.stderr;
    const file = view.record.files?.find((entry) => entry.key === key);
    if (file === undefined) throw new Error(`the run has no output ${key}`);
    return file;
  }

  private update(runId: string, change: Partial<Omit<SearchRunView, 'snapshot'>>): void {
    // A run that ends lets go of the saved counts before it is shown as ended, so that a Replace pressed then is not refused.
    if (change.status !== undefined && isTerminal(change.status)) this.releaseCounts(runId);
    this.setRuns(
      this.state.get().runs.map((run) =>
        run.snapshot.runId === runId && !isToolRun(run)
          ? { ...run, ...change, record: { ...run.record, ...change.record } }
          : run,
      ),
    );
  }

  private updateTool(runId: string, change: Partial<Omit<ToolRunView, 'snapshot'>>): void {
    if (change.status !== undefined && isTerminal(change.status)) this.releaseCounts(runId);
    this.setRuns(
      this.state.get().runs.map((run) =>
        run.snapshot.runId === runId && isToolRun(run)
          ? { ...run, ...change, record: { ...run.record, ...change.record } }
          : run,
      ),
    );
  }

  private setRuns(runs: readonly AnyRunView[]): void {
    this.state.set({ ...this.state.get(), runs: Object.freeze([...runs]) });
  }
}

/** A request that the engine validated, with its argv, the names of its inputs, and what the engine's start read of it. */
interface Prepared {
  readonly request: SearchRequest;
  /** The run's number and the time of its "Run LOSAT" (the names of its made inputs have them). */
  readonly number: number;
  readonly queuedAt: number;
  readonly argv: readonly string[];
  readonly queryName: string;
  readonly subjectName: string;
  /** The registry file of the argv's `-conffile`, read from the chosen file: what the snapshot keeps of it, and its bytes. */
  readonly conffile: { readonly snapshot: ConffileSnapshot; readonly file: Blob } | undefined;
  readonly dryRun: boolean;
  /** What the run prints in place of a search, as the engine's start says; null for a search. */
  readonly prints: PrintedKind | null;
  /** Whether the run reads its inputs with NCBI's data loaders on (the validation's verdict: `[BLAST] DATA_LOADERS` of its registry). */
  readonly dataLoaders: boolean;
}

/** What one enqueueAll made: sources of text and files, and the cached run inputs of dataset revisions. */
interface Made {
  readonly sources: string[];
  readonly runInputs: string[];
}

/**
 * The name of an input in the run's argv (line WQA's D22, DECISIONS 44): the user's file name for the user's file as given
 * (D21), and for every file that the page makes (pasted text, several files together, a file with records left out) the run's
 * own prefix and the role (`inputFileName`), which the page saves the file under too.
 */
function inputName(input: SequenceInput, role: 'query' | 'subject', run: Parameters<typeof inputFileName>[0]): string {
  if ('file' in input) return input.file.name;
  if ('dataset' in input && input.dataset.from === undefined) {
    if (input.dataset.name === undefined) throw new Error(`the ${role} is neither a file of the user's nor made of parts`);
    return input.dataset.name;
  }
  return inputFileName(run, role);
}

function engineInput(input: InputSnapshot): EngineInput {
  return { bytes: input.bytes, sha256: input.sha256, revisionIds: input.revisionIds, records: input.records };
}

function errorMessage(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}

/** Why the engine could not check the options of a request: a part of LOSAT that is unavailable is said as the page says it. */
function uncheckedOptions(error: unknown): string {
  return isUnavailable(error) ? failureText(error) : `The options could not be checked: ${errorMessage(error)}`;
}

/**
 * Whether a file chosen earlier can still be read: its last byte is read. A file of the browser's storage that was
 * written or removed since it was opened refuses to be read (or, where the browser reads the current file, comes
 * shorter).
 */
async function readable(file: Blob): Promise<boolean> {
  if (file.size === 0) return true;
  try {
    return (await file.slice(file.size - 1, file.size).arrayBuffer()).byteLength === 1;
  } catch {
    return false;
  }
}
