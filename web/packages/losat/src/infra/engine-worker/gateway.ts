// The Wasm EngineGateway (plan §3.1, §4.6, §4.7, §5.5): the main-thread side of the Engine
// worker. It chooses the threads of a search, starts an Engine worker (a runtime
// generation) when there is none, posts the search with the Data worker's output port,
// and relays the phases. A cancel ends the Engine worker with its thread workers; the next
// search starts a new generation, whose instance registers the subject again from the run
// snapshot's bytes (R1). Messages of an older generation are dropped. A worker that ends by
// itself fails its search: with an error event, or by the release of the Web Lock that it holds
// for as long as it lives, when it ends without one (worker-life.ts; line WQA's sweep, D2; where the browser refuses
// that lock, nothing tells it). A runtime is also
// renewed before a search when its instance ran too many searches or its memory grew too
// much (policy.ts). A preparation hint (prepare) starts the generation, its instance and its
// thread workers, or that renewal, before the search that would start them. `describe`, `start`
// and `validate` go to the Data worker's serial reactor, which can answer while a search runs
// (src/infra/reactor/control.ts). A tool run (`runTool`, ABI v2 version 13) is posted to the same worker, and the chunks and file
// frames that the worker posts as the tool writes them go to a ToolOutputRouter (tool-output.ts), which gives them to the
// caller's sinks or collects them as Blobs: the Gateway holds no output whole (D-WW3-6), and where the page can share memory
// with the worker the worker waits while the sinks are behind (tool-credit.ts, D-WW3-11). `validateTool` runs on the Engine
// worker too, with the run's directory (D-WW3-10): the Gateway keeps the checks that are pending and posts them again to
// the next worker when the one they went to ends.
import type { CountsFile } from '../../domain/counts.js';
import type { DirectoryFile, RegistryFile } from '../../domain/run-files.js';
import { parseToolDescription, type ToolDescription } from '../../domain/tool-description.js';
import {
  InputMismatchError,
  RunCancelledError,
  UnavailableError,
  withRunFiles,
  withRunOutputFailed,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type EnginePrepareHint,
  type EngineRunRequest,
  type ProgramDescription,
  type RuntimeInfo,
  type StartResult,
  type ToolEngine,
  type ToolRunOptions,
  type ToolRunResult,
  type ToolValidateOptions,
  type ToolValidation,
  type ValidationResult,
} from '../../ports/engine.js';
import { engineBuildName, type EngineAssets } from '../reactor/assets.js';
import type { EngineControl, ToolControl } from '../reactor/control.js';
import { loadReactor } from '../reactor/instance.js';
import { watchLife, webLocks, type LifeLocks } from '../worker-life.js';
import { InstanceUse } from './instance-use.js';
import { answerLink } from './link.js';
import { chooseThreads, DEFAULT_RENEWAL, limitsHolding, renewalReason, type RenewalLimits } from './policy.js';
import type {
  EngineCommand,
  EngineEvent,
  EngineInit,
  KeeperGrow,
  ModuleSource,
  WorkerError,
  WorkerCounts,
  WorkerInput,
  WorkerLink,
  WorkerRun,
  WorkerTool,
  WorkerToolValidate,
} from './protocol.js';
import { createToolCredit, ToolCreditPage } from './tool-credit.js';
import { ToolOutputRouter } from './tool-output.js';

/**
 * How a runtime ended (retire): at once (its threaded instance was not in use, or its worker had
 * ended), or once the keeper grew the shared memory.
 */
export type RuntimeEnd = 'at-once' | 'grown';

export interface WasmEngineOptions {
  readonly assets: EngineAssets;
  /** describe/start/validate (and optional describeTool): the Data worker's serial reactor (`reactorControl`). `validateTool` is not here: it runs on the Engine worker. */
  readonly control: EngineControl & Partial<ToolControl>;
  readonly hardwareConcurrency?: number;
  readonly renewal?: Partial<RenewalLimits>;
  /** Starts an Engine worker (tests may start a test build). */
  readonly createWorker?: () => Worker;
  /** Loads the two modules; by default they are fetched, checked and compiled (tests replace it). */
  readonly loadModules?: () => Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }>;
  /** Told how each runtime ended (tests record it). */
  readonly onRuntimeEnd?: (end: RuntimeEnd) => void;
  /** The Web Locks with which the Gateway learns that a worker ended (worker-life.ts); the browser's by default (tests replace them). */
  readonly locks?: LifeLocks;
}

/** Why a search fails whose Engine worker could not start: its script or its first message (line WQA's sweep, D1). */
export const ENGINE_NOT_STARTED = 'LOSAT could not start the part that runs searches';
/**
 * Why a search fails whose Engine worker ended without an error event (line WQA's sweep, D2). The next search starts
 * another worker, which the sentence says: it said only what stopped (the second review of D1 and D2, L2, 2026-11-19).
 */
export const ENGINE_ENDED = 'The part of LOSAT that runs searches stopped. Run the search again.';

interface Generation {
  readonly id: number;
  readonly worker: Worker;
  /** Set after a search whose instance must be renewed before the next one. */
  renew: string | undefined;
  /** The threads and argv of the last preparation posted since the last search (prepare). */
  prepared: string | undefined;
  /**
   * The link of the keeper of the threaded instance's shared memory (memory-keeper.ts, link.ts): the last that the Engine
   * worker posted, the keeper of the memory in use (runtime.ts, load). The main thread holds no shared memory: a
   * reference here would keep the memory until this thread collects it.
   */
  keeper: MessagePort | undefined;
  /**
   * The main thread's ends of the links of the generation's current ThreadHost: its keeper's and its thread workers'
   * (link.ts). They are closed when another host links its keeper, and when the generation ends.
   */
  readonly links: MessagePort[];
  /**
   * Whether the threaded instance is in use, read and marked synchronously (instance-use.ts,
   * retire). Undefined where the page cannot share memory: no threaded instance starts then.
   */
  readonly use: InstanceUse | undefined;
  /** The Web Lock that the worker holds for as long as it lives (EngineInit.lifeLock). */
  readonly lifeLock: string;
  /** Ends the wait for that lock when the Gateway ends the worker itself (end). */
  readonly life: AbortController;
  /** Whether the generation has ended (end): its worker, its links and the keeper are gone. */
  ended: boolean;
  /**
   * Whether the worker said that it started (`started`, the answer to `init`): an error before it is a worker that could
   * not start, and one after it a worker that stopped, also during its first search, whose error read as "could not
   * start" while the worker said so only once it held its lock (line WQA's sweep, D2; its review, M1, 2026-11-19).
   */
  started: boolean;
}

interface ActiveRun {
  readonly runId: string;
  readonly onPhase: (phase: EnginePhase) => void;
  /** Takes the run's `RuntimeInfo`, or a tool run's `ToolRunResult` (`tool` is set). */
  readonly resolve: (value: RuntimeInfo | ToolRunResult) => void;
  readonly reject: (error: Error) => void;
  /** A tool run: where the worker's chunks and frames go. */
  readonly tool?: ToolOutputRouter;
  /** A tool run where the page can share memory with the worker: the counter of the bytes that the page has not taken. */
  readonly credit?: { readonly page: ToolCreditPage; readonly buffer: SharedArrayBuffer };
  generation: Generation | undefined;
  settled: boolean;
  /** A tool run whose status arrived: it waits for its sinks only, and the end of its worker no longer takes its result. */
  finishing?: boolean;
}

/** A check of a tool's words that waits for its answer (`WasmEngine.validateTool`). */
interface PendingValidation {
  readonly id: string;
  readonly argv: readonly string[];
  readonly files: readonly DirectoryFile[] | undefined;
  readonly resolve: (verdict: ToolValidation) => void;
  readonly reject: (error: Error) => void;
  /** The generation that it was posted to; undefined before the post, and after that generation ended (it is posted again). */
  generation: number | undefined;
  /**
   * How many workers that held it ended through no run's doing (a stop of the instance, a crash with no run on the worker):
   * cancels, renewals and the ends that a run caused are not counted.
   */
  ends: number;
}

/** How many such ends a check survives: a worker that ends each time it is asked would otherwise be started again for good. */
const VALIDATION_ENDS = 3;

export class WasmEngine implements EngineGateway, ToolEngine {
  private current: Generation | undefined;
  private active: ActiveRun | undefined;
  private generations = 0;
  private modules: Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }> | undefined;
  private readonly renewal: RenewalLimits;
  private readonly locks: LifeLocks | undefined;
  private readonly validations = new Map<string, PendingValidation>();

  constructor(private readonly options: WasmEngineOptions) {
    this.renewal = { ...DEFAULT_RENEWAL, ...options.renewal };
    this.locks = options.locks ?? webLocks();
  }

  describe(program: ProgramDescription['program']): Promise<ProgramDescription> {
    return this.options.control.describe(program);
  }

  async describeTool(tool: string): Promise<ToolDescription> {
    if (this.options.control.describeTool !== undefined) return this.options.control.describeTool(tool);
    const { status, stdout, stderr } = await this.runTool([tool, '-xmlhelp']);
    if (status !== 0) throw new Error((await stderr.text()).trim() || `${tool} -xmlhelp exited with status ${status}`);
    return parseToolDescription(await stdout.text());
  }

  start(argv: readonly string[]): Promise<StartResult> {
    return this.options.control.start(argv);
  }

  validate(argv: readonly string[], registry?: RegistryFile): Promise<ValidationResult> {
    return this.options.control.validate(argv, registry);
  }

  run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    if (this.active !== undefined) return Promise.reject(new Error('the engine runs one search at a time'));
    return new Promise<RuntimeInfo>((resolve, reject) => {
      const active: ActiveRun = {
        runId: request.runId,
        onPhase,
        resolve: (info) => resolve(info as RuntimeInfo),
        reject,
        generation: undefined,
        settled: false,
      };
      this.active = active;
      this.launch(active, request, output).catch((error: unknown) => {
        output.close();
        this.settle(active, error instanceof Error ? error : new Error(String(error)));
      });
    });
  }

  /**
   * Checks the words of a tool (`windowmasker`) as the command line does before the tool opens a file (`ToolEngine.validateTool`),
   * on the Engine worker's loaded instance (the serial one when none is: it disturbs no runtime that the searches use), with
   * `options.directory` opened there as for `runTool` (the check reads the run's registry files). It takes its turn in the
   * worker like a run: behind the run in progress, if there is one. It makes the engine no busier: a run may start while it is
   * pending, and it may be asked during a run. When the worker that it went to ends (a cancel, a renewal, a stop), it is posted
   * again to the next one; a stop of its own instance ends the worker as a run's does (`validated`).
   */
  validateTool(argv: readonly string[], options: ToolValidateOptions = {}): Promise<ToolValidation> {
    return new Promise<ToolValidation>((resolve, reject) => {
      const id = `validate:${crypto.randomUUID()}`;
      this.validations.set(id, { id, argv, files: options.directory, resolve, reject, generation: undefined, ends: 0 });
      this.postValidations();
    });
  }

  /**
   * Runs a tool (`windowmasker`) on the Engine worker's serial module, in the run's working directory that
   * `options.directory` makes (`ToolEngine.runTool`). The worker posts the tool's standard output, standard error and files
   * as it writes them; each goes to `options.sinks` or is collected as a Blob (tool-output.ts). A search and a tool run take
   * turns: one at a time. `options.signal` cancels the run as `cancel` does a search.
   */
  runTool(argv: readonly string[], options: ToolRunOptions = {}): Promise<ToolRunResult> {
    if (this.active !== undefined) return Promise.reject(new Error('the engine runs one search at a time'));
    const runId = `tool:${crypto.randomUUID()}`;
    const { signal } = options;
    if (signal?.aborted === true) return Promise.reject(new RunCancelledError(runId));
    return new Promise<ToolRunResult>((resolve, reject) => {
      const stop = () => this.cancel(runId);
      const finish = () => signal?.removeEventListener('abort', stop);
      // Where the page can share memory with the worker, the worker waits while the sinks are behind (tool-credit.ts).
      const buffer = createToolCredit();
      const page = buffer === undefined ? undefined : new ToolCreditPage(buffer);
      // A sink or the frames that fail end the run with their error, as a cancel ends a search.
      const router: ToolOutputRouter = new ToolOutputRouter(options.sinks, (error) => this.stopTool(active, error), page);
      const active: ActiveRun = {
        runId,
        onPhase: () => undefined,
        resolve: (result) => {
          finish();
          resolve(result as ToolRunResult);
        },
        reject: (error) => {
          finish();
          reject(error);
        },
        generation: undefined,
        settled: false,
        tool: router,
        ...(buffer === undefined || page === undefined ? {} : { credit: { page, buffer } }),
      };
      this.active = active;
      signal?.addEventListener('abort', stop, { once: true });
      this.launchTool(active, argv, options.directory).catch((error: unknown) => this.settle(active, error instanceof Error ? error : new Error(String(error))));
    });
  }

  /**
   * Starts the runtime of a search like `hint` before it comes (EngineGateway.prepare): the
   * modules, the Engine worker, the instance and the thread workers that the search would
   * start, a renewal that the next search would make, and with the hint's argv the warm-up
   * of that search (warm-up.ts). Nothing happens while a search runs.
   */
  prepare(hint: EnginePrepareHint): void {
    if (this.active !== undefined) return;
    const threads = chooseThreads(
      hint.requestedThreads,
      hint.inputBytes,
      this.options.hardwareConcurrency ?? navigator.hardwareConcurrency,
    );
    this.compiled()
      .then((modules) => {
        // A search that started meanwhile prepares for itself.
        if (this.active !== undefined) return;
        if (this.current?.renew !== undefined) this.retire();
        const generation = this.current ?? this.startGeneration(modules);
        const prepared = [threads, ...(hint.argv ?? [])].join('\0');
        if (generation.prepared === prepared) return;
        generation.prepared = prepared;
        this.post(generation.worker, {
          type: 'prepare',
          generation: generation.id,
          threads,
          ...(hint.argv === undefined ? {} : { argv: hint.argv }),
        });
      })
      // The next search loads the modules (or starts the worker) again and reports the failure.
      .catch(() => undefined);
  }

  cancel(runId: string): void {
    const active = this.active;
    if (active === undefined || active.runId !== runId || active.settled) return;
    // Ending the worker ends its thread workers too; nothing of this generation is used again.
    if (active.generation !== undefined && active.generation === this.current) this.retire();
    this.settle(active, new RunCancelledError(runId));
  }

  private async launchTool(active: ActiveRun, argv: readonly string[], files: ToolRunOptions['directory']): Promise<void> {
    const modules = await this.compiled();
    if (active.settled) return;
    if (this.current?.renew !== undefined) this.retire();
    const generation = this.current ?? this.startGeneration(modules);
    active.generation = generation;
    generation.prepared = undefined;
    const command: WorkerTool = {
      type: 'tool',
      generation: generation.id,
      runId: active.runId,
      argv,
      ...(files === undefined ? {} : { files }),
      ...(active.credit === undefined ? {} : { credit: active.credit.buffer }),
    };
    this.post(generation.worker, command);
  }

  /** Posts the checks that wait for a worker (`validateTool`; those whose worker ended) to the current generation, starting one if needed. */
  private postValidations(): void {
    this.compiled().then(
      (modules) => {
        const waiting = [...this.validations.values()].filter((validation) => validation.generation === undefined);
        if (waiting.length === 0) return;
        let generation: Generation;
        try {
          generation = this.current ?? this.startGeneration(modules);
        } catch (error) {
          // A worker that cannot be started (a refused `new Worker`, a throwing `createWorker`): every check that waits for one
          // fails with the reason, as a run does (`launchTool`); a later call asks again.
          for (const validation of waiting) this.dropValidation(validation, toThrown(error));
          return;
        }
        for (const validation of waiting) {
          validation.generation = generation.id;
          const command: WorkerToolValidate = {
            type: 'tool-validate',
            generation: generation.id,
            id: validation.id,
            argv: validation.argv,
            ...(validation.files === undefined ? {} : { files: validation.files }),
          };
          try {
            this.post(generation.worker, command);
          } catch (error) {
            this.dropValidation(validation, toThrown(error));
          }
        }
      },
      (error: unknown) => {
        // The modules did not load: the checks that wait fail with the reason (a later call loads them again).
        for (const validation of [...this.validations.values()].filter((entry) => entry.generation === undefined)) {
          this.dropValidation(validation, toThrown(error));
        }
      },
    );
  }

  private dropValidation(validation: PendingValidation, error: Error): void {
    this.validations.delete(validation.id);
    validation.reject(error);
  }

  /**
   * The answer of the worker to a check: its verdict, or its error. A stopped instance ends the worker as it does after a run
   * (`done` with `stopped`, `fail`): the instance of a threaded runtime keeps its thread workers and keeper until the Gateway
   * ends the worker (retire), so a check and a run are handled alike. The run that waits behind the check on that worker
   * fails with ENGINE_ENDED (it never ran), and the other checks go to the next worker.
   */
  private validated(generation: Generation, event: Extract<EngineEvent, { type: 'tool-validated' | 'tool-validate-failed' }>): void {
    const validation = this.validations.get(event.id);
    if (validation === undefined || validation.generation !== generation.id) return;
    this.validations.delete(event.id);
    if (event.type === 'tool-validated') {
      validation.resolve(event.verdict);
      return;
    }
    if (event.error.stopped === true) this.fail(generation, new Error(ENGINE_ENDED), true);
    validation.reject(toError(event.error));
  }

  /** A tool run whose sink failed: its worker ends (the tool may still be running), and the run rejects with the sink's error. */
  private stopTool(active: ActiveRun, error: unknown): void {
    if (active.settled) return;
    if (active.generation !== undefined && active.generation === this.current) this.retire();
    this.settle(active, error instanceof Error ? error : new Error(String(error)));
  }

  private async launch(active: ActiveRun, request: EngineRunRequest, output: MessagePort): Promise<void> {
    // A run searches its subject and the files of its directory (the database of -db, the repeat database of
    // BLASTN's -filtering_db), by their size.
    const searched = (request.subject?.bytes.length ?? 0) + (request.files ?? []).reduce((total, { file }) => total + file.size, 0);
    const threads = chooseThreads(
      request.requestedThreads,
      (request.query?.bytes.length ?? 0) + searched,
      this.options.hardwareConcurrency ?? navigator.hardwareConcurrency,
    );
    const modules = await this.compiled();
    if (active.settled) {
      output.close();
      return;
    }
    if (this.current?.renew !== undefined) this.retire();
    const generation = this.current ?? this.startGeneration(modules);
    active.generation = generation;
    generation.prepared = undefined;
    const program = request.argv[0] ?? '';
    const command: WorkerRun = {
      type: 'run',
      generation: generation.id,
      runId: request.runId,
      program,
      argv: request.argv,
      threads,
      ...(request.query === undefined ? {} : { query: workerInput(program, request.query) }),
      ...(request.subject === undefined ? {} : { subject: workerInput(program, request.subject) }),
      ...(request.files === undefined ? {} : { files: request.files }),
      ...(request.windowMasker === undefined ? {} : { windowMasker: workerCounts(request.windowMasker) }),
      ...(request.conffile === undefined ? {} : { conffile: request.conffile }),
      output,
    };
    generation.worker.postMessage(command, [output]);
  }

  /** Fetches, checks and compiles both modules once; each generation gets the compiled ones. */
  private compiled(): Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }> {
    if (this.modules === undefined) {
      const modules =
        this.options.loadModules?.() ??
        Promise.all([loadReactor(this.options.assets.serial), loadReactor(this.options.assets.threads)]).then(([serial, threads]) => ({
          serial,
          threads,
        }));
      this.modules = modules;
      // A failed load is tried again by the next search.
      modules.catch(() => {
        if (this.modules === modules) this.modules = undefined;
      });
    }
    return this.modules;
  }

  private startGeneration(modules: { readonly serial: ModuleSource; readonly threads: ModuleSource }): Generation {
    const id = ++this.generations;
    const worker =
      this.options.createWorker?.() ??
      new Worker(new URL('./engine-worker.js', import.meta.url), { type: 'module', name: 'losat-engine' });
    // A page that cannot share memory with its workers has no SharedArrayBuffer, or cannot post one.
    const use = typeof SharedArrayBuffer === 'function' && globalThis.crossOriginIsolated !== false ? InstanceUse.create() : undefined;
    const generation: Generation = {
      id,
      worker,
      keeper: undefined,
      links: [],
      renew: undefined,
      prepared: undefined,
      use,
      lifeLock: `losat-web:engine-life:${crypto.randomUUID()}`,
      life: new AbortController(),
      ended: false,
      started: false,
    };
    worker.onmessage = (event: MessageEvent<EngineEvent>) => this.receive(generation, event.data);
    worker.onerror = (event) => {
      event.preventDefault();
      // A script that could not be loaded gives an error event without a message in Chromium, Firefox and WebKit, and the
      // search failed with "The engine worker stopped: an error without a message" (line WQA's sweep, D1).
      this.fail(
        generation,
        generation.started
          ? new Error(`The engine worker stopped: ${event.message || 'an error without a message'}`)
          : new UnavailableError(ENGINE_NOT_STARTED),
      );
    };
    worker.onmessageerror = () => this.fail(generation, new Error('The engine worker sent a message that could not be read'));
    const { assets } = this.options;
    const init: EngineInit = {
      type: 'init',
      generation: id,
      modules,
      builds: { serial: engineBuildName('serial', assets.serial), threads: engineBuildName('threads', assets.threads) },
      lifeLock: generation.lifeLock,
      ...(use === undefined ? {} : { use: use.buffer }),
    };
    try {
      this.post(worker, init);
    } catch {
      // A browser that cannot post a compiled module: the worker fetches and compiles.
      this.post(worker, { ...init, modules: { serial: { asset: assets.serial }, threads: { asset: assets.threads } } });
    }
    this.current = generation;
    return generation;
  }

  private post(worker: Worker, command: EngineCommand): void {
    worker.postMessage(command);
  }

  private receive(generation: Generation, event: EngineEvent): void {
    if (!('generation' in event) || event.generation !== generation.id) return;
    if (event.type === 'life') {
      // The worker holds its lock now, and no preparation or search of it has started: the page may wait for the lock,
      // also when the generation is already ending (a search cancelled before this message came).
      if (event.held && this.locks !== undefined && !generation.ended) {
        watchLife(this.locks, generation.lifeLock, generation.life.signal, () => this.lifeEnded(generation));
      }
      return;
    }
    if (event.type === 'link') {
      this.link(generation, event);
      return;
    }
    if (generation !== this.current) return;
    if (event.type === 'started') {
      generation.started = true;
      return;
    }
    if (event.type === 'fault') {
      this.fail(generation, new Error(event.message));
      return;
    }
    if (event.type === 'tool-validated' || event.type === 'tool-validate-failed') {
      this.validated(generation, event);
      return;
    }
    const active = this.active;
    if (active === undefined || active.settled || active.generation !== generation || active.runId !== event.runId) return;
    if (event.type === 'phase') {
      active.onPhase(event.phase);
      return;
    }
    if (event.type === 'tool-chunk') {
      if (event.stream === 'stdout') active.tool?.stdout(event.bytes);
      else active.tool?.stderr(event.bytes);
      return;
    }
    if (event.type === 'tool-file') {
      active.tool?.frame(event.frame);
      return;
    }
    if (event.type === 'tool-done') {
      const router = active.tool;
      if (router === undefined) return;
      active.finishing = true;
      const { status, memory, countsBytes } = event.result;
      generation.renew = renewalReason(memory, limitsHolding(this.renewal, countsBytes));
      // The result comes once the router has handed every chunk and frame that the worker posted before it to the sinks.
      router.finish().then(
        (collected) => this.settle(active, undefined, { status, ...collected }),
        (error: unknown) => this.settle(active, error instanceof Error ? error : new Error(String(error))),
      );
      return;
    }
    if (event.type !== 'done') return;
    if (!event.ok) {
      if (event.error.stopped === true) this.retire();
      this.settle(active, toError(event.error));
      return;
    }
    const { result } = event;
    generation.renew = renewalReason(result.memory, limitsHolding(this.renewal, result.countsBytes));
    this.settle(active, undefined, {
      path: result.path,
      threads: result.threads,
      engineBuild: result.engineBuild,
      ...(result.fallbackReason === undefined ? {} : { fallbackReason: result.fallbackReason }),
      runtimeGeneration: generation.id,
      memory: result.memory,
      subjectRetained: result.subjectRetained,
      messages: result.messages,
      ...(result.shownIds === undefined ? {} : { shownIds: result.shownIds }),
      ...(result.windowMasker === undefined ? {} : { windowMasker: result.windowMasker }),
      ...(result.runFiles === undefined ? {} : { runFiles: result.runFiles }),
      ...(result.printed === undefined ? {} : { printed: result.printed }),
      ...(result.readRecords === undefined ? {} : { readRecords: result.readRecords }),
    });
  }

  /**
   * Answers the main thread's end of a link of a thread worker or the keeper (link.ts): a thread that stopped fails the
   * search, and the keeper's answer to the grow ends the generation (retire). Also while the generation is ending, which
   * may wait for that answer; the links close when it has ended (end).
   *
   * A ThreadHost links its keeper first, and the runtime makes a host only once the previous one has ended with its
   * thread workers and keeper (runtime.ts, load: a search after one on the other module, a failed preparation). The
   * links before a keeper's lead to workers that are gone: they are closed then. They stayed open until the generation
   * ended, which with Auto threads alternating small and large searches kept up to 32 more for each switch for the
   * whole session (the review of line WQA's D8, L1, 2026-11-20).
   */
  private link(generation: Generation, event: WorkerLink): void {
    const { port } = event;
    if (event.role === 'keeper') {
      for (const earlier of generation.links.splice(0)) earlier.close();
      generation.keeper = port;
    }
    generation.links.push(port);
    answerLink(port, (message) => {
      if (message?.type === 'fault') this.fail(generation, new Error(message.message));
      else if (message?.type === 'grown') this.end(generation, 'grown');
    });
  }

  /**
   * The runtime stopped (a trap in a thread, the worker itself, or the instance of a check); its search fails. The end counts
   * against the checks that the worker held unless a run is on it (the run is the likelier cause) or `byCheck` says that a
   * check stopped it (VALIDATION_ENDS).
   */
  private fail(generation: Generation, error: Error, byCheck = false): void {
    if (generation !== this.current) return;
    // A worker that never started cannot answer a check: it fails with the reason, rather than asking the next worker for good.
    if (!generation.started) {
      for (const validation of [...this.validations.values()]) if (validation.generation === generation.id) this.dropValidation(validation, error);
    }
    const active = this.active;
    const run = active !== undefined && active.generation === generation ? active : undefined;
    this.retire(byCheck || run === undefined);
    if (run !== undefined && !run.finishing) this.settle(run, error);
  }

  /**
   * The worker ended (the release of its lock), and its thread workers and its keeper with it: no thread can grow the
   * memory, and no keeper answers, so the generation ends at once, whether it was current (its search fails) or
   * already ending, waiting for the keeper after a cancel, a renewal, a thread fault, a stopped instance or an error of
   * the worker. Only a current generation was told before: one already ending waited for the keeper for good, with its
   * two channels open (the second and third reviews of line WQA's D1 and D2, L1 and L2, 2026-11-19).
   */
  private lifeEnded(generation: Generation): void {
    this.end(generation, 'at-once');
    this.fail(generation, new Error(ENGINE_ENDED));
  }

  /**
   * Ends the current generation: its worker and, with it, its thread workers and the keeper.
   * While the threaded instance is in use (a call of the Engine worker into it, a search or a
   * warm-up; a thread of its search pool that runs; the grow with which the Engine worker ends
   * thread workers whose threads run), a thread may be growing the shared memory, and WebKit
   * crashes the page when workers with instances on a shared memory end while a thread grows it
   * (growToMaximum). The Gateway marks the generation retiring and reads whether the instance is
   * in use in one atomic operation (instance-use.ts): if it is not, the worker ends at once, and
   * no call or thread can start before it has ended. If it is, the generation ends once the
   * keeper has grown the memory to its maximum, where no thread can grow it. The keeper is itself
   * a thread that grows the memory, so the worker does not end before its answer, however long
   * the grow takes: on a busy machine it took up to 2 s, longer than the 1 s that the Gateway
   * once waited before it ended the worker without the grow. Once the generation is retiring, the
   * Engine worker ends none of its thread workers and not the keeper, also after the instance
   * stopped (ThreadHost.terminate, runtime.ts dispose), and the keeper answers also when the grow
   * fails (memory-keeper.ts). Until then the ended search may still run and write to its output
   * port, whose Data worker side the coordinator discards. The page waits for the worker's lock until
   * the generation has ended: a worker that ends meanwhile ends it at once (lifeEnded).
   */
  private retire(endsChecks = false): void {
    const generation = this.current;
    if (generation === undefined) return;
    this.current = undefined;
    // The checks of tools that went to this worker are posted again to the next (validateTool); `endsChecks`: the end counts
    // against them (it is the checks' doing, or no run's: a cancel or a renewal of a run does not).
    for (const validation of [...this.validations.values()]) {
      if (validation.generation !== generation.id) continue;
      validation.generation = undefined;
      if (endsChecks) validation.ends++;
      if (validation.ends >= VALIDATION_ENDS) this.dropValidation(validation, new Error(ENGINE_ENDED));
    }
    if (this.validations.size > 0) this.postValidations();
    // A worker that ended took its threads and its keeper with it (lifeEnded).
    if (generation.ended) return;
    generation.worker.onerror = null;
    generation.worker.onmessageerror = null;
    // A use starts only once the keeper's link has had this thread's answer (thread-host.ts, prepare): a runtime in use
    // has its keeper here.
    const { keeper } = generation;
    if (generation.use?.retire() !== true || keeper === undefined) {
      this.end(generation, 'at-once');
      return;
    }
    const grow: KeeperGrow = { type: 'grow' };
    keeper.postMessage(grow);
  }

  /** Ends a generation once: its worker, with its thread workers and the keeper, and its links. */
  private end(generation: Generation, how: RuntimeEnd): void {
    if (generation.ended) return;
    generation.ended = true;
    // The Gateway ends this worker itself: the release of its lock is no news.
    generation.life.abort();
    generation.worker.onmessage = null;
    generation.worker.onerror = null;
    generation.worker.onmessageerror = null;
    generation.worker.terminate();
    for (const port of generation.links) port.close();
    try {
      this.options.onRuntimeEnd?.(how);
    } catch {
      // An observer does not change how the runtime ends.
    }
  }

  private settle(active: ActiveRun, error: Error | undefined, info?: RuntimeInfo | ToolRunResult): void {
    if (active.settled) return;
    active.settled = true;
    // A worker that waits for the page to take bytes goes on (tool-credit.ts): the run is over.
    active.credit?.page.close();
    if (this.active === active) this.active = undefined;
    if (error !== undefined) {
      // A tool run that ends without its status: the sinks of its open files are told.
      void active.tool?.abort(error);
      active.reject(error);
    } else active.resolve(info!);
  }
}

/** The retained-subject key: the program, the SHA-256 and the dataset revisions (plan §4.6). */
export function inputKey(program: string, input: EngineInput): string {
  return JSON.stringify([program, input.sha256, input.revisionIds]);
}

/** The key of a counts file: its name, size and last-modified time (the SHA-256 is known only after the engine read it). */
export function countsKey(counts: CountsFile): string {
  return JSON.stringify([counts.name, counts.size, counts.lastModified]);
}

function workerCounts(counts: CountsFile): WorkerCounts {
  return { name: counts.name, blob: counts.file, key: countsKey(counts) };
}

function workerInput(program: string, input: EngineInput): WorkerInput {
  return {
    bytes: input.bytes,
    records: input.records,
    key: inputKey(program, input),
    ...(input.assembled === true ? { assembled: true } : {}),
  };
}

/** A thrown value as an Error. */
function toThrown(error: unknown): Error {
  return error instanceof Error ? error : new Error(String(error));
}

/**
 * Rebuilds an error of the worker with its name and message (as web/app/src/infra/data-worker/rpc.ts),
 * the files of a failed run (`runFilesOf`), and whether it ended the run's output with `fail`
 * (`runOutputFailedOf`).
 */
function toError(error: WorkerError): Error {
  if (error.name === 'InputMismatchError' && error.role !== undefined && error.detail !== undefined) {
    return new InputMismatchError(error.role, error.detail);
  }
  const rebuilt = new Error(error.message);
  rebuilt.name = error.name;
  const failed = error.outputFailed === true ? withRunOutputFailed(rebuilt) : rebuilt;
  return error.runFiles === undefined ? failed : withRunFiles(failed, error.runFiles);
}
