// The Wasm EngineGateway (plan §3.1, §4.6, §4.7, §5.5): the main-thread side of the Engine
// worker. It chooses the threads of a search, starts an Engine worker (a runtime
// generation) when there is none, posts the search with the Data worker's output port,
// and relays the phases. A cancel ends the Engine worker with its thread workers; the next
// search starts a new generation, whose instance registers the subject again from the run
// snapshot's bytes (R1). Messages of an older generation are dropped. A runtime is also
// renewed before a search when its instance ran too many searches or its memory grew too
// much (policy.ts). A preparation hint (prepare) starts the generation, its instance and its
// thread workers, or that renewal, before the search that would start them. `describe`, `start`
// and `validate` go to the Data worker's serial reactor, which can answer while a search runs
// (src/infra/reactor/control.ts).
import type { CountsFile } from '../../domain/counts.js';
import type { RegistryFile } from '../../domain/run-files.js';
import {
  InputMismatchError,
  RunCancelledError,
  withRunFiles,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type EnginePrepareHint,
  type EngineRunRequest,
  type ProgramDescription,
  type RuntimeInfo,
  type StartResult,
  type ValidationResult,
} from '../../ports/engine.js';
import { engineBuildName, type EngineAssets } from '../reactor/assets.js';
import type { EngineControl } from '../reactor/control.js';
import { loadReactor } from '../reactor/instance.js';
import { InstanceUse } from './instance-use.js';
import { chooseThreads, DEFAULT_RENEWAL, limitsHolding, renewalReason, type RenewalLimits } from './policy.js';
import type {
  EngineCommand,
  EngineEvent,
  EngineInit,
  KeeperGrow,
  KeeperGrown,
  ModuleSource,
  ThreadFault,
  WorkerError,
  WorkerCounts,
  WorkerInput,
  WorkerRun,
} from './protocol.js';

/**
 * How a runtime ended (retire): at once (its threaded instance was not in use), or once the
 * keeper grew the shared memory.
 */
export type RuntimeEnd = 'at-once' | 'grown';

export interface WasmEngineOptions {
  readonly assets: EngineAssets;
  /** describe/start/validate: the Data worker's serial reactor. */
  readonly control: EngineControl;
  readonly hardwareConcurrency?: number;
  readonly renewal?: Partial<RenewalLimits>;
  /** Starts an Engine worker (tests may start a test build). */
  readonly createWorker?: () => Worker;
  /** Loads the two modules; by default they are fetched, checked and compiled (tests replace it). */
  readonly loadModules?: () => Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }>;
  /** Told how each runtime ended (tests record it). */
  readonly onRuntimeEnd?: (end: RuntimeEnd) => void;
}

interface Generation {
  readonly id: number;
  readonly worker: Worker;
  readonly faults: BroadcastChannel;
  /** Set after a search whose instance must be renewed before the next one. */
  renew: string | undefined;
  /** The threads and argv of the last preparation posted since the last search (prepare). */
  prepared: string | undefined;
  /**
   * The channel to the keeper of the threaded instance's shared memory (memory-keeper.ts).
   * The main thread holds no shared memory: a reference here would keep the memory until this
   * thread collects it.
   */
  readonly keeper: BroadcastChannel;
  /**
   * Whether the threaded instance is in use, read and marked synchronously (instance-use.ts,
   * retire). Undefined where the page cannot share memory: no threaded instance starts then.
   */
  readonly use: InstanceUse | undefined;
}

interface ActiveRun {
  readonly runId: string;
  readonly onPhase: (phase: EnginePhase) => void;
  readonly resolve: (info: RuntimeInfo) => void;
  readonly reject: (error: Error) => void;
  generation: Generation | undefined;
  settled: boolean;
}

export class WasmEngine implements EngineGateway {
  private current: Generation | undefined;
  private active: ActiveRun | undefined;
  private generations = 0;
  private modules: Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }> | undefined;
  private readonly renewal: RenewalLimits;

  constructor(private readonly options: WasmEngineOptions) {
    this.renewal = { ...DEFAULT_RENEWAL, ...options.renewal };
  }

  describe(program: ProgramDescription['program']): Promise<ProgramDescription> {
    return this.options.control.describe(program);
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
      const active: ActiveRun = { runId: request.runId, onPhase, resolve, reject, generation: undefined, settled: false };
      this.active = active;
      this.launch(active, request, output).catch((error: unknown) => {
        output.close();
        this.settle(active, error instanceof Error ? error : new Error(String(error)));
      });
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

  private async launch(active: ActiveRun, request: EngineRunRequest, output: MessagePort): Promise<void> {
    const threads = chooseThreads(
      request.requestedThreads,
      request.query.bytes.length + request.subject.bytes.length,
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
      query: workerInput(program, request.query),
      subject: workerInput(program, request.subject),
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
    const faults = new BroadcastChannel(`losat-web:engine-faults:${crypto.randomUUID()}`);
    const keeper = new BroadcastChannel(`losat-web:engine-keeper:${crypto.randomUUID()}`);
    // A page that cannot share memory with its workers has no SharedArrayBuffer, or cannot post one.
    const use = typeof SharedArrayBuffer === 'function' && globalThis.crossOriginIsolated !== false ? InstanceUse.create() : undefined;
    const generation: Generation = { id, worker, faults, keeper, renew: undefined, prepared: undefined, use };
    worker.onmessage = (event: MessageEvent<EngineEvent>) => this.receive(generation, event.data);
    worker.onerror = (event) => {
      event.preventDefault();
      this.fail(generation, `The engine worker stopped: ${event.message || 'an error without a message'}`);
    };
    worker.onmessageerror = () => this.fail(generation, 'The engine worker sent a message that could not be read');
    faults.onmessage = (event: MessageEvent<ThreadFault>) => this.fail(generation, event.data.message);
    const { assets } = this.options;
    const init: EngineInit = {
      type: 'init',
      generation: id,
      modules,
      builds: { serial: engineBuildName('serial', assets.serial), threads: engineBuildName('threads', assets.threads) },
      faultChannel: faults.name,
      keeperChannel: keeper.name,
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
    if (generation !== this.current || !('generation' in event) || event.generation !== generation.id) return;
    const active = this.active;
    if (active === undefined || active.settled || active.generation !== generation || active.runId !== event.runId) return;
    if (event.type === 'phase') {
      active.onPhase(event.phase);
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
      ...(result.windowMasker === undefined ? {} : { windowMasker: result.windowMasker }),
      ...(result.runFiles === undefined ? {} : { runFiles: result.runFiles }),
    });
  }

  /** The runtime stopped (a trap in a thread, or the worker itself); its search fails. */
  private fail(generation: Generation, message: string): void {
    if (generation !== this.current) return;
    this.retire();
    const active = this.active;
    if (active !== undefined && active.generation === generation) this.settle(active, new Error(message));
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
   * port, whose Data worker side the coordinator discards.
   */
  private retire(): void {
    const generation = this.current;
    if (generation === undefined) return;
    this.current = undefined;
    generation.worker.onmessage = null;
    generation.worker.onerror = null;
    generation.worker.onmessageerror = null;
    const { keeper } = generation;
    let ended = false;
    const end = (how: RuntimeEnd) => {
      if (ended) return;
      ended = true;
      generation.worker.terminate();
      generation.faults.close();
      keeper.close();
      try {
        this.options.onRuntimeEnd?.(how);
      } catch {
        // An observer does not change how the runtime ends.
      }
    };
    if (generation.use?.retire() !== true) {
      end('at-once');
      return;
    }
    keeper.onmessage = (event: MessageEvent<KeeperGrown>) => {
      if (event.data?.type === 'grown') end('grown');
    };
    const grow: KeeperGrow = { type: 'grow' };
    keeper.postMessage(grow);
  }

  private settle(active: ActiveRun, error: Error | undefined, info?: RuntimeInfo): void {
    if (active.settled) return;
    active.settled = true;
    if (this.active === active) this.active = undefined;
    if (error !== undefined) active.reject(error);
    else active.resolve(info!);
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
  return { bytes: input.bytes, records: input.records, key: inputKey(program, input) };
}

/**
 * Rebuilds an error of the worker with its name and message (as web/app/src/infra/data-worker/rpc.ts),
 * and the files of a failed run (`runFilesOf`).
 */
function toError(error: WorkerError): Error {
  if (error.name === 'InputMismatchError' && error.role !== undefined && error.detail !== undefined) {
    return new InputMismatchError(error.role, error.detail);
  }
  const rebuilt = new Error(error.message);
  rebuilt.name = error.name;
  return error.runFiles === undefined ? rebuilt : withRunFiles(rebuilt, error.runFiles);
}
