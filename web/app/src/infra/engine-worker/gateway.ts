// The Wasm EngineGateway (plan §3.1, §4.6, §4.7, §5.5): the main-thread side of the Engine
// worker. It chooses the threads of a search, starts an Engine worker (a runtime
// generation) when there is none, posts the search with the Data worker's output port,
// and relays the phases. A cancel ends the Engine worker with its thread workers; the next
// search starts a new generation, whose instance registers the subject again from the run
// snapshot's bytes (R1). Messages of an older generation are dropped. A runtime is also
// renewed before a search when its instance ran too many searches or its memory grew too
// much (policy.ts). `describe` and `validate` go to the Data worker's serial reactor, which
// can answer while a search runs (src/infra/reactor/control.ts).
import {
  InputMismatchError,
  RunCancelledError,
  type EngineGateway,
  type EngineInput,
  type EnginePhase,
  type EngineRunRequest,
  type ProgramDescription,
  type RuntimeInfo,
  type ValidationResult,
} from '../../ports/engine';
import { engineBuildName, type EngineAssets, type ReactorAsset } from '../reactor/assets';
import type { EngineControl } from '../reactor/control';
import { compileReactor, fetchReactor } from '../reactor/instance';
import { chooseThreads, DEFAULT_RENEWAL, renewalReason, type RenewalLimits } from './policy';
import type { EngineCommand, EngineEvent, EngineInit, ModuleSource, ThreadFault, WorkerError, WorkerInput, WorkerRun } from './protocol';

export interface WasmEngineOptions {
  readonly assets: EngineAssets;
  /** describe/validate: the Data worker's serial reactor. */
  readonly control: EngineControl;
  readonly hardwareConcurrency?: number;
  readonly renewal?: Partial<RenewalLimits>;
  /** Starts an Engine worker (tests may start a test build). */
  readonly createWorker?: () => Worker;
  /** Loads the two modules; by default they are fetched, checked and compiled (tests replace it). */
  readonly loadModules?: () => Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }>;
}

interface Generation {
  readonly id: number;
  readonly worker: Worker;
  readonly faults: BroadcastChannel;
  /** Set after a search whose instance must be renewed before the next one. */
  renew: string | undefined;
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

  validate(argv: readonly string[]): Promise<ValidationResult> {
    return this.options.control.validate(argv);
  }

  run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo> {
    if (this.active !== undefined) return Promise.reject(new Error('the engine runs one search at a time'));
    return new Promise<RuntimeInfo>((resolve, reject) => {
      const active: ActiveRun = { runId: request.runId, onPhase, resolve, reject, generation: undefined, settled: false };
      this.active = active;
      this.start(active, request, output).catch((error: unknown) => {
        output.close();
        this.settle(active, error instanceof Error ? error : new Error(String(error)));
      });
    });
  }

  cancel(runId: string): void {
    const active = this.active;
    if (active === undefined || active.runId !== runId || active.settled) return;
    // Ending the worker ends its thread workers too; nothing of this generation is used again.
    if (active.generation !== undefined && active.generation === this.current) this.retire();
    this.settle(active, new RunCancelledError(runId));
  }

  private async start(active: ActiveRun, request: EngineRunRequest, output: MessagePort): Promise<void> {
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
      output,
    };
    generation.worker.postMessage(command, [output]);
  }

  /** Fetches, checks and compiles both modules once; each generation gets the compiled ones. */
  private compiled(): Promise<{ readonly serial: ModuleSource; readonly threads: ModuleSource }> {
    if (this.modules === undefined) {
      const load = async (asset: ReactorAsset): Promise<ModuleSource> => compileReactor(await fetchReactor(asset));
      const modules =
        this.options.loadModules?.() ??
        Promise.all([load(this.options.assets.serial), load(this.options.assets.threads)]).then(([serial, threads]) => ({
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
      new Worker(new URL('./engine-worker.ts', import.meta.url), { type: 'module', name: 'losat-engine' });
    const faults = new BroadcastChannel(`losat-web:engine-faults:${crypto.randomUUID()}`);
    const generation: Generation = { id, worker, faults, renew: undefined };
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
    generation.renew = renewalReason(result.memory, this.renewal);
    this.settle(active, undefined, {
      path: result.path,
      threads: result.threads,
      engineBuild: result.engineBuild,
      ...(result.fallbackReason === undefined ? {} : { fallbackReason: result.fallbackReason }),
      runtimeGeneration: generation.id,
      memory: result.memory,
      subjectRetained: result.subjectRetained,
    });
  }

  /** The runtime stopped (a trap in a thread, or the worker itself); its search fails. */
  private fail(generation: Generation, message: string): void {
    if (generation !== this.current) return;
    this.retire();
    const active = this.active;
    if (active !== undefined && active.generation === generation) this.settle(active, new Error(message));
  }

  /** Ends the current generation: its worker and, with it, its thread workers. */
  private retire(): void {
    const generation = this.current;
    if (generation === undefined) return;
    this.current = undefined;
    generation.worker.onmessage = null;
    generation.worker.onerror = null;
    generation.worker.onmessageerror = null;
    generation.worker.terminate();
    generation.faults.close();
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

function workerInput(program: string, input: EngineInput): WorkerInput {
  return { bytes: input.bytes, records: input.records, key: inputKey(program, input) };
}

/** Rebuilds an error of the worker with its name and message (as src/infra/data-worker/rpc.ts). */
function toError(error: WorkerError): Error {
  if (error.name === 'InputMismatchError' && error.role !== undefined && error.detail !== undefined) {
    return new InputMismatchError(error.role, error.detail);
  }
  const rebuilt = new Error(error.message);
  rebuilt.name = error.name;
  return rebuilt;
}
