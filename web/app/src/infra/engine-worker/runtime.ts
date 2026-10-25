// The engine inside the Engine worker (plan §3.1, §4.6, §4.7): one warm reactor instance,
// serial or threaded, the subject it holds (R1), and one search at a time.
//
// A search: choose the module (threaded for more than one thread, unless threaded
// searches are not possible here, then serial with one thread and the reason), append
// -num_threads and validate that argv (ABI v2 §7), register the query and, unless the
// instance holds it, the subject, compare the records that the engine read with the
// record tables (plan §5.4; a difference stops the run before any output), run, and send
// the output streams to the Data worker's port.
import { recordMismatch, type RecordKey } from '../../domain/dataset';
import { InputMismatchError, type EnginePhase } from '../../ports/engine';
import { OUTPUT_STREAMS, type OutputStream } from '../../ports/run-output';
import { EngineCallError, ROLE_QUERY, ROLE_SUBJECT } from '../reactor/abi';
import type { ReactorArtifact } from '../reactor/artifact';
import { instantiateSerial, instantiateThreadedMain, type ReactorInstance } from '../reactor/instance';
import { createSharedMemory } from '../reactor/shared-memory';
import { RunOutputWriter } from '../run-output/writer';
import type { RuntimePath, WorkerRun, WorkerRunResult } from './protocol';
import { ThreadHost } from './thread-host';
import { warmUp } from './warm-up';

export interface RuntimeOptions {
  readonly serial: () => Promise<WebAssembly.Module>;
  readonly threads: () => Promise<{ readonly module: WebAssembly.Module; readonly artifact: ReactorArtifact }>;
  readonly builds: { readonly serial: string; readonly threads: string };
  readonly faultChannel: string;
}

interface Retained {
  readonly key: string;
  readonly handle: number;
  readonly records: readonly RecordKey[];
}

interface Loaded {
  readonly path: RuntimePath;
  readonly instance: ReactorInstance;
  readonly threads?: ThreadHost;
  runs: number;
  retained: Retained | undefined;
  /** The argvs (without -num_threads) of the searches and warm-ups that ran on this instance. */
  readonly ran: Set<string>;
}

/** Why this worker cannot run the threaded module, or undefined if it may. */
export function threadedSupport(): string | undefined {
  if (globalThis.crossOriginIsolated !== true) {
    return 'this page is not cross-origin isolated, so the browser does not allow memory shared between threads';
  }
  if (typeof SharedArrayBuffer !== 'function') return 'this browser has no SharedArrayBuffer';
  if (typeof Atomics?.wait !== 'function') return 'this browser cannot wait on shared memory in a worker';
  if (typeof Worker !== 'function') return 'this browser cannot start workers from a worker';
  return undefined;
}

export class EngineRuntime {
  private loaded: Loaded | undefined;
  private threadedUnavailable: string | undefined;

  constructor(private readonly options: RuntimeOptions) {}

  async run(request: WorkerRun, phase: (phase: EnginePhase) => void): Promise<WorkerRunResult> {
    phase('preparing');
    // A search of more than one thread that cannot run threaded runs on the serial module
    // with -num_threads 1 and the reason (plan §4.7).
    let fallbackReason: string | undefined;
    let loaded: Loaded | undefined;
    if (request.threads > 1) fallbackReason = this.threadedUnavailable;
    if (request.threads > 1 && fallbackReason === undefined) {
      try {
        loaded = await this.load('threaded');
      } catch (error) {
        // Not possible in this browser or page: later searches of this runtime run serially too.
        this.threadedUnavailable = `the threaded engine could not start: ${messageOf(error)}`;
        fallbackReason = this.threadedUnavailable;
      }
      try {
        await loaded?.threads!.prepare(request.threads - 1);
      } catch (error) {
        // A thread worker that did not become ready: the next search tries with new ones.
        fallbackReason = `the threaded engine could not start: ${messageOf(error)}`;
        this.dispose();
        loaded = undefined;
      }
    }
    loaded ??= await this.load('serial');
    const threads = loaded.path === 'threaded' ? request.threads : 1;
    const argv = [...request.argv, '-num_threads', String(threads)];
    const abi = loaded.instance.abi;
    const invalid = abi.validate(argv);
    if (invalid !== undefined) throw new EngineCallError(invalid);

    const query = abi.register(request.program, ROLE_QUERY, request.query.bytes);
    let queryReleased = false;
    const releaseQuery = () => {
      if (queryReleased || abi.stoppedBy !== undefined) return;
      queryReleased = true;
      abi.release(query.handle);
    };
    try {
      const queryMismatch = recordMismatch(request.query.records, query.records);
      if (queryMismatch !== undefined) throw new InputMismatchError('query', queryMismatch);
      const subjectRetained = loaded.retained?.key === request.subject.key;
      const subject = this.subject(loaded, request);
      const subjectMismatch = recordMismatch(request.subject.records, subject.records);
      if (subjectMismatch !== undefined) throw new InputMismatchError('subject', subjectMismatch);

      const writer = new RunOutputWriter(request.output);
      const linearBytesBefore = abi.memoryBytes();
      // The engine's standard error of this search only, for a stop message.
      loaded.instance.wasi.clear();
      phase('running');
      abi.run(argv, query.handle, subject.handle, (stream, bytes) => {
        if (OUTPUT_STREAMS.includes(stream as OutputStream)) writer.write(stream as OutputStream, bytes);
      });
      const linearBytesAfter = abi.memoryBytes();
      loaded.runs++;
      loaded.ran.add(request.argv.join('\0'));
      // The last call into the instance comes before `end`: a failure after it must not
      // follow an output that the data layer may already commit.
      releaseQuery();
      phase('finalizing');
      writer.end();
      const path = loaded.path;
      return {
        path,
        threads,
        engineBuild: path === 'threaded' ? this.options.builds.threads : this.options.builds.serial,
        ...(fallbackReason === undefined ? {} : { fallbackReason }),
        memory: { linearBytesBefore, linearBytesAfter, instanceRuns: loaded.runs },
        subjectRetained,
      };
    } finally {
      if (abi.stoppedBy === undefined) releaseQuery();
      else this.dispose();
    }
  }

  /**
   * Starts what a search of `threads` threads starts before it runs: the module's instance
   * and, for the threaded module, its thread workers. With the search's `argv` (program
   * first, without -num_threads), it then runs that search once on one thread on built-in
   * sequences (warm-up.ts) unless one like it ran on the instance. It never ends a working
   * instance of the other module: a queued search may still use it and the subject it holds.
   * A failure is thrown; `run` loads and prepares again and handles it as before. The Engine
   * worker calls `prepare` and `run` one after another, never together.
   */
  async prepare(threads: number, argv?: readonly string[]): Promise<void> {
    const path: RuntimePath = threads > 1 && this.threadedUnavailable === undefined ? 'threaded' : 'serial';
    const current = this.loaded;
    if (current !== undefined && current.path !== path && current.instance.abi.stoppedBy === undefined) return;
    if (path === 'serial') {
      const loaded = await this.load('serial');
      if (argv !== undefined) this.warmUp(loaded, argv);
      return;
    }
    let loaded: Loaded;
    try {
      loaded = await this.load('threaded');
    } catch (error) {
      // As in run: later searches of this runtime run serially, with this reason.
      this.threadedUnavailable = `the threaded engine could not start: ${messageOf(error)}`;
      throw error;
    }
    try {
      await loaded.threads!.prepare(threads - 1);
    } catch (error) {
      // The search starts new thread workers rather than meeting the failed one.
      this.dispose();
      throw error;
    }
    if (argv !== undefined) this.warmUp(loaded, argv);
  }

  /** The warm-up of `argv`, unless a search or warm-up like it ran on the instance. */
  private warmUp(loaded: Loaded, argv: readonly string[]): void {
    const key = argv.join('\0');
    if (loaded.ran.has(key)) return;
    loaded.ran.add(key);
    warmUp(loaded.instance.abi, argv);
    loaded.instance.wasi.clear();
    // The next search loads the module again.
    if (loaded.instance.abi.stoppedBy !== undefined) this.dispose();
  }

  /** Ends the instance and its thread workers. */
  dispose(): void {
    this.loaded?.threads?.terminate();
    this.loaded = undefined;
  }

  /** The subject handle: the retained one if its key matches, else a new registration. */
  private subject(loaded: Loaded, request: WorkerRun): Retained {
    const abi = loaded.instance.abi;
    if (loaded.retained?.key === request.subject.key) return loaded.retained;
    if (loaded.retained !== undefined) {
      abi.release(loaded.retained.handle);
      loaded.retained = undefined;
    }
    const registered = abi.register(request.program, ROLE_SUBJECT, request.subject.bytes);
    loaded.retained = { key: request.subject.key, handle: registered.handle, records: registered.records };
    return loaded.retained;
  }

  private async load(path: RuntimePath): Promise<Loaded> {
    const current = this.loaded;
    if (current?.path === path && current.instance.abi.stoppedBy === undefined) return current;
    this.dispose();
    if (path === 'serial') {
      const instance = await instantiateSerial(await this.options.serial());
      this.loaded = { path, instance, runs: 0, retained: undefined, ran: new Set() };
      return this.loaded;
    }
    const unsupported = threadedSupport();
    if (unsupported !== undefined) throw new Error(unsupported);
    const { module, artifact } = await this.options.threads();
    // The maximum that the module declares, or a smaller one where the browser cannot
    // reserve it; never a larger one (plan TD-7).
    const { memory } = createSharedMemory(artifact.memory);
    const threads = new ThreadHost({ module, memory, faultChannel: this.options.faultChannel });
    try {
      const instance = await instantiateThreadedMain(module, memory, threads.spawn);
      this.loaded = { path, instance, threads, runs: 0, retained: undefined, ran: new Set() };
      return this.loaded;
    } catch (error) {
      threads.terminate();
      throw error;
    }
  }
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
