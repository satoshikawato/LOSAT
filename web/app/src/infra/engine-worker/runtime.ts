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
import { RunOutputWriter } from '../run-output/writer';
import type { RuntimePath, WorkerRun, WorkerRunResult } from './protocol';
import { ThreadHost } from './thread-host';

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
      this.loaded = { path, instance: await instantiateSerial(await this.options.serial()), runs: 0, retained: undefined };
      return this.loaded;
    }
    const unsupported = threadedSupport();
    if (unsupported !== undefined) throw new Error(unsupported);
    const { module, artifact } = await this.options.threads();
    let memory: WebAssembly.Memory;
    try {
      // Exactly the limits that the module declares; the host never enlarges them (plan TD-7).
      memory = new WebAssembly.Memory({
        initial: artifact.memory.initial,
        maximum: artifact.memory.maximum!,
        shared: true,
      });
    } catch (error) {
      const gib = (artifact.memory.maximum! * 65_536) / 2 ** 30;
      throw new Error(`the browser could not reserve ${gib} GiB of shared memory (${messageOf(error)})`, { cause: error });
    }
    const threads = new ThreadHost({ module, memory, faultChannel: this.options.faultChannel });
    try {
      const instance = await instantiateThreadedMain(module, memory, threads.spawn);
      this.loaded = { path, instance, threads, runs: 0, retained: undefined };
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
