// The engine inside the Engine worker (plan §3.1, §4.6, §4.7): one warm reactor instance,
// serial or threaded, the subject it holds (R1), and one search at a time.
//
// A search: choose the module (threaded for more than one thread, unless threaded
// searches are not possible here, then serial with one thread and the reason), append
// -num_threads and validate that argv (ABI v2 §7), register the query and, unless the
// instance holds it, the subject, compare the records that the engine read with the
// record tables (plan §5.4; a difference stops the run before any output), run, and send
// the output streams to the Data worker's port. A BLASTN search with -window_masker_db also
// gives the engine the counts file, which the instance holds for later BLASTN searches.
import { recordMismatch, type RecordKey } from '../../domain/dataset';
import { InputMismatchError, type EnginePhase } from '../../ports/engine';
import { OUTPUT_STREAMS, type OutputStream } from '../../ports/run-output';
import type { CountsRecord } from '../../domain/run';
import { EngineCallError, FILE_WINDOW_MASKER_COUNTS, ROLE_QUERY, ROLE_SUBJECT, type ChunkBuffer, type FileInfo } from '../reactor/abi';
import type { ReactorArtifact } from '../reactor/artifact';
import { instantiateSerial, instantiateThreadedMain, type ReactorInstance } from '../reactor/instance';
import { createSharedMemory } from '../reactor/shared-memory';
import { RunOutputWriter } from '../run-output/writer';
import type { InstanceUse } from './instance-use';
import type { RuntimePath, WorkerCounts, WorkerRun, WorkerRunResult } from './protocol';
import { ThreadHost } from './thread-host';
import { warmUp } from './warm-up';

export interface RuntimeOptions {
  readonly serial: () => Promise<WebAssembly.Module>;
  readonly threads: () => Promise<{ readonly module: WebAssembly.Module; readonly artifact: ReactorArtifact }>;
  readonly builds: { readonly serial: string; readonly threads: string };
  readonly faultChannel: string;
  /** The channel of the keeper of the threaded instance's shared memory (memory-keeper.ts). */
  readonly keeperChannel: string;
  /**
   * Whether the threaded instance is in use, shared with the Gateway (instance-use.ts): the
   * calls into it and the threads of its search pool count, and the Gateway ends this worker
   * during a use only after the keeper's grow (gateway.ts, retire). Undefined where the page
   * cannot share memory with a worker: the threaded module does not start then.
   */
  readonly use: InstanceUse | undefined;
  /** Test builds only: called before each output of a search, which it may stop by throwing. */
  readonly beforeOutput?: () => void;
  /** Test builds only: called at the start of each use of the threaded instance (a search, a warm-up). */
  readonly onUse?: () => void;
}

interface Retained {
  readonly key: string;
  readonly handle: number;
  readonly records: readonly RecordKey[];
}

/** The WindowMasker counts file that the instance holds: at most one (DESIGN_ww1 §5). */
interface HeldCounts {
  readonly key: string;
  readonly handle: number;
  readonly info: FileInfo;
}

/** Bytes of a counts file given to the engine at a time, through one buffer of engine memory. */
export const COUNTS_CHUNK_BYTES = 8 * 1024 * 1024;

interface Loaded {
  readonly path: RuntimePath;
  readonly instance: ReactorInstance;
  readonly threads?: ThreadHost;
  runs: number;
  retained: Retained | undefined;
  /** The counts file of BLASTN's -window_masker_db that the instance holds. */
  counts: HeldCounts | undefined;
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
    // The calls into the threaded instance are a use from the first until the last has
    // returned; its threads count themselves (ThreadHost.spawn). The search is synchronous.
    const use = loaded.path === 'threaded' ? this.options.use : undefined;
    // The counts file first: `validate` (in `search`) checks it where the command line does
    // (after -dust, before the query range and the other options; ABI v2 §4), and the command
    // line reads it before the queries. A counts file that is released or loaded here is not
    // held together with the inputs of this run, and a replaced one is not held with its
    // successor. Its calls into the instance are uses of their own, with the file's chunks read
    // between them.
    let windowMasker: CountsRecord | undefined;
    try {
      windowMasker = await this.counts(loaded, request, use);
    } catch (error) {
      // A trap while the file was given to the engine: the next search loads the module again
      // (a stopped threaded instance is ended by the Gateway, as in `search`).
      if (loaded.instance.abi.stoppedBy !== undefined && loaded.path === 'serial') this.dispose();
      throw error;
    }
    if (use !== undefined && !use.enter()) throw new Error(ENDING);
    try {
      if (use !== undefined) this.options.onUse?.();
      return this.search(loaded, request, fallbackReason, phase, windowMasker);
    } finally {
      use?.leave();
    }
  }

  /** The search of `run` on the instance that `run` chose. */
  private search(
    loaded: Loaded,
    request: WorkerRun,
    fallbackReason: string | undefined,
    phase: (phase: EnginePhase) => void,
    windowMasker: CountsRecord | undefined,
  ): WorkerRunResult {
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
        this.options.beforeOutput?.();
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
        ...(windowMasker === undefined ? {} : { windowMasker }),
        ...(loaded.counts === undefined ? {} : { countsBytes: loaded.counts.info.size }),
      };
    } finally {
      if (abi.stoppedBy === undefined) releaseQuery();
      // A stopped threaded instance keeps its thread workers and the keeper: threads of the
      // search may still run and grow the memory. They count as a use, so the Gateway, which
      // ends the runtime on the stop, ends this worker once the keeper has grown the memory
      // (gateway.ts, retire).
      else if (loaded.path === 'serial') this.dispose();
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
    // Its calls into the threaded instance are a use, as a search's are (run).
    const use = loaded.path === 'threaded' ? this.options.use : undefined;
    if (use !== undefined && !use.enter()) return;
    try {
      if (use !== undefined) this.options.onUse?.();
      loaded.ran.add(key);
      warmUp(loaded.instance.abi, argv);
      loaded.instance.wasi.clear();
    } finally {
      use?.leave();
    }
    // The next search loads the module again.
    if (loaded.instance.abi.stoppedBy !== undefined) this.dispose();
  }

  /**
   * Ends the instance and its thread workers, and returns true; false, ending nothing, once the
   * Gateway ends this runtime: it may be waiting for the answer of the keeper, which ends with
   * the thread workers (ThreadHost.terminate).
   */
  dispose(): boolean {
    const threads = this.loaded?.threads;
    if (threads !== undefined ? !threads.terminate() : this.options.use?.retiring() === true) return false;
    this.loaded = undefined;
    return true;
  }

  /**
   * The counts file of the search, by the rules of retention (DESIGN_ww1 §5): the instance holds
   * at most one; the same one (key) is used again; it is released before a search of another
   * program (the memory goes back to that search), before a different one is loaded (never two at
   * once), and with the instance (`dispose`). A BLASTN search without -window_masker_db keeps
   * it. A file is given to the engine in the 'preparing' phase only. Returns what the search did
   * with the file, if it has one.
   */
  private async counts(loaded: Loaded, request: WorkerRun, use: InstanceUse | undefined): Promise<CountsRecord | undefined> {
    const wanted = request.windowMasker;
    if (request.program !== 'blastn') {
      this.releaseCounts(loaded, use);
      return undefined;
    }
    if (wanted === undefined) return undefined;
    const key = `${wanted.key} ${await edgeDigest(wanted.blob)}`;
    const held = loaded.counts;
    if (held?.key === key) {
      return { name: wanted.name, size: held.info.size, sha256: held.info.sha256, retained: true, loadMs: 0 };
    }
    this.releaseCounts(loaded, use);
    const started = performance.now();
    const info = await this.loadCounts(loaded, wanted, key, use);
    return { name: wanted.name, size: info.size, sha256: info.sha256, retained: false, loadMs: Math.round(performance.now() - started) };
  }

  private releaseCounts(loaded: Loaded, use: InstanceUse | undefined): void {
    const held = loaded.counts;
    if (held === undefined) return;
    loaded.counts = undefined;
    inUse(use, () => loaded.instance.abi.fileRelease(held.handle));
  }

  /** Gives the file to the engine in chunks, ends it, and holds it. */
  private async loadCounts(loaded: Loaded, wanted: WorkerCounts, key: string, use: InstanceUse | undefined): Promise<FileInfo> {
    const { abi } = loaded.instance;
    const size = wanted.blob.size;
    const handle = inUse(use, () => abi.fileBegin(FILE_WINDOW_MASKER_COUNTS, wanted.name, size));
    let buffer: ChunkBuffer | undefined;
    try {
      buffer = inUse(use, () => abi.chunkBuffer(Math.max(1, Math.min(COUNTS_CHUNK_BYTES, size))));
      for (let offset = 0; offset < size; offset += COUNTS_CHUNK_BYTES) {
        const bytes = new Uint8Array(await wanted.blob.slice(offset, offset + COUNTS_CHUNK_BYTES).arrayBuffer());
        inUse(use, () => abi.fileChunk(handle, bytes, buffer!));
      }
      const info = inUse(use, () => abi.fileEnd(handle));
      loaded.counts = { key, handle, info };
      return info;
    } catch (error) {
      // The engine drops a file whose chunk or end it refused; one that it still holds (the file
      // could not be read) is released, so that the name is free for the next search. Not when
      // the runtime is ending: the instance goes with it.
      if (abi.stoppedBy === undefined) {
        try {
          inUse(use, () => abi.fileRelease(handle));
        } catch {
          // Already dropped by the engine, or the runtime is ending.
        }
      }
      throw error;
    } finally {
      if (buffer !== undefined) {
        const taken = buffer;
        try {
          inUse(use, () => abi.freeChunkBuffer(taken));
        } catch {
          // The runtime is ending: the instance goes with it.
        }
      }
    }
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
    // A second keeper on the Gateway's channel could answer for the other memory.
    if (!this.dispose()) throw new Error(ENDING);
    if (path === 'serial') {
      const instance = await instantiateSerial(await this.options.serial());
      this.loaded = { path, instance, runs: 0, retained: undefined, counts: undefined, ran: new Set() };
      return this.loaded;
    }
    const unsupported = threadedSupport();
    if (unsupported !== undefined) throw new Error(unsupported);
    const { use } = this.options;
    if (use === undefined) throw new Error('this page cannot share memory with the engine worker');
    const { module, artifact } = await this.options.threads();
    // The maximum that the module declares, or a smaller one where the browser cannot
    // reserve it; never a larger one (plan TD-7).
    const { memory, maximumPages } = createSharedMemory(artifact.memory);
    const { faultChannel, keeperChannel } = this.options;
    const threads = new ThreadHost({ module, memory, maximumPages, faultChannel, keeperChannel, use });
    try {
      const instance = await instantiateThreadedMain(module, memory, threads.spawn);
      this.loaded = { path, instance, threads, runs: 0, retained: undefined, counts: undefined, ran: new Set() };
      return this.loaded;
    } catch (error) {
      threads.terminate();
      throw error;
    }
  }
}

/** A runtime that the Gateway ends starts nothing and makes no call into its threaded instance. */
const ENDING = 'the runtime is ending';

/**
 * Calls into the instance without a wait between them, counted as a use of the threaded
 * instance (instance-use.ts): `use` is undefined for the serial instance.
 */
function inUse<T>(use: InstanceUse | undefined, body: () => T): T {
  if (use !== undefined && !use.enter()) throw new Error(ENDING);
  try {
    return body();
  } finally {
    use?.leave();
  }
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}

/** Bytes at each end of a counts file that its key includes (`edgeDigest`). */
export const COUNTS_EDGE_BYTES = 64 * 1024;

/**
 * The SHA-256 of the first and the last `COUNTS_EDGE_BYTES` of a file, in hex: with the name, size
 * and last-modified time, the key that tells a held counts file from another one without reading
 * it whole (the SHA-256 of the whole file is known only after the engine read it). The ends hold
 * the counts' metadata and parameters, and the end of their last table.
 */
export async function edgeDigest(blob: Blob): Promise<string> {
  const head = blob.slice(0, COUNTS_EDGE_BYTES);
  const tail = blob.slice(Math.max(COUNTS_EDGE_BYTES, blob.size - COUNTS_EDGE_BYTES));
  const bytes = await new Blob([head, tail]).arrayBuffer();
  const digest = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
  return Array.from(digest, (byte) => byte.toString(16).padStart(2, '0')).join('');
}
