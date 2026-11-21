// The engine inside the Engine worker (plan §3.1, §4.6, §4.7): one warm reactor instance,
// serial or threaded, the subject it holds (R1), and one search at a time.
//
// A search: choose the module (threaded for more than one thread, unless threaded
// searches are not possible here, then serial with one thread and the reason), append
// -num_threads and validate that argv (ABI v2 §7), register the query and, unless the
// instance holds it, the subject, compare the records that the engine read with the
// record tables (plan §5.4; a difference stops the run before any output), run, and send
// the output streams to the Data worker's port, with the engine's signal of each output that
// failed (stream 9: the Data worker drops what it received of that output). A BLASTN search
// with -window_masker_db also gives the engine the counts file, which the instance holds for
// later BLASTN searches. A search with -conffile gives the engine the registry file for that
// search only, and the files that a run makes besides its outputs (stream 8: the file of
// -logfile) come back with its result, or with its error when it fails (line RG E1w-3b); they
// never go to the Data worker. Where the run's messages go (stream 10: the diagnostics, the file of
// -logfile, or nowhere) comes back with its result too (line WQA's sweep, D7, 2026-11-20), and so does
// the text that a run prints in place of a search (stream 11: the command line's help or version), which
// registers no input (line WQA's D14, 2026-11-20). The warnings of a run (stream 3) go to the Data worker;
// a run that the engine failed ends its output with `fail`, so that the Data worker keeps those of a failed
// run (line WQA's sweep, D3, 2026-11-19).
import { recordMismatch, type RecordKey } from '../../domain/dataset.js';
import { readStoredZip, type MessagePlace, type OutputPlace, type RunFile } from '../../domain/run-files.js';
import { InputMismatchError, withRunFiles, withRunOutputFailed, type EnginePhase } from '../../ports/engine.js';
import { DIAGNOSTICS_STREAM, isRunOutput, OUTPUT_STREAMS, type OutputStream } from '../../ports/run-output.js';
import type { CountsUse } from '../../domain/counts.js';
import type { RunOutput } from '../../domain/output-format.js';
import {
  discardSignal,
  DISCARD_STREAM,
  EngineCallError,
  FILE_REGISTRY,
  FILE_WINDOW_MASKER_COUNTS,
  loadFile,
  MESSAGES_STREAM,
  messagesSignal,
  PRINTED_STREAM,
  ROLE_QUERY,
  ROLE_SUBJECT,
  RUN_FILES_STREAM,
  type FileInfo,
  type RegisteredInputs,
} from '../reactor/abi.js';
import type { ReactorArtifact } from '../reactor/artifact.js';
import { instantiateSerial, instantiateThreadedMain, type ReactorInstance } from '../reactor/instance.js';
import { createSharedMemory } from '../reactor/shared-memory.js';
import { RunOutputWriter } from '../run-output/writer.js';
import type { InstanceUse } from './instance-use.js';
import type { MainLinks } from './link.js';
import type { RuntimePath, WorkerCounts, WorkerRun, WorkerRunResult } from './protocol.js';
import { ThreadHost } from './thread-host.js';
import { warmUp } from './warm-up.js';

export interface RuntimeOptions {
  readonly serial: () => Promise<WebAssembly.Module>;
  readonly threads: () => Promise<{ readonly module: WebAssembly.Module; readonly artifact: ReactorArtifact }>;
  readonly builds: { readonly serial: string; readonly threads: string };
  /** The links of the thread workers and the keeper of the threaded instance's shared memory to the main thread (link.ts). */
  readonly main: MainLinks;
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
  /** Whether the engine read the subject with NCBI's data loaders on (a run with them on takes no subject read with them off). */
  readonly dataLoaders: boolean;
  readonly handle: number;
  readonly records: readonly RecordKey[];
}

/** The WindowMasker counts file that the instance holds: at most one (DESIGN_ww1 §5). */
interface HeldCounts {
  readonly key: string;
  readonly handle: number;
  readonly info: FileInfo;
}

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
    let windowMasker: CountsUse | undefined;
    let registry: number | undefined;
    try {
      windowMasker = await this.counts(loaded, request, use);
      // The registry file of -conffile, which the start of `validate` and `run` reads (ABI v2 §4),
      // for this search only (`search` releases it): the start of every later call would read a file
      // named `.ncbirc` or `<program>.ini` as the CLI reads one in its directory. After the counts
      // file, whose name it cannot take while the instance holds that file.
      if (request.conffile !== undefined) {
        const { name, file } = request.conffile;
        registry = (await loadFile(loaded.instance.abi, FILE_REGISTRY, name, file, (body) => inUse(use, body))).handle;
      }
    } catch (error) {
      // A trap while a file was given to the engine: the next search loads the module again
      // (a stopped threaded instance is ended by the Gateway, as in `search`).
      if (loaded.instance.abi.stoppedBy !== undefined && loaded.path === 'serial') this.dispose();
      throw error;
    }
    if (use !== undefined && !use.enter()) throw new Error(ENDING);
    try {
      if (use !== undefined) this.options.onUse?.();
      return this.search(loaded, request, fallbackReason, phase, windowMasker, registry);
    } finally {
      use?.leave();
    }
  }

  /** The search of `run` on the instance that `run` chose, with the handle of the registry file that `run` gave the engine. */
  private search(
    loaded: Loaded,
    request: WorkerRun,
    fallbackReason: string | undefined,
    phase: (phase: EnginePhase) => void,
    windowMasker: CountsUse | undefined,
    registryHandle: number | undefined,
  ): WorkerRunResult {
    const threads = loaded.path === 'threaded' ? request.threads : 1;
    const argv = [...request.argv, '-num_threads', String(threads)];
    const abi = loaded.instance.abi;
    // What the engine holds for this search alone: the registry file of -conffile, the query, and the
    // subject of inputs read together for this argv (`registerInputs`). Each is released once, before
    // `end` or when the search fails, unless the instance stopped.
    let registry = registryHandle;
    let query: number | undefined;
    let ownSubject: number | undefined;
    const release = () => {
      if (abi.stoppedBy !== undefined) return;
      if (query !== undefined) {
        const handle = query;
        query = undefined;
        abi.release(handle);
      }
      if (ownSubject !== undefined) {
        const handle = ownSubject;
        ownSubject = undefined;
        abi.release(handle);
      }
      if (registry !== undefined) {
        const handle = registry;
        registry = undefined;
        abi.fileRelease(handle);
      }
    };
    try {
      const validation = abi.validate(argv);
      if (!validation.ok) throw new EngineCallError(validation.message);
      // The inputs are read with the data loaders of the argv's registry, as the CLI reads them (line WQA's D16, 2026-11-20).
      const { dataLoaders } = validation;
      const started = abi.start(argv);
      if (!started.ok) throw new EngineCallError(started.message);

      // A run that prints the command line's help or version in place of a search reads no input, as the command line
      // reads no file for it (ABI v2 §4): nothing is registered, and the run takes the handles 0 (line WQA's D14).
      const prints = started.start.prints !== null;
      // -parse_deflines (line WQA's sweep, D7, 2026-11-19): the CLI reads the subjects and then the
      // query of this argv in one scope of IDs, so the engine reads both inputs together for it, and
      // the two handles are this search's alone. The retained subject, read without the IDs, is
      // released first: the instance does not hold two copies of the subject.
      let together: RegisteredInputs | undefined;
      if (!prints && started.start.registerInputs) {
        this.releaseRetained(loaded);
        together = abi.registerInputs(argv, request.query.bytes, request.subject.bytes);
        ownSubject = together.subject.handle;
      }
      let handles: readonly [number, number] = [0, 0];
      let subjectRetained = false;
      if (!prints) {
        const registered = together?.query ?? abi.register(request.program, ROLE_QUERY, dataLoaders, request.query.bytes);
        query = registered.handle;
        const queryMismatch = recordMismatch(request.query.records, registered.records);
        if (queryMismatch !== undefined) throw new InputMismatchError('query', queryMismatch);
        subjectRetained = together === undefined && this.retains(loaded, request, dataLoaders);
        const subject = together?.subject ?? this.subject(loaded, request, dataLoaders);
        const subjectMismatch = recordMismatch(request.subject.records, subject.records);
        if (subjectMismatch !== undefined) throw new InputMismatchError('subject', subjectMismatch);
        handles = [registered.handle, subject.handle];
      }
      const shownIds = together?.shownIds;

      const writer = new RunOutputWriter(request.output);
      const linearBytesBefore = abi.memoryBytes();
      // The engine's standard error of this search only, for a stop message.
      loaded.instance.wasi.clear();
      phase('running');
      // The run's files (stream 8) stay here, apart from the outputs, and so does where its messages go (stream 10).
      // The text that the run prints in place of a search (stream 11) stays here too.
      const archive: Uint8Array[] = [];
      const printed: Uint8Array[] = [];
      let messages: MessagePlace | undefined;
      try {
        abi.run(argv, ...handles, (stream, bytes) => {
          this.options.beforeOutput?.();
          if (stream === RUN_FILES_STREAM) archive.push(bytes.slice());
          else if (stream === PRINTED_STREAM) printed.push(bytes.slice());
          else if (stream === MESSAGES_STREAM) messages = messagesPlace(bytes);
          else if (stream === DISCARD_STREAM) writer.discard(...discardedOutput(bytes));
          else if (OUTPUT_STREAMS.includes(stream as OutputStream)) writer.write(stream as OutputStream, bytes);
        });
      } catch (error) {
        // A run that the engine failed gives the files that it made (TBLASTX's file of -logfile holds
        // its error), and ends its output with `fail`, so that the Data worker keeps what it wrote to
        // the CLI's standard error before its error (stream 3; line WQA's sweep, D3, 2026-11-19); a
        // stop does neither (the chunks before a trap are discarded, ABI v2 §4).
        if (abi.stoppedBy !== undefined) throw error;
        writer.fail();
        throw failedWithFiles(error, archive);
      }
      const runFiles = readRunFiles(archive);
      if (messages === undefined) throw new Error(`the engine did not say where the messages of the run go (stream ${MESSAGES_STREAM})`);
      const linearBytesAfter = abi.memoryBytes();
      loaded.runs++;
      loaded.ran.add(request.argv.join('\0'));
      // The last call into the instance comes before `end`: a failure after it must not
      // follow an output that the data layer may already commit.
      release();
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
        messages,
        ...(shownIds === undefined ? {} : { shownIds }),
        ...(windowMasker === undefined ? {} : { windowMasker }),
        ...(loaded.counts === undefined ? {} : { countsBytes: loaded.counts.info.size }),
        ...(runFiles.length === 0 ? {} : { runFiles }),
        ...(prints ? { printed: joined(printed) } : {}),
      };
    } finally {
      if (abi.stoppedBy === undefined) release();
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
  private async counts(loaded: Loaded, request: WorkerRun, use: InstanceUse | undefined): Promise<CountsUse | undefined> {
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

  /** Gives the file to the engine in chunks (`loadFile`), ends it, and holds it. */
  private async loadCounts(loaded: Loaded, wanted: WorkerCounts, key: string, use: InstanceUse | undefined): Promise<FileInfo> {
    const info = await loadFile(loaded.instance.abi, FILE_WINDOW_MASKER_COUNTS, wanted.name, wanted.blob, (body) => inUse(use, body));
    loaded.counts = { key, handle: info.handle, info };
    return info;
  }

  /**
   * Whether the instance retains the run's subject: the same key, read with the data loaders on, or
   * off as the run has them (a subject read with them on is the one that a reader without them reads).
   */
  private retains(loaded: Loaded, request: WorkerRun, dataLoaders: boolean): boolean {
    const retained = loaded.retained;
    return retained !== undefined && retained.key === request.subject.key && (retained.dataLoaders || !dataLoaders);
  }

  /** The subject handle: the retained one if it is the run's (`retains`), else a new registration. */
  private subject(loaded: Loaded, request: WorkerRun, dataLoaders: boolean): Retained {
    const abi = loaded.instance.abi;
    if (this.retains(loaded, request, dataLoaders)) return loaded.retained!;
    this.releaseRetained(loaded);
    const registered = abi.register(request.program, ROLE_SUBJECT, dataLoaders, request.subject.bytes);
    loaded.retained = { key: request.subject.key, dataLoaders, handle: registered.handle, records: registered.records };
    return loaded.retained;
  }

  /** Releases the subject that the instance retains, if it retains one. */
  private releaseRetained(loaded: Loaded): void {
    if (loaded.retained === undefined) return;
    const handle = loaded.retained.handle;
    loaded.retained = undefined;
    loaded.instance.abi.release(handle);
  }

  private async load(path: RuntimePath): Promise<Loaded> {
    const current = this.loaded;
    if (current?.path === path && current.instance.abi.stoppedBy === undefined) return current;
    // The Gateway asks the keeper whose link came last (gateway.ts): the keeper of the other memory ends first.
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
    const threads = new ThreadHost({ module, memory, maximumPages, main: this.options.main, use });
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

/**
 * The output that the engine's signal of stream 9 names (ABI v2 §5), an output stream of the run
 * that failed, and where its error is: stream 3 (the diagnostics), stream 8 (the run's files: the
 * file of -logfile), nowhere, or the error that ended the run (an output that a failed run does not keep, version 7). Another
 * number breaks the ABI.
 */
function discardedOutput(bytes: Uint8Array): [RunOutput, OutputPlace] {
  const { stream, errorStream, ended } = discardSignal(bytes);
  if (!isRunOutput(stream)) throw new Error(`the engine discarded stream ${stream}, which is not an output of a run`);
  const error: OutputPlace | undefined = ended ? 'run' : errorStream === undefined ? 'nowhere' : PLACES.get(errorStream);
  if (error === undefined) throw new Error(`the engine put the error of output stream ${stream} on stream ${errorStream}, which holds no error`);
  return [stream, error];
}

/**
 * Where the run's messages go, as the engine's signal of stream 10 says (ABI v2 §5): stream 3 (the diagnostics), stream 8
 * (the run's files: the file of -logfile) or nowhere. Another number breaks the ABI (line WQA's sweep, D7, 2026-11-20).
 */
function messagesPlace(bytes: Uint8Array): MessagePlace {
  const stream = messagesSignal(bytes);
  const place = stream === undefined ? 'nowhere' : PLACES.get(stream);
  if (place === undefined) throw new Error(`the engine put the messages of the run on stream ${stream}, which holds no messages`);
  return place;
}

/** The place of messages that a signal of stream 9 or 10 gives as a stream (ABI v2 §5). */
const PLACES: ReadonlyMap<number, MessagePlace> = new Map([
  [DIAGNOSTICS_STREAM, 'diagnostics'],
  [RUN_FILES_STREAM, 'logfile'],
]);

/** The chunks of a stream as one array. */
function joined(chunks: readonly Uint8Array[]): Uint8Array {
  const bytes = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
  let offset = 0;
  for (const chunk of chunks) {
    bytes.set(chunk, offset);
    offset += chunk.length;
  }
  return bytes;
}

/** The files of a run's archive (stream 8, in chunks); none when the run sent none. */
function readRunFiles(chunks: readonly Uint8Array[]): RunFile[] {
  if (chunks.length === 0) return [];
  const archive = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
  let offset = 0;
  for (const chunk of chunks) {
    archive.set(chunk, offset);
    offset += chunk.length;
  }
  return readStoredZip(archive);
}

/**
 * The error of a run that the engine failed after it ended the run's output with `fail`
 * (`withRunOutputFailed`), with the files that the run made (`withRunFiles`). An archive that cannot
 * be read is told after the run's error, which stays the run's.
 */
function failedWithFiles(error: unknown, chunks: readonly Uint8Array[]): unknown {
  if (!(error instanceof Error)) return error;
  const failed = withRunOutputFailed(error);
  let runFiles: RunFile[];
  try {
    runFiles = readRunFiles(chunks);
  } catch (unreadable) {
    failed.message = `${failed.message} (${messageOf(unreadable)})`;
    return failed;
  }
  return runFiles.length === 0 ? failed : withRunFiles(failed, runFiles);
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
