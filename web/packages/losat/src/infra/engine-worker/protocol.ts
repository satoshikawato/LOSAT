// Messages between the main thread (gateway.ts), the Engine worker (engine-worker.ts), its
// thread workers (thread-worker.ts) and the keeper of their shared memory (memory-keeper.ts).
// Every message of the Engine worker carries the runtime generation it belongs to; the main
// thread drops messages of an older one.
import type { RecordKey, ShownIds } from '../../domain/dataset.js';
import type { CountsUse } from '../../domain/counts.js';
import type { DirectoryFile, MessagePlace, RegistryFile, RunFile } from '../../domain/run-files.js';
import type { ToolFrame } from '../../domain/tool-frames.js';
import type { EnginePhase, ToolValidation } from '../../ports/engine.js';
import type { OutputStream } from '../../ports/run-output.js';
import type { ReactorArtifact } from '../reactor/artifact.js';
import type { ReactorAsset } from '../reactor/assets.js';

export type RuntimePath = 'serial' | 'threaded';

/**
 * A module for the worker: compiled once by the main thread and posted, or, where a
 * browser cannot post a compiled module, fetched and compiled by the worker.
 */
export type ModuleSource =
  | { readonly module: WebAssembly.Module; readonly artifact: ReactorArtifact }
  | { readonly asset: ReactorAsset };

export interface EngineInit {
  readonly type: 'init';
  readonly generation: number;
  readonly modules: { readonly serial: ModuleSource; readonly threads: ModuleSource };
  /** Engine build names for the run records. */
  readonly builds: { readonly serial: string; readonly threads: string };
  /**
   * The Web Lock that the worker holds for as long as it lives (worker-life.ts): the Gateway learns from its release
   * that the worker ended without an error event (line WQA's sweep, D2). The worker answers `started` at once, and
   * `life` once it holds the lock, before any preparation or search of the worker starts.
   */
  readonly lifeLock: string;
  /**
   * Whether the threaded instance is in use (instance-use.ts), shared with the Gateway. Absent
   * where the page cannot share memory with a worker: no threaded instance starts then.
   */
  readonly use?: SharedArrayBuffer;
}

export interface WorkerInput {
  readonly bytes: Uint8Array;
  readonly records: readonly RecordKey[];
  /**
   * Program, SHA-256 and dataset revisions of the input. A subject whose key equals the
   * one the worker holds is not registered again (plan §4.6 R1).
   */
  readonly key: string;
  /** Whether `bytes` was assembled from parts (`EngineInput.assembled`, line WQA's D17). */
  readonly assembled?: boolean;
}

/**
 * The WindowMasker counts file of a run (`-window_masker_db` or `-window_masker_taxid`). A File posted to the worker is
 * shared, not copied. A file whose `key` equals the one the instance holds is not given to the
 * engine again (plan §5, DESIGN_ww1 §5).
 */
export interface WorkerCounts {
  /** The value of `-window_masker_db` in the argv; for `-window_masker_taxid`, windowmasker/<taxid>/wmasker.obinary. */
  readonly name: string;
  readonly blob: Blob;
  /** Name, size and last-modified time of the file (counts-key); the runtime adds the SHA-256 of its ends (`edgeDigest`). */
  readonly key: string;
}

export interface WorkerRun {
  readonly type: 'run';
  readonly generation: number;
  readonly runId: string;
  readonly program: string;
  /** The run snapshot's argv, without -num_threads. */
  readonly argv: readonly string[];
  /** Threads to use; 1 runs the serial module. */
  readonly threads: number;
  /** None for a run of `-in_pssm` or `-in_msa`, which reads its query from `files` (EngineRunRequest.query). */
  readonly query?: WorkerInput;
  /** None for a run of `-db`, which searches the database in `files` (EngineRunRequest.subject). */
  readonly subject?: WorkerInput;
  /** The files of the run's working directory (EngineRunRequest.files); Files are shared, not copied. */
  readonly files?: readonly DirectoryFile[];
  /** The folder that the run's `BLASTDB` names (EngineRunRequest.blastdb). */
  readonly blastdb?: string;
  readonly windowMasker?: WorkerCounts;
  /**
   * The registry file of `-conffile` (EngineRunRequest.conffile): registered before the argv is
   * validated and released after the run; no later search of the instance reads it.
   */
  readonly conffile?: RegistryFile;
  readonly output: MessagePort;
}

/**
 * A tool run (`Engine.runTool`, ABI v2 `losat_web2_tool`, version 13): the words of the command line, `argv[0]` the tool's
 * name, and the files of the run's working directory (Files are shared, not copied). It runs on the serial module. Its
 * standard output, standard error and files come back as `tool-chunk` and `tool-file` events as the tool writes them, and its
 * exit status in `tool-done`, which comes last; an error comes as `done` with `ok: false`.
 */
export interface WorkerTool {
  readonly type: 'tool';
  readonly generation: number;
  readonly runId: string;
  readonly argv: readonly string[];
  readonly files?: readonly DirectoryFile[];
  /**
   * The counter of the bytes that the worker has posted and the page has not taken (tool-credit.ts): the worker waits in the
   * engine's call while it is above the limit. Absent where the page cannot share memory with the worker: nothing waits then.
   */
  readonly credit?: SharedArrayBuffer;
}

/**
 * The check of a tool's words (`ToolEngine.validateTool`, ABI v2 `losat_web2_validate` with a tool's name as `argv[0]`): it
 * runs on the instance that the worker has loaded (the serial one when none is), in the worker's turn like a run (behind a run in progress), with the files of
 * the run's directory opened there as for a run, because the check reads the run's registry files. The answer is
 * `tool-validated` with the same `id`.
 */
export interface WorkerToolValidate {
  readonly type: 'tool-validate';
  readonly generation: number;
  readonly id: string;
  readonly argv: readonly string[];
  readonly files?: readonly DirectoryFile[];
}

/**
 * Starts what a search of `threads` threads would start (the instance and, for more than one
 * thread, its thread workers) before the search comes (EngineGateway.prepare). A failure is
 * not reported: the search that follows starts them again and reports it.
 */
export interface WorkerPrepare {
  readonly type: 'prepare';
  readonly generation: number;
  readonly threads: number;
  /** The argv of the coming search (EnginePrepareHint.argv), for the warm-up. */
  readonly argv?: readonly string[];
}

/**
 * Test builds only (`__LOSAT_TEST_HOOKS__`): the BroadcastChannel on which a page arms the
 * Engine worker's next search to stop at its first output ({type: 'arm'}, answered by
 * {type: 'armed'}).
 */
export const STOP_AT_OUTPUT_CHANNEL = 'losat-web:test:stop-at-output';

/**
 * Test builds only (`__LOSAT_TEST_HOOKS__`): the BroadcastChannel on which a page arms the
 * Engine worker's next use of the threaded instance (a search or a warm-up) to hold at its
 * start until the Gateway ends the runtime ({type: 'arm'}, answered by {type: 'armed'}; the
 * worker posts {type: 'holding'} when it holds, and {type: 'not-held', reason} after a
 * preparation with an argv that did not hold).
 */
export const HOLD_USE_CHANNEL = 'losat-web:test:hold-use';

/** Test builds only (`__LOSAT_TEST_HOOKS__`): drives the run output writer of the worker. */
export type WorkerTestCommand =
  | { readonly type: 'test-open'; readonly id: number; readonly port: MessagePort }
  | { readonly type: 'test-write'; readonly id: number; readonly writer: number; readonly stream: OutputStream; readonly bytes: Uint8Array }
  | { readonly type: 'test-end'; readonly id: number; readonly writer: number }
  | { readonly type: 'test-post'; readonly id: number; readonly writer: number; readonly message: unknown };

export type EngineCommand = EngineInit | WorkerRun | WorkerTool | WorkerToolValidate | WorkerPrepare | WorkerTestCommand;

export interface MemoryRecord {
  /** Linear memory of the instance before and after the search, in bytes. */
  readonly linearBytesBefore: number;
  readonly linearBytesAfter: number;
  /** Searches this instance has run, including this one. */
  readonly instanceRuns: number;
}

export interface WorkerRunResult {
  readonly path: RuntimePath;
  readonly threads: number;
  readonly engineBuild: string;
  readonly fallbackReason?: string;
  readonly memory: MemoryRecord;
  /** Whether the search used the subject that the instance held (plan §4.6 R1). */
  readonly subjectRetained: boolean;
  /** Where the run's messages go, as the engine said (stream 10; `RuntimeInfo.messages`). */
  readonly messages: MessagePlace;
  /** The IDs that the outputs show for the records, where the engine read them from the deflines (-parse_deflines). */
  readonly shownIds?: ShownIds;
  /** What the run did with its WindowMasker counts file, when it had one. */
  readonly windowMasker?: CountsUse;
  /** The bytes of the counts file that the instance holds after the run (the renewal's high-water mark allows for them); absent when it holds none. */
  readonly countsBytes?: number;
  /** The files that the run made besides its outputs (stream 8), when it made any. */
  readonly runFiles?: readonly RunFile[];
  /** The text that the run printed in place of a search (stream 11; `RuntimeInfo.printed`). */
  readonly printed?: Uint8Array;
  /** Records read at `register` for an assembled input when they differ from the per-file table (`RuntimeInfo.readRecords`). */
  readonly readRecords?: {
    readonly query?: readonly RecordKey[];
    readonly subject?: readonly RecordKey[];
  };
}

/** The end of a tool run that the engine ran to its exit status (`tool-done`). */
export interface WorkerToolResult {
  /** The command line's exit status. */
  readonly status: number;
  /** Linear memory of the instance around the run, as for a search: the renewal policy reads it. */
  readonly memory: MemoryRecord;
  /** The bytes of the counts file that the instance holds after the run; absent when it holds none. */
  readonly countsBytes?: number;
}

export interface WorkerError {
  readonly name: string;
  readonly message: string;
  /** InputMismatchError: the role and the difference. */
  readonly role?: 'query' | 'subject';
  readonly detail?: string;
  /** The instance stopped; the worker must not be used again. */
  readonly stopped?: boolean;
  /** The files that the failed run made besides its outputs (stream 8; `runFilesOf`). */
  readonly runFiles?: readonly RunFile[];
  /** The failed run ended its output with `fail` (`runOutputFailedOf`). */
  readonly outputFailed?: true;
}

export type EngineEvent =
  /**
   * The worker started: its script loaded and took `init`. An error before this is a worker that could not start (line
   * WQA's sweep, D1); an error after it, a worker that stopped.
   */
  | { readonly type: 'started'; readonly generation: number }
  /**
   * Whether the worker holds the lock `lifeLock` (EngineInit), posted before any preparation or search of the worker
   * starts: the Gateway waits for that lock from then on, also during the worker's first search (line WQA's sweep, D2).
   */
  | { readonly type: 'life'; readonly generation: number; readonly held: boolean }
  | WorkerLink
  /** A fault of the ThreadHost itself (a thread that did not start in time), which ends the runtime as a thread's does. */
  | ({ readonly generation: number } & ThreadFault)
  | { readonly type: 'phase'; readonly generation: number; readonly runId: string; readonly phase: EnginePhase }
  | { readonly type: 'done'; readonly generation: number; readonly runId: string; readonly ok: true; readonly result: WorkerRunResult }
  /**
   * A chunk of a tool's standard output or standard error, posted as the tool writes it, with its buffer transferred. The
   * events of a run reach the main thread in the order that the worker posted them: all the chunks and frames of a tool run
   * come before its `tool-done`.
   */
  | {
      readonly type: 'tool-chunk';
      readonly generation: number;
      readonly runId: string;
      readonly stream: 'stdout' | 'stderr';
      readonly bytes: Uint8Array;
    }
  /** A frame of a file that a tool writes (stream 12), posted as it arrives, with its bytes transferred. */
  | { readonly type: 'tool-file'; readonly generation: number; readonly runId: string; readonly frame: ToolFrame }
  /**
   * The answer to a `tool-validate`: the verdict, or the error that stopped the check (a stopped instance sets `stopped`: the
   * Gateway ends the worker and posts the checks that are still pending to the next one).
   */
  | { readonly type: 'tool-validated'; readonly generation: number; readonly id: string; readonly verdict: ToolValidation }
  | { readonly type: 'tool-validate-failed'; readonly generation: number; readonly id: string; readonly error: WorkerError }
  | { readonly type: 'tool-done'; readonly generation: number; readonly runId: string; readonly result: WorkerToolResult }
  | { readonly type: 'done'; readonly generation: number; readonly runId: string; readonly ok: false; readonly error: WorkerError }
  | { readonly type: 'test-reply'; readonly id: number; readonly ok: boolean; readonly value?: number; readonly error?: WorkerError };

// --- thread workers -------------------------------------------------------------------

/** Slot states, in the slot's shared Int32Array (index 0). */
export const SLOT_PREPARING = 0;
export const SLOT_READY = 1;
export const SLOT_RUNNING = 2;
export const SLOT_FAILED = -1;
/** The ThreadHost ended the slot's worker while its thread ran (ThreadHost.terminate). */
export const SLOT_ENDED = -2;

export interface ThreadPrepare {
  readonly type: 'prepare';
  readonly module: WebAssembly.Module;
  readonly memory: WebAssembly.Memory;
  /** Int32Array(1) on a SharedArrayBuffer: the slot state. */
  readonly slot: SharedArrayBuffer;
  /** The worker's link to the main thread (link.ts), on which it reports a thread that stopped (ThreadFault). */
  readonly faults: MessagePort;
  /** The use count of the instance (instance-use.ts): a thread counts itself off when it returns. */
  readonly use: SharedArrayBuffer;
  /** The table of the descriptors of the run's working directory that the instance's threads share (run-directory.ts). */
  readonly directory: SharedArrayBuffer;
}

/** Start states, in the start message's shared Int32Array (index 0). */
export const START_PENDING = 0;
export const START_TAKEN = 1;
/** The thread worker cannot start the thread, or the host stopped waiting for it. */
export const START_ABANDONED = -1;

export interface ThreadStart {
  readonly type: 'start';
  readonly tid: number;
  readonly startArg: number;
  /** Int32Array(1) on a SharedArrayBuffer: START_PENDING, then START_TAKEN or START_ABANDONED. */
  readonly started: SharedArrayBuffer;
  /** The files of the run's working directory (run-directory.ts), which the thread reads as the run's other threads do. */
  readonly files: readonly DirectoryFile[];
}

export type ThreadCommand = ThreadPrepare | ThreadStart;

/** Posted on a thread worker's link when its thread stops (for example, a trap). */
export interface ThreadFault {
  readonly type: 'fault';
  readonly tid: number;
  readonly message: string;
}

// --- the memory keeper ------------------------------------------------------------------

/**
 * From the ThreadHost to the keeper of its shared memory (memory-keeper.ts): the memory and
 * the keeper's link to the main thread (link.ts), on which the Gateway asks it to grow the memory.
 */
export interface KeeperCommand {
  readonly type: 'keep';
  readonly memory: WebAssembly.Memory;
  readonly maximumPages: number;
  readonly port: MessagePort;
}

/**
 * From the keeper to the ThreadHost: it holds the memory, and its link carries messages (ThreadHost.prepare waits for it).
 */
export interface KeeperReady {
  readonly type: 'ready';
}

/** From the Gateway to the keeper on its link: grow the memory to its maximum. */
export interface KeeperGrow {
  readonly type: 'grow';
}

/** The keeper's answer on its link, after the grow. */
export interface KeeperGrown {
  readonly type: 'grown';
}

// --- links ------------------------------------------------------------------------------

/** Whose link a port is: a thread worker's, which reports faults, or the keeper's, which the Gateway asks to grow the memory. */
export type LinkRole = 'thread' | 'keeper';

/**
 * From the Engine worker to the main thread: the main thread's end of the link of a thread worker or the keeper (link.ts),
 * posted before the worker or the keeper gets the other end.
 */
export interface WorkerLink {
  readonly type: 'link';
  readonly generation: number;
  readonly role: LinkRole;
  readonly port: MessagePort;
}

/** The first message of a link, from the worker's end, and the main thread's answer to it (link.ts). */
export interface LinkHello {
  readonly type: 'hello';
}

/** A message on a link. */
export type LinkMessage = LinkHello | ThreadFault | KeeperGrow | KeeperGrown;
