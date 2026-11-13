// Messages between the main thread (gateway.ts), the Engine worker (engine-worker.ts), its
// thread workers (thread-worker.ts) and the keeper of their shared memory (memory-keeper.ts).
// Every message of the Engine worker carries the runtime generation it belongs to; the main
// thread drops messages of an older one.
import type { RecordKey } from '../../domain/dataset';
import type { CountsUse } from '../../domain/counts';
import type { EnginePhase } from '../../ports/engine';
import type { OutputStream } from '../../ports/run-output';
import type { ReactorArtifact } from '../reactor/artifact';
import type { ReactorAsset } from '../reactor/assets';

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
  /** BroadcastChannel on which thread workers report a thread that stopped. */
  readonly faultChannel: string;
  /** BroadcastChannel on which the Gateway asks the keeper of the shared memory to grow it (memory-keeper.ts). */
  readonly keeperChannel: string;
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
  readonly query: WorkerInput;
  readonly subject: WorkerInput;
  readonly windowMasker?: WorkerCounts;
  readonly output: MessagePort;
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

export type EngineCommand = EngineInit | WorkerRun | WorkerPrepare | WorkerTestCommand;

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
  /** What the run did with its WindowMasker counts file, when it had one. */
  readonly windowMasker?: CountsUse;
  /** The bytes of the counts file that the instance holds after the run (the renewal's high-water mark allows for them); absent when it holds none. */
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
}

export type EngineEvent =
  | { readonly type: 'phase'; readonly generation: number; readonly runId: string; readonly phase: EnginePhase }
  | { readonly type: 'done'; readonly generation: number; readonly runId: string; readonly ok: true; readonly result: WorkerRunResult }
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
  readonly faultChannel: string;
  /** The use count of the instance (instance-use.ts): a thread counts itself off when it returns. */
  readonly use: SharedArrayBuffer;
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
}

export type ThreadCommand = ThreadPrepare | ThreadStart;

/** Posted on the fault BroadcastChannel when a thread stops (for example, a trap). */
export interface ThreadFault {
  readonly type: 'fault';
  readonly tid: number;
  readonly message: string;
}

// --- the memory keeper ------------------------------------------------------------------

/**
 * From the ThreadHost to the keeper of its shared memory (memory-keeper.ts): the memory and
 * the BroadcastChannel on which the Gateway asks the keeper to grow it.
 */
export interface KeeperCommand {
  readonly type: 'keep';
  readonly memory: WebAssembly.Memory;
  readonly maximumPages: number;
  readonly channel: string;
}

/** From the keeper to the ThreadHost: it holds the memory and listens (ThreadHost.prepare waits for it). */
export interface KeeperReady {
  readonly type: 'ready';
}

/** From the Gateway to the keeper on the channel: grow the memory to its maximum. */
export interface KeeperGrow {
  readonly type: 'grow';
}

/** The keeper's answer on the channel, after the grow. */
export interface KeeperGrown {
  readonly type: 'grown';
}
