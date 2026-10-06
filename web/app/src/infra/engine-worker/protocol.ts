// Messages between the main thread (gateway.ts), the Engine worker (engine-worker.ts) and
// its thread workers (thread-worker.ts). Every message of the Engine worker carries the
// runtime generation it belongs to; the main thread drops messages of an older one.
import type { RecordKey } from '../../domain/dataset';
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
  readonly output: MessagePort;
}

/** Test builds only (`__LOSAT_TEST_HOOKS__`): drives the run output writer of the worker. */
export type WorkerTestCommand =
  | { readonly type: 'test-open'; readonly id: number; readonly port: MessagePort }
  | { readonly type: 'test-write'; readonly id: number; readonly writer: number; readonly stream: OutputStream; readonly bytes: Uint8Array }
  | { readonly type: 'test-end'; readonly id: number; readonly writer: number }
  | { readonly type: 'test-post'; readonly id: number; readonly writer: number; readonly message: unknown };

export type EngineCommand = EngineInit | WorkerRun | WorkerTestCommand;

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

export interface ThreadPrepare {
  readonly type: 'prepare';
  readonly module: WebAssembly.Module;
  readonly memory: WebAssembly.Memory;
  /** Int32Array(1) on a SharedArrayBuffer: the slot state. */
  readonly slot: SharedArrayBuffer;
  readonly faultChannel: string;
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
