// Messages between the harness page and its workers.
import type { CaseResult } from '../../contract/contract';

export type BlockStoreWorkerMessage =
  | { readonly type: 'exhaust' | 'restore'; readonly id: number }
  | { readonly type: 'progress'; readonly result: CaseResult }
  | { readonly type: 'results'; readonly results: CaseResult[] };

export type BlockStorePageMessage =
  | { readonly type: 'start'; readonly backend: 'opfs' | 'memory' }
  | { readonly type: 'done'; readonly id: number };

export type EngineDoubleCommand =
  | { readonly id: number; readonly type: 'open'; readonly port: MessagePort }
  | { readonly id: number; readonly type: 'write'; readonly writer: number; readonly stream: number; readonly bytes: Uint8Array }
  | { readonly id: number; readonly type: 'end'; readonly writer: number }
  | { readonly id: number; readonly type: 'post'; readonly writer: number; readonly message: unknown };

export type EngineDoubleReply =
  | { readonly id: number; readonly ok: true; readonly value?: number }
  | { readonly id: number; readonly ok: false; readonly name: string; readonly message: string };

/** One search of the high-memory worker: the engine's argv without -num_threads, and its inputs. */
export interface HighMemorySearch {
  readonly program: string;
  readonly argv: readonly string[];
  readonly query: Uint8Array;
  readonly subject: Uint8Array;
}

export interface HighMemoryRequest {
  readonly threads: number;
  readonly faultChannel: string;
  readonly cases: readonly HighMemorySearch[];
}

/** SHA-256 and length of every stream that a search emitted (0, 6, 7, 1 and 3). */
export interface HighMemoryCase {
  readonly program: string;
  readonly sha256: Readonly<Record<number, string>>;
  readonly bytes: Readonly<Record<number, number>>;
}

export interface HighMemoryMode {
  /** The maximum of the shared memory that the browser provided, in pages. */
  readonly maximumPages: number;
  /** Linear memory when the instance started and after the searches, in bytes. */
  readonly initialBytes: number;
  readonly linearBytes: number;
  /** The first allocation above 2 GiB after the memory below was used up; 0 in the normal mode. */
  readonly firstHighPointer: number;
  readonly outputs: readonly HighMemoryCase[];
}

export interface HighMemoryResult {
  /** The maximum that the threaded module declares, in pages. */
  readonly declaredPages: number;
  readonly normal: HighMemoryMode;
  readonly high: HighMemoryMode;
}
