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
