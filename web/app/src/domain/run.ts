// Run records (plan §5.2). A RunSnapshot is fixed when a job is queued and never
// changes; what actually happened is appended to the RunRecord.
import type { RecordKey } from './dataset';
import type { ProgramId } from './programs';

export type RunStatus =
  | 'queued'
  | 'preparing'
  | 'running'
  | 'finalizing'
  | 'completed'
  | 'cancelled'
  | 'failed';

export const TERMINAL_STATUSES: ReadonlySet<RunStatus> = new Set(['completed', 'cancelled', 'failed']);

export interface InputSnapshot {
  /** Name passed as -query / -subject (plan §5.3). */
  readonly name: string;
  /** Exact FASTA bytes given to the engine. */
  readonly bytes: Uint8Array;
  /** Lower-case hex SHA-256 of `bytes`. */
  readonly sha256: string;
  /** Dataset revisions whose included records make up `bytes`, in order. */
  readonly revisionIds: readonly string[];
  /** ID and length of each record in `bytes`; the engine checks them at `register`. */
  readonly records: readonly RecordKey[];
}

export interface RunSnapshot {
  readonly runId: string;
  /** 1-based position in this working session, used for display and file names. */
  readonly number: number;
  readonly program: ProgramId;
  readonly argv: readonly string[];
  readonly query: InputSnapshot;
  readonly subject: InputSnapshot;
  readonly requestedThreads: number | 'auto';
  readonly queuedAt: number;
}

export interface RunRecord {
  readonly runtimePath?: 'threaded' | 'serial' | 'fake';
  readonly threads?: number;
  readonly fallbackReason?: string;
  readonly engineBuild?: string;
  /** The engine runtime (Engine worker) that ran the search; it changes after a cancel. */
  readonly runtimeGeneration?: number;
  /** Linear memory of the engine instance around the search (plan §5.5). */
  readonly memory?: {
    readonly linearBytesBefore: number;
    readonly linearBytesAfter: number;
    readonly instanceRuns: number;
  };
  /** The engine searched the subject that it held from an earlier search (R1). */
  readonly subjectRetained?: boolean;
  readonly startedAt?: number;
  /** When each engine phase started. */
  readonly phaseTimes?: { readonly preparing?: number; readonly running?: number; readonly finalizing?: number };
  readonly endedAt?: number;
  readonly error?: string;
}

export function isTerminal(status: RunStatus): boolean {
  return TERMINAL_STATUSES.has(status);
}
