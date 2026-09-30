// Run records (plan §5.2). A RunSnapshot is fixed when a job is queued and never
// changes; what actually happened is appended to the RunRecord.
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
  readonly startedAt?: number;
  readonly endedAt?: number;
  readonly error?: string;
}

export function isTerminal(status: RunStatus): boolean {
  return TERMINAL_STATUSES.has(status);
}
