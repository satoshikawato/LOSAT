// Run records (plan §5.2). A RunSnapshot is fixed when a job is queued and never
// changes; what actually happened is appended to the RunRecord.
import type { CountsOrigin } from './counts';
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
  /**
   * Dataset revisions whose included records made up `bytes`, in order: ids only (the engine's
   * retained-subject key). The data layer may have dropped them since (SURVEY #15); `bytes` and
   * `records` are the run's input.
   */
  readonly revisionIds: readonly string[];
  /** ID and length of each record in `bytes`; the engine checks them at `register`. */
  readonly records: readonly RecordKey[];
}

/**
 * The WindowMasker unit counts file of a run (`-window_masker_db`, or `-window_masker_taxid`): the file
 * that the run's engine reads under the name below. It lives in the working session only, as the
 * run does (the file of a taxonomy ID is kept by the browser, and the run reads it from there).
 */
export interface CountsSnapshot {
  /** The value of `-window_masker_db` in the argv; for `-window_masker_taxid`, the name that the engine looks for (windowmasker/<taxid>/wmasker.obinary). */
  readonly name: string;
  readonly file: Blob;
  readonly size: number;
  readonly lastModified: number;
  /** For `-window_masker_taxid`: where the file was downloaded from. */
  readonly origin?: CountsOrigin;
}

/** What a run did with its counts file (RunRecord.windowMasker). */
export interface CountsRecord {
  readonly name: string;
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes that the engine read. */
  readonly sha256: string;
  /** The engine used the counts that it held from an earlier run. */
  readonly retained: boolean;
  /** Milliseconds spent giving the file to the engine; 0 when it was retained. */
  readonly loadMs: number;
  /** For `-window_masker_taxid`: the taxonomy ID, the URL and the Last-Modified of the file that was read. */
  readonly origin?: CountsOrigin;
}

export interface RunSnapshot {
  readonly runId: string;
  /** 1-based position in this working session, used for display and file names. */
  readonly number: number;
  readonly program: ProgramId;
  /**
   * The run's name (NCBI's "Job Title"), without the white space around it; absent when
   * none was given. It is never part of the argv.
   */
  readonly title?: string;
  readonly argv: readonly string[];
  readonly query: InputSnapshot;
  readonly subject: InputSnapshot;
  readonly requestedThreads: number | 'auto';
  /** The counts file of `-window_masker_db`; absent without the option (BLASTN only). */
  readonly windowMasker?: CountsSnapshot;
  readonly queuedAt: number;
  /**
   * The group of runs queued together as separate searches (plan §5.2, Separate): one
   * run for each file. Undefined for a run queued alone.
   */
  readonly group?: RunGroup;
}

export interface RunGroup {
  readonly groupId: string;
  /** 1-based position of the run in its group, and the group's size. */
  readonly position: number;
  readonly size: number;
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
  /** What the run did with its WindowMasker counts file. */
  readonly windowMasker?: CountsRecord;
  readonly startedAt?: number;
  /** When each engine phase started. */
  readonly phaseTimes?: { readonly preparing?: number; readonly running?: number; readonly finalizing?: number };
  readonly endedAt?: number;
  readonly error?: string;
}

export function isTerminal(status: RunStatus): boolean {
  return TERMINAL_STATUSES.has(status);
}
