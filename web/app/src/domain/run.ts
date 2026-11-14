// Run records (plan §5.2). A RunSnapshot is fixed when a job is queued and never
// changes; what actually happened is appended to the RunRecord.
import type { CountsFile, CountsUse, ProgramId, RecordKey, RunFile } from 'losat/domain';
import type { CountsOrigin } from './counts';

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
 * that the run's engine reads under its name. It lives in the working session only, as the
 * run does (the file of a taxonomy ID is kept by the browser, and the run reads it from there).
 */
export interface CountsSnapshot extends CountsFile {
  /** For `-window_masker_taxid`: where the file was downloaded from. */
  readonly origin?: CountsOrigin;
}

/** What a run did with its counts file (RunRecord.windowMasker). */
export interface CountsRecord extends CountsUse {
  /** For `-window_masker_taxid`: the taxonomy ID, the URL and the Last-Modified of the file that was read. */
  readonly origin?: CountsOrigin;
}

/**
 * The registry file of `-conffile` (Additional options; line RG E1w-3b), read when the run was
 * submitted: the run registers its bytes under the name that the argv gives, whatever the chosen
 * file is called. It lives in the working session only, as the run does.
 */
export interface ConffileSnapshot {
  /** The name that the argv's `-conffile` gives (the engine's start reads it). */
  readonly name: string;
  /** The name of the file that the user chose. */
  readonly fileName: string;
  readonly bytes: Uint8Array;
  /** Lower-case hex SHA-256 of `bytes`. */
  readonly sha256: string;
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
  /** The registry file of the argv's `-conffile`, when one was chosen. */
  readonly conffile?: ConffileSnapshot;
  /** The argv has `-dryrun` (the engine's start says so): the run reads the arguments and the settings and searches nothing. */
  readonly dryRun?: true;
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
  /**
   * The files that the run made besides its outputs (the file of `-logfile`), as the engine wrote
   * them: kept for a completed run and for a failed one (TBLASTX's file then holds NCBI's error).
   */
  readonly runFiles?: readonly RunFile[];
  readonly startedAt?: number;
  /** When each engine phase started. */
  readonly phaseTimes?: { readonly preparing?: number; readonly running?: number; readonly finalizing?: number };
  readonly endedAt?: number;
  readonly error?: string;
}

export function isTerminal(status: RunStatus): boolean {
  return TERMINAL_STATUSES.has(status);
}
