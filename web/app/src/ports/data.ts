// Data port (plan §5.4, §5.6): sources and their record tables, the outputs of runs, and
// the temporary storage that holds them. The implementation runs in the Data worker
// (src/infra/data-worker), so every argument and result is structured-cloneable. The
// outputs of a run arrive over the MessagePort that `openRun` returns
// (ports/run-output.ts), not through these methods.
import type { DatasetRevision, FastaParserKind, RecordKey } from '../domain/dataset';
import type { HspTable } from '../domain/hsp-table';
import type { OutputFormat } from '../domain/output-format';
import type { InputRole, ProgramId } from '../domain/programs';
import type { HspRecord } from './engine';
import type { InputCheck } from './input-check';

export interface SourceRef {
  readonly sourceId: string;
  /** File name; pasted text arrives as a File named `query.fa` or `subject.fa` (plan §5.3). */
  readonly name: string;
  readonly size: number;
}

/** The engine input of one role of a run. */
export interface RunInput {
  /** The included original records of the revisions, in order (PD-LOSAT-WEB-APP-BOUNDARY 2.1). */
  readonly bytes: Uint8Array;
  /** Lower-case hex SHA-256 of `bytes`. */
  readonly sha256: string;
  /** ID and length of each record in `bytes`. */
  readonly records: readonly RecordKey[];
}

export interface DatasetStore {
  /**
   * Keeps a reference to `file`. Its bytes are read with `File.slice` when needed and are
   * never copied into storage (design §6.1).
   */
  addSource(file: File): Promise<SourceRef>;
  /**
   * Builds the record table of a source with the index scan of `parser` and the SHA-256
   * of every record. Rejects with the scan's error.
   */
  indexSource(sourceId: string, parser: FastaParserKind): Promise<DatasetRevision>;
  /** A new revision of the same record table that leaves out the records `excluded`. */
  reviseDataset(revisionId: string, excluded: readonly number[]): Promise<DatasetRevision>;
  /**
   * The engine input made of the included records of the revisions, in order. A revision
   * that includes every record contributes its source unchanged; a newline is added after
   * a source that does not end with one when another follows (plan §5.3).
   */
  buildRunInput(revisionIds: readonly string[]): Promise<RunInput>;
  /**
   * The engine's reading of the run input of the revisions for `program` and `role`
   * (ports/input-check.ts). Resolves with the engine's verdict; a refusal whose message
   * names a line of the run input (NCBI's line numbers) also gives the position of the
   * record that holds the line. Rejects only when the check itself cannot run.
   */
  checkInput(program: ProgramId, role: InputRole, revisionIds: readonly string[]): Promise<InputCheck>;
  /** The first `maxBytes` bytes of a source, for its preview. */
  previewSource(sourceId: string, maxBytes: number): Promise<Uint8Array>;
}

export interface ResultSetRef {
  readonly runId: string;
  readonly byteLengths: Readonly<Record<OutputFormat, number>>;
  readonly hitCount: number;
}

/**
 * Outputs of runs. A run is staged from `openRun` until `commitRun`; only committed runs
 * can be read. A staged run that is discarded (cancel, failure) leaves nothing behind.
 */
export interface RunStore {
  /** Stages a run and returns the port that its engine writes to. */
  openRun(runId: string): Promise<MessagePort>;
  /**
   * Waits for the `end` of the run's output and commits it. Rejects, and discards the run,
   * if the output is incomplete or the storage ran out (the error then says so).
   */
  commitRun(runId: string): Promise<ResultSetRef>;
  /** Drops a staged run. Committed and unknown runs are left alone. */
  discardRun(runId: string): Promise<void>;
  readOutput(runId: string, format: OutputFormat): Promise<Uint8Array>;
  /** Bytes [start, end) of one output of a committed run (an HSP's row or section). */
  readOutputRange(runId: string, format: OutputFormat, start: number, end: number): Promise<Uint8Array>;
  readHits(runId: string): Promise<readonly HspRecord[]>;
  /**
   * The HSP records of a committed run as columns, without the aligned sequences
   * (domain/hsp-table.ts). The Data worker builds it, so the records never cross to the
   * UI thread as objects (design §10.2).
   */
  readHitTable(runId: string): Promise<HspTable>;
  readDiagnostics(runId: string): Promise<string>;
  deleteRun(runId: string): Promise<void>;
}

export type StorageBackend = 'opfs' | 'memory';

/** Removal of the temporary data that closed tabs left behind (plan §5.6). */
export type CleanupState =
  | { readonly state: 'pending' }
  | { readonly state: 'done'; readonly removedSessions: number }
  | { readonly state: 'unavailable'; readonly reason: string };

export interface StorageInfo {
  readonly backend: StorageBackend;
  /** Why OPFS is not used, when the backend is memory. */
  readonly fallbackReason?: string;
  /** Bytes that this working session holds in temporary storage. */
  readonly sessionBytes: number;
  /** The browser's estimate for the whole site (`navigator.storage.estimate`); not a guarantee. */
  readonly estimate?: { readonly usage: number; readonly quota: number };
  readonly cleanup: CleanupState;
}

export interface DataGateway extends DatasetStore, RunStore {
  storageInfo(): Promise<StorageInfo>;
}
