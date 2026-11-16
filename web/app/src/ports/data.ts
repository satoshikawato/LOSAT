// Data port (plan §5.4, §5.6): sources and their record tables, the outputs of runs, and
// the temporary storage that holds them. The implementation runs in the Data worker
// (src/infra/data-worker), so every argument and result is structured-cloneable. The
// outputs of a run arrive over the MessagePort that `openRun` returns
// (losat/ports: run-output.ts), not through these methods.
import type { CUSTOM_OUTPUT, FastaParserKind, InputRole, OutputFormat, ProgramId, RecordKey, RunOutput } from 'losat/domain';
import type { HspRecord, InputCheck } from 'losat/ports';
import type { DatasetRevision } from '../domain/dataset';
import type { HspTable } from '../domain/hsp-table';

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
   * never copied into storage (design §6.1). A File named `.gz`, `.bz2` or `.zst` is inflated first
   * by the engine's decoders (losat/ports: inflate.ts) and kept under its name without the suffix; its ref
   * carries that name and the inflated size. A file whose bytes are not what its name says is
   * refused with the engine's sentence.
   */
  addSource(file: File): Promise<SourceRef>;
  /**
   * Builds the record table of a source with the index scan of `parser` and the SHA-256
   * of every record. Rejects with the parser's error.
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
   * (losat/ports: input-check.ts). Resolves with the engine's verdict; rejects only when the
   * check itself cannot run.
   */
  checkInput(program: ProgramId, role: InputRole, revisionIds: readonly string[]): Promise<InputCheck>;
  /** The first `maxBytes` bytes of a source, for its preview. */
  previewSource(sourceId: string, maxBytes: number): Promise<Uint8Array>;
  /**
   * Forgets a source and every revision of it, giving back their memory; later calls that
   * name them reject as for an unknown id. The caller drops only what no one uses
   * (src/application/source-refs.ts). An unknown id is ignored.
   */
  dropSource(sourceId: string): Promise<void>;
}

export interface ResultSetRef {
  readonly runId: string;
  /** The bytes of each output of the run: the formats that every run writes, and the output of its -outfmt where it has one. */
  readonly byteLengths: Readonly<Record<OutputFormat, number> & Partial<Record<typeof CUSTOM_OUTPUT, number>>>;
  readonly hitCount: number;
}

/** The committed runs whose stored outputs cannot be read back (`RunStore.checkStoredResults`), in the order of the runs. */
export interface StoredResultsCheck {
  /** The browser removed their stored outputs: it answers that they are not found. */
  readonly removed: readonly string[];
  /** Their stored outputs failed to read for another reason. */
  readonly unreadable: readonly string[];
}

/**
 * Outputs of runs. A run is staged from `openRun` until `commitRun`; only committed runs
 * can be read. A staged run that is discarded (cancel, failure) leaves nothing behind.
 */
export interface RunStore {
  /**
   * Stages a run that writes `outputs` (`runOutputs` of its argv; the formats that every run
   * writes always) and returns the port that its engine writes to. An output that the run
   * does not write fails it at the commit.
   */
  openRun(runId: string, outputs?: readonly RunOutput[]): Promise<MessagePort>;
  /**
   * Waits for the `end` of the run's output and commits it. Rejects, and discards the run,
   * if the output is incomplete or the storage ran out (the error then says so).
   */
  commitRun(runId: string): Promise<ResultSetRef>;
  /** Drops a staged run. Committed and unknown runs are left alone. */
  discardRun(runId: string): Promise<void>;
  readOutput(runId: string, output: RunOutput): Promise<Uint8Array>;
  /** Bytes [start, end) of one output of a committed run (an HSP's row or section). */
  readOutputRange(runId: string, output: RunOutput, start: number, end: number): Promise<Uint8Array>;
  readHits(runId: string): Promise<readonly HspRecord[]>;
  /**
   * The HSP records of a committed run as columns, without the aligned sequences
   * (domain/hsp-table.ts). The Data worker builds it, so the records never cross to the
   * UI thread as objects (design §10.2).
   */
  readHitTable(runId: string): Promise<HspTable>;
  readDiagnostics(runId: string): Promise<string>;
  deleteRun(runId: string): Promise<void>;
  /**
   * Whether the stored outputs of every committed run are still there: it opens the HSP records
   * of each run (the first block that opening the run asks for) and reads its first byte, a few
   * runs at a time, so it stays quick for hundreds of runs. The browser can remove them under an
   * open tab (WA2-SC-05); the check after the page was hidden (application/attention.ts) asks.
   */
  checkStoredResults(): Promise<StoredResultsCheck>;
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
