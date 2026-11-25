// Data port (plan §5.4, §5.6): sources and their record tables, the outputs of runs, and
// the temporary storage that holds them. The implementation runs in the Data worker
// (src/infra/data-worker), so every argument and result is structured-cloneable. The
// outputs of a run arrive over the MessagePort that `openRun` returns
// (losat/ports: run-output.ts), not through these methods.
import type { CUSTOM_OUTPUT, FastaReader, OutputFormat, OutputPlace, RecordKey, RunOutput } from 'losat/domain';
import type { HspRecord, InputCheck } from 'losat/ports';
import type { DatasetRevision } from '../domain/dataset';
import type { HspTable } from '../domain/hsp-table';
import { withUnit } from '../domain/units';

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
   * Builds the record table of a source with the index scan of the engine's reader `reader` (its
   * program, role and data loaders) and the SHA-256 of every record. Rejects with the reader's
   * error, or its refusal in the CLI's words.
   */
  indexSource(sourceId: string, reader: FastaReader): Promise<DatasetRevision>;
  /** A new revision of the same record table that leaves out the records `excluded`. */
  reviseDataset(revisionId: string, excluded: readonly number[]): Promise<DatasetRevision>;
  /**
   * The engine input made of the included records of the revisions, in order. A revision
   * that includes every record contributes its source unchanged; a newline is added after
   * a source that does not end with one when another follows (plan §5.3; a CR-only stream
   * keeps CR so NCBI's line reader stays in CR mode, line WQA's D17). When assembled from
   * multiple sources or a record selection, `bytes` is scanned with `reader` (or the first
   * revision's reader) so `records` agrees with what the engine and the CLI read from `bytes`.
   */
  buildRunInput(revisionIds: readonly string[], reader?: FastaReader): Promise<RunInput>;
  /**
   * The engine's reading of the run input of the revisions with the reader `reader` (its program,
   * role and data loaders; losat/ports: input-check.ts). Resolves with the engine's verdict;
   * rejects only when the check itself cannot run.
   */
  checkInput(reader: FastaReader, revisionIds: readonly string[]): Promise<InputCheck>;
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
  /**
   * The outputs that the run could not write, by stream number: those that the engine said failed (its signal of
   * stream 9, docs/web/abi_v2.md §5: a format that failed alone, or the formats of a search that failed in a run of
   * two searches), each with where its error is, as the engine said. Each is empty (line WQA's sweep, D7, 2026-11-19
   * and 2026-11-20).
   */
  readonly failed: readonly FailedOutput[];
}

/** An output that a run could not write, and where its error is (`ResultSetRef.failed`). */
export interface FailedOutput {
  readonly output: RunOutput;
  readonly error: OutputPlace;
}

/** What `RunStore.keepFailedRun` keeps of a failed run: its outputs, and its warnings. */
export interface KeptRun {
  /** The outputs that the run wrote before its error; undefined where no kept output has bytes, and the run is dropped. */
  readonly result?: ResultSetRef;
  /** What the run wrote to the CLI's standard error before its error; '' where it wrote nothing. */
  readonly diagnostics: string;
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
  /**
   * Keeps a staged run that the engine failed after it began to write (its rejection says so: `runOutputFailedOf`),
   * once its `fail` arrived (losat/ports run-output.ts): the outputs that the engine did not discard, the bytes that the
   * CLI leaves in the output of the argv before its error, are committed as a completed run's are, to be read and
   * exported (`readOutput`, `readOutputRange`), and the warnings that the run wrote to stream 3 before its error come
   * back (UTF-8, as `readDiagnostics` reads a committed run's; docs/web/abi_v2.md §5). The run has no HSP records. Line
   * WQA's sweep, D15, 2026-11-20 (before, a failed run was dropped and only its warnings were read: D3, 2026-11-19).
   * Resolves to undefined, and drops the run, where its output arrived incomplete or could not be stored, and for a run
   * that is not staged.
   */
  keepFailedRun(runId: string): Promise<KeptRun | undefined>;
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

/**
 * Results kept in memory (no OPFS) may take at most this many bytes in a tab (S09). Memory
 * is then the only copy of the results, so running out of it would end the tab and lose
 * every result; the budget makes it the same recoverable failure as a full OPFS instead.
 * The engine's own memory comes on top (up to 4 GiB less one page, plan TD-7).
 */
export const MEMORY_RESULTS_CAPACITY_BYTES = 512 * 1024 * 1024;
/**
 * MEMORY_RESULTS_CAPACITY_BYTES as the page says it, "512 MB" (of 1024 × 1024 bytes): in the
 * message when it is reached, and in the note that says results are kept in memory. The number and
 * its unit are joined by a no-break space, so that the note never breaks "at most 512" / "MB" on a
 * narrow screen (line WQA's sweep, D9, 2026-11-19).
 */
export const MEMORY_RESULTS_CAPACITY = withUnit(MEMORY_RESULTS_CAPACITY_BYTES / (1024 * 1024), 'MB');

/** Removal of the temporary data that closed tabs left behind (plan §5.6). */
export type CleanupState =
  | { readonly state: 'pending' }
  | { readonly state: 'done'; readonly removedSessions: number }
  | {
      readonly state: 'unavailable';
      /** In LOSAT's words, never the browser's message. */
      readonly reason: string;
      /** The error behind it, "<name>: <message>", for developers (data-error). */
      readonly error?: string;
    };

export interface StorageInfo {
  readonly backend: StorageBackend;
  /** Why OPFS is not used, when the backend is memory: LOSAT's words, never the browser's message. */
  readonly fallbackReason?: string;
  /** The error behind fallbackReason, "<name>: <message>", for developers (data-error). */
  readonly fallbackError?: string;
  /** Bytes that this working session holds in temporary storage. */
  readonly sessionBytes: number;
  /** The browser's estimate for the whole site (`navigator.storage.estimate`); not a guarantee. */
  readonly estimate?: { readonly usage: number; readonly quota: number };
  readonly cleanup: CleanupState;
}

/** Inflated session file stored in session OPFS (`tmp/`) without keeping inflated bytes in page memory. */
export interface SessionFileRef {
  /** Storage token for `dropSessionFile`; undefined when the input file was plain and kept as-is. */
  readonly token?: string;
  readonly file: File;
  readonly name: string;
  readonly size: number;
  /** Original compressed file name when inflated (for example `genome.fa.gz`). */
  readonly inflatedFrom?: string;
}

/** One sealed output kept for a tool run in the Data worker's session storage. */
export interface StoredToolFile {
  readonly key: string;
  readonly name: string;
  readonly size: number;
}

/** Result of sealing a tool run's outputs in the Data worker. */
export interface ToolOutputsSummary {
  readonly stdoutSize: number;
  readonly stderrSize: number;
  readonly files: readonly StoredToolFile[];
}

/** Messages sent over the dedicated `MessagePort` returned by `DataGateway.openToolRun`. */
export type ToolOutputCommand =
  | { readonly id: number; readonly op: 'stdout'; readonly bytes: Uint8Array }
  | { readonly id: number; readonly op: 'stderr'; readonly bytes: Uint8Array }
  | { readonly id: number; readonly op: 'open'; readonly key: string; readonly name: string }
  | { readonly id: number; readonly op: 'write'; readonly key: string; readonly offset: number; readonly bytes: Uint8Array }
  | { readonly id: number; readonly op: 'truncate'; readonly key: string; readonly length: number }
  | { readonly id: number; readonly op: 'close'; readonly key: string }
  | { readonly id: number; readonly op: 'commit' }
  | { readonly id: number; readonly op: 'discard' };

export type ToolOutputReply =
  | { readonly id: number; readonly ok: true; readonly summary?: ToolOutputsSummary }
  | { readonly id: number; readonly ok: false; readonly error: { readonly name: string; readonly message: string } };

export interface ToolRunStore {
  /**
   * Stages a tool run in session storage (`tmp/`) and returns a dedicated `MessagePort` that accepts
   * `ToolOutputCommand` messages (with transferred buffers for zero-copy streaming).
   */
  openToolRun(runId: string): Promise<MessagePort>;
  /** Reads bytes `[start, end)` of a sealed output (`'stdout'`, `'stderr'`, or `'file:N'`) of a tool run. */
  readToolOutputRange(runId: string, key: string, start: number, end: number): Promise<Uint8Array>;
  /** Returns a disk/OPFS-backed `File` for a sealed tool output (for Download, "Use in BLASTN", or Mask). */
  getToolOutputFile(runId: string, key: string, fileName: string): Promise<File>;
  /**
   * Inflates a compressed file (`.gz`, `.bz2`, `.zst`) directly into session OPFS (`tmp/`) and returns
   * an OPFS-backed `File` under its plain name. A plain file is returned as-is without copying.
   */
  inflateToSessionFile(file: File): Promise<SessionFileRef>;
  /** Drops a session file created by `inflateToSessionFile`. Unknown tokens are ignored. */
  dropSessionFile(token: string): Promise<void>;
}

export interface DataGateway extends DatasetStore, RunStore, ToolRunStore {
  storageInfo(): Promise<StorageInfo>;
}
