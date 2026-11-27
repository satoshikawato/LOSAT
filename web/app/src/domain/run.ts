// Run records (plan §5.2). A RunSnapshot is fixed when a job is queued and never
// changes; what actually happened is appended to the RunRecord.
import type { CountsFile, CountsUse, MessagePlace, PrintedKind, ProgramId, RecordKey, RunFile, SequenceKind, ShownIds } from 'losat/domain';
import type { CountsOrigin } from './counts';
import type { ToolId, ToolStage } from './tools';

export type RunStatus =
  | 'queued'
  | 'preparing'
  | 'running'
  | 'finalizing'
  | 'completed'
  | 'cancelled'
  | 'failed';

export const TERMINAL_STATUSES: ReadonlySet<RunStatus> = new Set(['completed', 'cancelled', 'failed']);

/**
 * What a file of the user's gave to an input that is not that file as given (line WQA's D22, 2026-11-23): its name (`pasted text`
 * for pasted text), the records of it that the input holds, and the first IDs of those it leaves out with their number.
 */
export interface InputPart {
  readonly name: string;
  /** The records of the part that the input holds. */
  readonly records: number;
  /** The IDs of the first records that the input leaves out, in the order of the file (at most MOST_LEFT_OUT_IDS). */
  readonly leftOut: readonly string[];
  /** The number of records that the input leaves out. */
  readonly leftOutCount: number;
}

export interface InputSnapshot {
  /**
   * Name passed as -query / -subject (plan §5.3): the user's own file name where the input is that file as given, else the
   * name of a file that the page offers with the run (`from`), whose bytes are exactly `bytes`.
   */
  readonly name: string;
  /**
   * What the input is made of when it is not one file of the user's as given: pasted text, several inputs searched together,
   * or a file with records left out. Absent where `name` is the user's file, whose bytes `bytes` are (line WQA's D22).
   */
  readonly from?: readonly InputPart[];
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
  /** ID and length of each record from the per-file record tables; the engine checks them at `register`. */
  readonly records: readonly RecordKey[];
  /**
   * Whether `bytes` was assembled from parts (several sources joined, or a source with records
   * left out; line WQA's D17). Absent for a single unmodified source.
   */
  readonly assembled?: boolean;
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
 * file is called. The coordinator holds the bytes until the run ends; the snapshot keeps what tells
 * the file: its names, its size and its SHA-256.
 */
export interface ConffileSnapshot {
  /** The name that the argv's `-conffile` gives (the engine's start reads it). */
  readonly name: string;
  /** The name of the file that the user chose. */
  readonly fileName: string;
  /** Bytes of the file that the run read. */
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes that the run read. */
  readonly sha256: string;
}

/** Metadata of one database filter list file (`-seqidlist`, `-taxidlist`, etc.) mounted for a `-db` run. */
export interface DatabaseFilterFileSnapshot {
  readonly flag: string;
  readonly name: string;
  readonly size: number;
}

/** Metadata of a local BLAST database (`-db` or BLASTN `-filtering_db`) fixed when a run is queued. */
export interface DatabaseSnapshot {
  /** Exact relative path passed to `-db` or `-filtering_db` (e.g. `mydb/nt_small`). */
  readonly name: string;
  readonly folder?: string;
  readonly kind: SequenceKind;
  readonly fileCount: number;
  readonly totalBytes: number;
  readonly taxdb: boolean;
  readonly sqlite: boolean;
  /** Relative database folder whose taxonomy files were mounted at the run directory root (`BLASTDB=<dir>` on CLI). */
  readonly blastdbEnv?: string;
  readonly filterFiles?: readonly DatabaseFilterFileSnapshot[];
}

/** Metadata of a file-based query (`-in_pssm` or `-in_msa`) fixed when a run is queued. */
export interface QueryFileSnapshot {
  readonly flag: '-in_pssm' | '-in_msa';
  readonly name: string;
  readonly size: number;
}

export interface RunSnapshot {
  readonly kind?: 'search';
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
  /** Present when the run used a PSSM (`-in_pssm`) or MSA (`-in_msa`) file query instead of FASTA `-query`. */
  readonly queryFile?: QueryFileSnapshot;
  readonly subject: InputSnapshot;
  /** Present when the run searched a local BLAST database (`-db`); absent for `-subject` FASTA runs. */
  readonly database?: DatabaseSnapshot;
  /** Present when a BLASTN run used a local repeats database (`-filtering_db`). */
  readonly filteringDb?: DatabaseSnapshot;
  readonly requestedThreads: number | 'auto';
  /** The counts file of `-window_masker_db`; absent without the option (BLASTN only). */
  readonly windowMasker?: CountsSnapshot;
  /** The registry file of the argv's `-conffile`, when one was chosen. */
  readonly conffile?: ConffileSnapshot;
  /** The argv has `-dryrun` (the engine's start says so): the run reads the arguments and the settings and searches nothing. */
  readonly dryRun?: true;
  /**
   * The run prints the command line's help or version in place of a search (the engine's start says so; line WQA's D14,
   * 2026-11-20): it reads no input (its query and subject are empty) and its text is `RunRecord.printed`.
   */
  readonly prints?: PrintedKind;
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
  /**
   * Where the run's messages went, as the engine said (`RuntimeInfo.messages`): its Warnings (stream 3), or, for TBLASTX
   * with -logfile, the file of -logfile or nowhere (line WQA's sweep, D7, 2026-11-20). Set when the run completes.
   */
  readonly messages?: MessagePlace;
  /**
   * The IDs that the run's outputs show for its records, where the engine read them from the
   * deflines (-parse_deflines; line WQA's sweep, D7, 2026-11-19): the names of the records in the results.
   */
  readonly shownIds?: ShownIds;
  /** What the run did with its WindowMasker counts file. */
  readonly windowMasker?: CountsRecord;
  /**
   * The files that the run made besides its outputs (the file of `-logfile`), as the engine wrote
   * them: kept for a completed run and for a failed one (TBLASTX's file then holds NCBI's error).
   */
  readonly runFiles?: readonly RunFile[];
  /** The text that the run printed in place of a search (`RunSnapshot.prints`), as the engine wrote it. */
  readonly printed?: Uint8Array;
  readonly startedAt?: number;
  /** When each engine phase started. */
  readonly phaseTimes?: { readonly preparing?: number; readonly running?: number; readonly finalizing?: number };
  readonly endedAt?: number;
  readonly error?: string;
  /**
   * What a failed run wrote to the command line's standard error before its error (engine stream 3:
   * its warnings, as NCBI's `Subject sequence contains no data` before `The average subject length is
   * too short`), as the engine wrote it; a completed run's are in its results (line WQA's sweep, D3,
   * 2026-11-19).
   */
  readonly diagnostics?: string;
  /**
   * The records that the engine actually read at `register` for an assembled input
   * (`InputSnapshot.assembled`) when reading the assembled bytes as one file differs from the
   * per-file record tables (`InputSnapshot.records`; line WQA's D17). Absent when they agree.
   */
  readonly readRecords?: {
    readonly query?: readonly RecordKey[];
    readonly subject?: readonly RecordKey[];
  };
}

/**
 * The records that a run actually searched for `role`: the engine's own reading when an assembled
 * input read differently as one file (`RunRecord.readRecords`), otherwise the snapshot's record table.
 */
export function runInputRecords(
  run: { readonly snapshot: RunSnapshot; readonly record: RunRecord },
  role: 'query' | 'subject',
): readonly RecordKey[] {
  return run.record.readRecords?.[role] ?? run.snapshot[role].records;
}

export function isTerminal(status: RunStatus): boolean {
  return TERMINAL_STATUSES.has(status);
}

// --- Tool runs (windowmasker, WW3b DESIGN §3-§4) -------------------------------------------

export type ToolFileOrigin =
  | { readonly kind: 'picked' }
  | { readonly kind: 'pasted' }
  | { readonly kind: 'inflated'; readonly from: string }
  | { readonly kind: 'combined'; readonly parts: readonly string[] }
  | { readonly kind: 'made'; readonly runNumber: number }
  | { readonly kind: 'archive'; readonly runNumber: number }
  | { readonly kind: 'saved'; readonly taxid: number };

/** Metadata of one file placed in a tool run's working directory (for Run details). */
export interface ToolDirectoryEntry {
  readonly name: string;
  readonly size: number;
  readonly origin: ToolFileOrigin;
}

/** One output stream or file kept in session storage for a tool run. */
export interface ToolOutputEntry {
  /** Key under which the Data worker keeps this block (`'stdout'`, `'stderr'`, or `'file:0'`, ...). */
  readonly key: string;
  /** Name as written in the argv (for example `wmasker.ascii`, `run.log`), or `'stdout'` / `'stderr'`. */
  readonly name: string;
  readonly size: number;
  /** Option(s) that named the file (`'-out'`, `'-logfile'`, `'-out, -logfile'`), if applicable. */
  readonly option?: string;
  /** True for the `-out` file of `-mk_counts` and `-convert` (offered to BLASTN and Mask). */
  readonly isCounts?: true;
}

export interface ToolRunSnapshot {
  readonly kind: 'tool';
  readonly runId: string;
  /** 1-based position in this working session, shared with search runs. */
  readonly number: number;
  readonly tool: ToolId;
  readonly stage: ToolStage;
  readonly argv: readonly string[];
  readonly directory: readonly ToolDirectoryEntry[];
  readonly prints?: PrintedKind;
  /** A masking run's `-outfmt` default from `-xmlhelp`, kept for the run's mask data (`domain/tool-argv.ts` `maskDataName`). */
  readonly outfmtDefault?: string;
  readonly queuedAt: number;
  readonly group?: RunGroup;
}

export interface ToolRunRecord {
  readonly exitStatus?: number;
  readonly stdout?: ToolOutputEntry;
  readonly stderr?: ToolOutputEntry;
  readonly files?: readonly ToolOutputEntry[];
  readonly runtimePath?: 'serial' | 'fake';
  readonly engineBuild?: string;
  readonly memory?: {
    readonly linearBytesBefore: number;
    readonly linearBytesAfter: number;
    readonly instanceRuns: number;
  };
  readonly startedAt?: number;
  readonly phaseTimes?: { readonly preparing?: number; readonly running?: number; readonly finalizing?: number };
  readonly endedAt?: number;
  readonly error?: string;
}

export type AnyRunSnapshot = RunSnapshot | ToolRunSnapshot;
export type AnyRunRecord = RunRecord | ToolRunRecord;

export function isToolSnapshot(snapshot: AnyRunSnapshot): snapshot is ToolRunSnapshot {
  return snapshot.kind === 'tool';
}
