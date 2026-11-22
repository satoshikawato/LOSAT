// Engine port: the application's only view of the search engine. The Wasm implementation
// (src/infra/engine-worker/gateway.ts) follows docs/web/abi_v2.md; LOSAT Web's development
// build has a FakeEngine of its own (web/app/src/infra/fake/fake-engine.ts).
import type { RecordKey, ShownIds } from '../domain/dataset.js';
import type { OutputFormat } from '../domain/output-format.js';
import type { ResolvedDefault } from '../domain/parameters.js';
import type { ProgramId } from '../domain/programs.js';
import type { CountsFile, CountsUse } from '../domain/counts.js';
import type { MessagePlace, PrintedKind, RegistryFile, RunFile } from '../domain/run-files.js';

export interface ParameterDescription {
  readonly flag: string;
  readonly help: string;
  readonly takesValue: boolean;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
  /** The default that the engine resolves when the option is left out (`describe`'s `defaults`). */
  readonly resolvedDefault?: ResolvedDefault;
}

export interface ProgramDescription {
  readonly program: ProgramId;
  /** Output formats the program supports; a run writes exactly these. */
  readonly formats: readonly OutputFormat[];
  readonly parameters: readonly ParameterDescription[];
  /** Genetic codes the engine accepts for -query_gencode, if the program has it. */
  readonly queryGencodes?: readonly number[];
  /** Genetic codes the engine accepts for -db_gencode, if the program has it. */
  readonly subjectGencodes?: readonly number[];
}

/**
 * The verdict of `validate`. A valid argv also says whether its run reads the inputs with NCBI's
 * data loaders on (`dataLoaders`, ABI v2 *validation*, version 8): off under `[BLAST]
 * DATA_LOADERS=none` of its `-conffile`, where the reader reads a first line that may be a Seq-id as
 * data, as the CLI does (line WQA's D16, 2026-11-20). A refused argv says it too when the engine's
 * start read the registry before the option that it refuses (line WQA's D13 round 3, L1).
 */
export type ValidationResult =
  | { readonly ok: true; readonly dataLoaders: boolean }
  | { readonly ok: false; readonly message: string; readonly dataLoaders?: boolean };

/**
 * What NCBI's application reads of an argv before its arguments (docs/web/abi_v2.md §4
 * `losat_web2_start`, §9 *start*; line RG item E1w-3b): the name of the last `-conffile` (the
 * registry file that `validate` and `run` read, which the host registers under that name), the
 * name of `-logfile` (NCBI's scan of it: a run gives the file of a name that makes a file among its
 * `runFiles`), the name of `-export_search_strategy` (the search strategy file, which a run gives
 * among its `runFiles` too; line RG item F1w), and whether the argv has `-dryrun` (a run reads the
 * arguments and the settings and searches nothing). Null where the argv names none.
 */
export interface StartInfo {
  readonly conffile: string | null;
  readonly logFile: string | null;
  readonly strategyFile: string | null;
  readonly dryRun: boolean;
  /**
   * Whether a run of the argv takes its query and subject read together for it
   * (`losat_web2_register_inputs`, ABI v2 §4: an argv with -parse_deflines, which reads the ids of
   * the deflines in one scope for both, and is not a dry run), not registered one at a time. The
   * Engine worker registers a run's inputs by it.
   */
  readonly registerInputs: boolean;
  /**
   * What a run of the argv prints in place of a search, as the command line does (ABI v2 §5 stream 11): its help or
   * its version; null for a search. Such a run reads no inputs (line WQA's D14, 2026-11-20).
   */
  readonly prints: PrintedKind | null;
}

/** The engine's *start* of an argv, or the CLI's message for an argv that it cannot start (a program that LOSAT Web does not run). */
export type StartResult = { readonly ok: true; readonly start: StartInfo } | { readonly ok: false; readonly message: string };

/**
 * One HSP of a finished run (docs/web/abi_v2.md "HSP record"). Field names follow the
 * Rust `Hit` / `PairwiseHit` fields so that no mapping layer is needed.
 */
export interface HspRecord {
  /** 0-based position in the final hit list of the run; the same HSP in every format. */
  readonly index: number;
  readonly q_idx: number;
  readonly s_idx: number;
  /** 0-based position among the HSPs of its query. */
  readonly rank: number;
  readonly raw_score: number;
  readonly bit_score: number;
  readonly e_value: number;
  readonly q_start: number;
  readonly q_end: number;
  readonly s_start: number;
  readonly s_end: number;
  readonly query_frame: number | null;
  readonly subject_frame: number | null;
  readonly subject_length: number | null;
  readonly query_aligned: string | null;
  readonly subject_aligned: string | null;
  /** Byte range [start, end) of this HSP's row in the outfmt 6 text, or null if not shown. */
  readonly out6: readonly [number, number] | null;
  /** Byte range [start, end) of this HSP's section in the outfmt 0 text, or null if not shown. */
  readonly out0: readonly [number, number] | null;
  /** Byte range [start, end) of the heading of this HSP's subject in the outfmt 0 text, or null. */
  readonly out0_subject: readonly [number, number] | null;
}

export type EnginePhase = 'preparing' | 'running' | 'finalizing';

/** Linear memory of the engine instance around one search (plan §5.5, G12). */
export interface EngineMemory {
  readonly linearBytesBefore: number;
  readonly linearBytesAfter: number;
  /** Searches the instance has run, including this one. */
  readonly instanceRuns: number;
}

export interface RuntimeInfo {
  readonly path: 'threaded' | 'serial' | 'fake';
  readonly threads: number;
  readonly engineBuild: string;
  /** Why the search ran on the serial module although more threads were requested. */
  readonly fallbackReason?: string;
  /**
   * The runtime (Engine worker) that ran the search. It changes when the runtime is ended
   * (cancel) or renewed; messages of an older runtime are dropped.
   */
  readonly runtimeGeneration?: number;
  readonly memory?: EngineMemory;
  /** The engine searched the subject that it held from an earlier search (plan §4.6 R1). */
  readonly subjectRetained?: boolean;
  /**
   * Where the run's messages go, as the engine said (ABI v2 §5, stream 10): the messages that the CLI writes to its
   * standard error without -logfile. The diagnostics, or, for TBLASTX with -logfile, the file of -logfile or nowhere,
   * as on the CLI (line WQA's sweep, D7, 2026-11-20).
   */
  readonly messages: MessagePlace;
  /**
   * The IDs that the run's outputs show for its records, where the engine read them from the
   * deflines (-parse_deflines; line WQA's sweep, D7, 2026-11-19).
   */
  readonly shownIds?: ShownIds;
  /** The WindowMasker counts file of the run, when the argv has `-window_masker_db`. */
  readonly windowMasker?: CountsUse;
  /**
   * The files that the run made besides its outputs (ABI v2 stream 8: the file of `-logfile`), when
   * it made any. A run that fails gives its files with its rejection (`runFilesOf`).
   */
  readonly runFiles?: readonly RunFile[];
  /**
   * The text that the run printed in place of a search (ABI v2 stream 11: the command line's standard output), for an
   * argv whose start `prints` (line WQA's D14, 2026-11-20).
   */
  readonly printed?: Uint8Array;
}

/** One input of a run: the exact FASTA bytes of the run snapshot and their record table. */
export interface EngineInput {
  readonly bytes: Uint8Array;
  /** Lower-case hex SHA-256 of `bytes`. With `revisionIds`, it identifies a retained subject. */
  readonly sha256: string;
  /** The dataset revisions whose records made up `bytes` (RunSnapshot): ids only, which the data layer may have dropped. */
  readonly revisionIds: readonly string[];
  /**
   * ID and length of each record in `bytes`, from the data layer's index scan. After
   * `register`, the engine compares the records its own parser read with these and fails
   * the run with InputMismatchError, before searching, if they differ (plan §5.4).
   */
  readonly records: readonly RecordKey[];
}

export interface EngineRunRequest {
  readonly runId: string;
  readonly argv: readonly string[];
  readonly query: EngineInput;
  readonly subject: EngineInput;
  readonly requestedThreads: number | 'auto';
  /** The file of `-window_masker_db`, when the argv has it. */
  readonly windowMasker?: CountsFile;
  /** The registry file of `-conffile`, under the name that the engine's start reads of the argv (`StartInfo.conffile`). */
  readonly conffile?: RegistryFile;
}

/** Thrown (as a rejection) by `run` when `register` read other records than the record table. */
export class InputMismatchError extends Error {
  constructor(
    readonly role: 'query' | 'subject',
    detail: string,
  ) {
    super(`The ${role} records that the engine read differ from the record table: ${detail}`);
    this.name = 'InputMismatchError';
  }
}

/**
 * Thrown (as a rejection) by `validate` or `run` when the file given with it (the registry file of
 * `-conffile`, a WindowMasker counts file) could not be read, as a browser refuses a file chosen from the
 * disk that changed after it was chosen: the message is the browser's. An error that crosses a worker
 * keeps its name only (`isFileReadError`).
 */
export class FileReadError extends Error {
  constructor(message: string, options?: { cause?: unknown }) {
    super(message, options);
    this.name = 'FileReadError';
  }
}

/** Whether an error is a FileReadError, also one rebuilt from its name after it crossed a worker. */
export function isFileReadError(error: unknown): error is Error {
  return error instanceof Error && error.name === 'FileReadError';
}

/**
 * A part of LOSAT could not be loaded, or stopped: an engine file that could not be downloaded, a worker that could not
 * start or that ended (line WQA's sweep, D1 and D2, 2026-11-19). It is no verdict on the inputs or the options, and a
 * page says it as LOSAT's own failure. Its message names the part and why, as a clause. An error that crosses a worker
 * keeps its name only (`isUnavailable`).
 */
export class UnavailableError extends Error {
  constructor(message: string, options?: { cause?: unknown }) {
    super(message, options);
    this.name = 'UnavailableError';
  }
}

/** Whether an error is an UnavailableError, also one rebuilt from its name after it crossed a worker. */
export function isUnavailable(error: unknown): error is Error {
  return error instanceof Error && error.name === 'UnavailableError';
}

/** Thrown (as a rejection) by `run` after `cancel` was called for that run. */
export class RunCancelledError extends Error {
  constructor(readonly runId: string) {
    super(`run ${runId} was cancelled`);
    this.name = 'RunCancelledError';
  }
}

/**
 * The files that a failed run made besides its outputs (ABI v2 §5: stream 8 comes also when a run
 * fails; TBLASTX's file of `-logfile` then holds NCBI's error), which the rejection of
 * `EngineGateway.run` carries (`withRunFiles`); undefined for an error without files, as a trap's or
 * a cancel's.
 */
export function runFilesOf(error: unknown): readonly RunFile[] | undefined {
  const files = typeof error === 'object' && error !== null ? (error as { readonly runFiles?: unknown }).runFiles : undefined;
  return Array.isArray(files) ? (files as readonly RunFile[]) : undefined;
}

/** Gives the error of a failed run the files that the run made (`runFilesOf`); its name and message stay. */
export function withRunFiles<E extends Error>(error: E, runFiles: readonly RunFile[]): E {
  return Object.assign(error, { runFiles });
}

/**
 * Whether the engine ended the output of a failed run with `fail` (ports/run-output.ts), which the
 * rejection of `EngineGateway.run` then says (`withRunOutputFailed`): the data layer has what the run
 * wrote to the CLI's standard error before its error (ABI v2 §5: stream 3, which the host keeps with a
 * failed run as it keeps stream 8), as NCBI's `Subject sequence contains no data` before `The average
 * subject length is too short` (line WQA's sweep, D3, 2026-11-19), and reads it before it discards
 * the run. False for a run that failed before it began to write, and for a trap's or a cancel's error.
 */
export function runOutputFailedOf(error: unknown): boolean {
  return typeof error === 'object' && error !== null && (error as { readonly outputFailed?: unknown }).outputFailed === true;
}

/** Says that the engine ended the output of the failed run with `fail` (`runOutputFailedOf`); the error's name and message stay. */
export function withRunOutputFailed<E extends Error>(error: E): E {
  return Object.assign(error, { outputFailed: true as const });
}

export interface EngineGateway {
  describe(program: ProgramId): Promise<ProgramDescription>;
  /**
   * The engine's *start* of an argv (`StartInfo`): a host calls it before `validate`, to learn which
   * registry file the argv reads, the name of its `-logfile` file and whether it is a dry run.
   */
  start(argv: readonly string[]): Promise<StartResult>;
  /**
   * Checks an argv as the command line does. `registry` is the file of the argv's `-conffile`
   * (`StartInfo.conffile`), which the check reads; it is given to the engine for this check only.
   * Without it, an argv with `-conffile` gets the CLI's error for a file that is not there.
   */
  validate(argv: readonly string[], registry?: RegistryFile): Promise<ValidationResult>;
  /**
   * Runs one search and sends its outputs to `output` (ports/run-output.ts): the chunks of
   * every stream, then `end`, all posted before the promise resolves. `output` comes from
   * `DataGateway.openRun` and may be transferred to a worker. Rejects on failure or cancel
   * without sending `end`; the data layer then discards what arrived. The files that a failed run
   * made besides its outputs come with the rejection (`runFilesOf`). A run that the engine failed
   * after it began to write sends `fail` in place of `end` before it rejects, and its rejection says
   * so (`runOutputFailedOf`): the data layer reads the warnings that it wrote before its error from
   * what arrived, then discards it.
   */
  run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo>;
  /** Stops the run as soon as possible. Unknown or finished runs are ignored. */
  cancel(runId: string): void;
  /**
   * A hint that a search like `hint` may come soon: the engine may start what its first
   * search would start (its runtime and threads) now instead. Optional; it changes no
   * search and reports nothing (a failure shows in the search).
   */
  prepare?(hint: EnginePrepareHint): void;
}

/** What `EngineGateway.prepare` may prepare for. */
export interface EnginePrepareHint {
  readonly requestedThreads: number | 'auto';
  /** About the query and subject bytes of the search. */
  readonly inputBytes: number;
  /**
   * The argv of the search (program first, without -num_threads), when known: the engine may
   * run it once on small sequences of its own, so that the browser compiles the engine code
   * that the search calls before it comes.
   */
  readonly argv?: readonly string[];
}
