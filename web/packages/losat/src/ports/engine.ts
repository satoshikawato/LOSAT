// Engine port: the application's only view of the search engine. The Wasm implementation
// (src/infra/engine-worker/gateway.ts) follows docs/web/abi_v2.md; LOSAT Web's development
// build has a FakeEngine of its own (web/app/src/infra/fake/fake-engine.ts).
import type { RecordKey } from '../domain/dataset.js';
import type { OutputFormat } from '../domain/output-format.js';
import type { ResolvedDefault } from '../domain/parameters.js';
import type { ProgramId } from '../domain/programs.js';
import type { CountsFile, CountsUse } from '../domain/counts.js';
import type { RegistryFile, RunFile } from '../domain/run-files.js';

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

export type ValidationResult = { readonly ok: true } | { readonly ok: false; readonly message: string };

/**
 * What NCBI's application reads of an argv before its arguments (docs/web/abi_v2.md §4
 * `losat_web2_start`, §9 *start*; line RG item E1w-3b): the name of the last `-conffile` (the
 * registry file that `validate` and `run` read, which the host registers under that name), the
 * name of `-logfile` (NCBI's scan of it: a run gives the file of a name that makes a file among its
 * `runFiles`), and whether the argv has `-dryrun` (a run reads the arguments and the settings and
 * searches nothing). Null where the argv names none.
 */
export interface StartInfo {
  readonly conffile: string | null;
  readonly logFile: string | null;
  readonly dryRun: boolean;
}

/** The engine's *start* of an argv, or the CLI's message for an argv that it cannot start (a version word). */
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
  /** The WindowMasker counts file of the run, when the argv has `-window_masker_db`. */
  readonly windowMasker?: CountsUse;
  /**
   * The files that the run made besides its outputs (ABI v2 stream 8: the file of `-logfile`), when
   * it made any. A run that fails gives its files with its rejection (`runFilesOf`).
   */
  readonly runFiles?: readonly RunFile[];
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
   * made besides its outputs come with the rejection (`runFilesOf`).
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
