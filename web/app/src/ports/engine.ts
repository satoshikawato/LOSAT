// Engine port: the application's only view of the search engine. Implementations live
// in src/infra; the Wasm implementation follows docs/web/abi_v2.md.
import type { RecordKey } from '../domain/dataset';
import type { OutputFormat } from '../domain/output-format';
import type { ProgramId } from '../domain/programs';

export interface ParameterDescription {
  readonly flag: string;
  readonly help: string;
  readonly takesValue: boolean;
  readonly defaultValue?: string;
  readonly choices?: readonly string[];
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
}

/** One input of a run: the exact FASTA bytes of the run snapshot and their record table. */
export interface EngineInput {
  readonly bytes: Uint8Array;
  /** Lower-case hex SHA-256 of `bytes`. With `revisionIds`, it identifies a retained subject. */
  readonly sha256: string;
  /** The dataset revisions whose records make up `bytes` (RunSnapshot). */
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

/** Thrown (as a rejection) by `run` after `cancel` was called for that run. */
export class RunCancelledError extends Error {
  constructor(readonly runId: string) {
    super(`run ${runId} was cancelled`);
    this.name = 'RunCancelledError';
  }
}

export interface EngineGateway {
  describe(program: ProgramId): Promise<ProgramDescription>;
  validate(argv: readonly string[]): Promise<ValidationResult>;
  /**
   * Runs one search and sends its outputs to `output` (ports/run-output.ts): the chunks of
   * every stream, then `end`, all posted before the promise resolves. `output` comes from
   * `DataGateway.openRun` and may be transferred to a worker. Rejects on failure or cancel
   * without sending `end`; the data layer then discards what arrived.
   */
  run(request: EngineRunRequest, output: MessagePort, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo>;
  /** Stops the run as soon as possible. Unknown or finished runs are ignored. */
  cancel(runId: string): void;
}
