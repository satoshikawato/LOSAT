// Engine port: the application's only view of the search engine. Implementations live
// in src/infra; the Wasm implementation follows docs/web/abi_v2.md.
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

export interface RuntimeInfo {
  readonly path: 'threaded' | 'serial' | 'fake';
  readonly threads: number;
  readonly engineBuild: string;
  readonly fallbackReason?: string;
}

export interface EngineRunRequest {
  readonly runId: string;
  readonly argv: readonly string[];
  readonly query: Uint8Array;
  readonly subject: Uint8Array;
  readonly requestedThreads: number | 'auto';
}

/** Receives the outputs of one run. Supplied by the data layer. */
export interface RunSink {
  write(format: OutputFormat, chunk: Uint8Array): void;
  hits(records: readonly HspRecord[]): void;
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
  /** Resolves after every output was written to `sink`; rejects on failure or cancel. */
  run(request: EngineRunRequest, sink: RunSink, onPhase: (phase: EnginePhase) => void): Promise<RuntimeInfo>;
  /** Stops the run as soon as possible. Unknown or finished runs are ignored. */
  cancel(runId: string): void;
}
