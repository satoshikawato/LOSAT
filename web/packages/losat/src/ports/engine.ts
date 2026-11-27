// Engine port: the application's only view of the search engine. The Wasm implementation
// (src/infra/engine-worker/gateway.ts) follows docs/web/abi_v2.md; LOSAT Web's development
// build has a FakeEngine of its own (web/app/src/infra/fake/fake-engine.ts).
import type { RecordKey, ShownIds } from '../domain/dataset.js';
import type { OutputFormat } from '../domain/output-format.js';
import type { ResolvedDefault } from '../domain/parameters.js';
import type { ProgramId } from '../domain/programs.js';
import type { CountsFile, CountsUse } from '../domain/counts.js';
import type { DirectoryFile, MessagePlace, PrintedKind, RegistryFile, RunFile } from '../domain/run-files.js';
import type { ToolDescription } from '../domain/tool-description.js';

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
 * data, as the CLI does (line WQA's D16, 2026-11-20), and whether `-parse_deflines` is active
 * (`parseDeflines`, ABI v2 *validation*, version 15; line RG-XT item X7-5, 2026-11-25). A refused
 * argv says `dataLoaders` when the engine's start read the registry before the option that it
 * refuses (line WQA's D13 round 3, L1), and says `parseDeflines` only when the argv parsed
 * (omitted on parse errors so the host keeps its previous value).
 */
export type ValidationResult =
  | { readonly ok: true; readonly dataLoaders: boolean; readonly parseDeflines: boolean }
  | { readonly ok: false; readonly message: string; readonly dataLoaders?: boolean; readonly parseDeflines?: boolean };

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
  /**
   * The command line's normalized file name under which a run gives the `-logfile` file on stream 8
   * (`logs` for `logs/`, `run.log` for `./run.log` and ` run.log`, `/tmp/run.log` outside the
   * working tree, `/dev/stdout` for `/dev/stdout`; null where `-logfile` never makes a file).
   */
  readonly logFileName?: string | null;
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
  /**
   * Whether a run of the argv searches the database of `-db`, which the engine opens in the run's working directory
   * (`EngineRunRequest.files`), and takes no subject (line WQA D12a, 2026-11-20).
   */
  readonly database: boolean;
  /**
   * The query file name (`-in_pssm` or `-in_msa`) in the run's working directory (`EngineRunRequest.files`),
   * when the run reads its query from a file instead of a FASTA query input (query handle 0; PSW1).
   */
  readonly queryFile?: string | null;
  /** The normalized host file name for `-out_pssm` on stream 8, if any (PSW1). */
  readonly outPssm?: string | null;
  /** The normalized host file name for `-out_ascii_pssm` on stream 8, if any (PSW1). */
  readonly outAsciiPssm?: string | null;
  /** Whether `-save_each_pssm` numbers PSSM output files (`F.1`, `F.2`, ...; PSW1). */
  readonly saveEachPssm?: boolean;
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
  /** Query sequence length reported by the engine (used for `-in_pssm`/`-in_msa` runs without FASTA query records). */
  readonly query_length?: number | null;
  readonly subject_length: number | null;
  readonly query_aligned: string | null;
  readonly subject_aligned: string | null;
  /** Byte range [start, end) of this HSP's row in the outfmt 6 text, or null if not shown. */
  readonly out6: readonly [number, number] | null;
  /** Byte range [start, end) of this HSP's section in the outfmt 0 text, or null if not shown. */
  readonly out0: readonly [number, number] | null;
  /** Byte range [start, end) of the heading of this HSP's subject in the outfmt 0 text, or null. */
  readonly out0_subject: readonly [number, number] | null;
  /** 1-based PSI-BLAST iteration round number, present for multi-round `psiblast` searches. */
  readonly round?: number;
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
  /**
   * The records that the engine read at `register` for an assembled input (`EngineInput.assembled`)
   * when they differ from the request's per-file record table (`EngineInput.records`; line WQA's D17).
   * Absent when the engine read the same records as the table.
   */
  readonly readRecords?: {
    readonly query?: readonly RecordKey[];
    readonly subject?: readonly RecordKey[];
  };
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
   * the run with InputMismatchError, before searching, if they differ (plan §5.4), unless
   * `assembled` is true.
   */
  readonly records: readonly RecordKey[];
  /**
   * Whether `bytes` was assembled from parts (several sources joined, or a source with records
   * left out; line WQA's D17). Because NCBI's line reader keeps line-end state across parts,
   * reading `bytes` as one file can differ from the per-file tables in `records`: the engine
   * accepts what `register` read and returns it in `RuntimeInfo.readRecords`.
   */
  readonly assembled?: boolean;
}

export interface EngineRunRequest {
  readonly runId: string;
  readonly argv: readonly string[];
  /** The query; none for a run of `-in_pssm` or `-in_msa` (`StartInfo.queryFile`), which reads its query from `files`. */
  readonly query?: EngineInput;
  /** The subject; none for a run of `-db` (`StartInfo.database`), which searches the database of `files`. */
  readonly subject?: EngineInput;
  /**
   * The files of the run's working directory (docs/web/abi_v2.md §3), which the engine opens by the names that the argv
   * gives, as the command line opens them in its directory: the files of the BLAST database of `-db`, `taxdb.*`, the
   * lists of its filters (`-seqidlist`, `-taxidlist`, `-gilist` and their negatives), and BLASTN's database of repeats
   * (`-filtering_db`, with `-db` or a `subject`; line RG item D4). The engine reads each in place,
   * the bytes that it asks for at a time; none is copied whole (line WQA D12a, 2026-11-20).
   */
  readonly files?: readonly DirectoryFile[];
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

// --- tools (ABI v2 version 13, `losat_web2_tool`; line NP WW3a, 2026-11-22) -------------------------------------

/**
 * Where a file that a tool run writes goes, as it is written (D-WW3-6): a tool's output may be larger than the
 * browser's memory (human counts are 724 MB, a genome masked as FASTA 3.1 GB), so the engine hands it over in
 * chunks. The calls come in order, each after the promise of the one before has settled: the file's writes and
 * truncations in the order that the run made them, then `close`. A `chunk` is the sink's own (it is not used again by
 * the host). The engine sends a file's bytes in frames of up to 1 MiB, joining writes that follow each other: a file that
 * grows in small writes (a `-logfile`'s messages) reaches the sink when a frame fills, before the run writes to another
 * file or cuts one, and when the run ends, so a small file may arrive only then. The standard error is sent at each write.
 *
 * A file has POSIX semantics: most writes are appends (`offset` is the length so far), but a run that opens one file
 * twice (`-out X -logfile X`) writes one file through two handles, each with its own position, so a write may come at
 * an offset before the end, and one past the end leaves zero bytes between (a hole). A `FileSystemWritableFileStream`
 * does both: `write({ type: 'write', position: offset, data: chunk })` and `truncate(length)`.
 */
export interface ToolFileSink {
  /** Writes `chunk` at `offset` of the file (bytes from 0); a write past the end leaves zero bytes between. */
  write(chunk: Uint8Array, offset: number): void | Promise<void>;
  /**
   * The file's length becomes `length`: the file is cut, or, when it is shorter, zero bytes are added. The run asks it
   * when it opens a file again for writing, which empties it (`length` 0).
   */
  truncate(length: number): void | Promise<void>;
  /** The run closed the file: all its bytes were written. */
  close(): void | Promise<void>;
  /**
   * The file was left open because the run failed, was cancelled or a sink failed: its bytes are not all there. Optional;
   * a sink that holds a resource (a writable stream of a file) releases it here.
   */
  abort?(reason: unknown): void | Promise<void>;
}

/**
 * Where the output of a tool run goes as the run writes it. A stream without a sink is collected: its chunks are kept
 * as a Blob, never joined into one array (`ToolRunResult`). With a sink, that output is not kept: its Blob is empty,
 * and the file is not in `files`. Each call of a sink waits for the promise of the one before.
 *
 * A sink slower than the engine does not make the page hold the whole output: where the page can share memory with
 * the Engine worker (`SharedArrayBuffer`, which a cross-origin isolated page and Node have), the worker stops the engine
 * while more than 64 MiB of the bytes that it posted wait for the sinks, and goes on as the sinks
 * take them. Where there is no `SharedArrayBuffer` (a page that is not cross-origin isolated), nothing can wait: the
 * chunks queue on the page until the sink takes them, so a sink slower than the engine holds what the engine has
 * written and the sink has not taken (a platform-condition path; the output of a run to a fast sink, or one that a
 * Blob collects, is unaffected).
 */
export interface ToolSinks {
  /** A chunk of the tool's standard output (ABI v2 stream 11); the chunk is the sink's own. */
  readonly onStdout?: (chunk: Uint8Array) => void | Promise<void>;
  /** A chunk of the tool's standard error (stream 3), as the tool writes it. */
  readonly onStderr?: (chunk: Uint8Array) => void | Promise<void>;
  /**
   * The run opens a file named `name`: returns the sink of its bytes, or undefined to have the file collected as a Blob
   * in `files`. Called once for each file that the run opens, also when it opens one file twice (`-out X -logfile X` is
   * one file). The name is the argv's spelling and may hold `/` (`sub/x`, for a folder of the run's directory).
   *
   * `makeblastdb` writes into the run's directory, as the command line writes into its working directory (ABI v2
   * version 17, line NP WW4b1): its files (the database's, a `-logfile`) are offered here when the run has ended, each once,
   * in the byte order of their names (UTF-8), each with its bytes from offset 0 in order and then `close`. They are the
   * files that the run made or changed and that are there when it ends: a file that it removed or renamed away is not
   * offered, and a file of `directory` that it did not change is not offered either. The name is the file's path in the
   * directory (`db.pin`; `sub/dir/db.pin` for `-out sub/dir/db`, whose folders the run makes as the command line does).
   */
  readonly onFile?: (name: string) => ToolFileSink | undefined | Promise<ToolFileSink | undefined>;
}

/** A file that a tool run wrote and the host collected: its name as the argv wrote it, and its bytes as a Blob. */
export interface ToolFile {
  readonly name: string;
  readonly file: Blob;
}

export interface ToolRunOptions {
  /**
   * The files of the run's working directory that the words of the tool name (`-in`, `-ustat`, `-fa_list` and the files it
   * lists, `-ids`, `-exclude_ids`, `-conffile`, a BLAST database of `-infmt blastdb`; for `makeblastdb`, `-in`, `-mask_data`,
   * `-gi_mask` files, `-taxid_map`, a BLAST database of `-input_type blastdb`, `.ncbirc`), as `EngineRunRequest.files` gives
   * them to a search: the engine reads each in place, never copied whole. A `ToolFile` of an earlier run may be given here.
   * A tool's run may make, change, rename and remove files in this directory as on Linux (a `makeblastdb` that replaces a
   * database given here): what it writes is kept in the Engine worker's memory until the run ends, and the given Blobs are
   * never changed.
   */
  readonly directory?: readonly DirectoryFile[];
  readonly sinks?: ToolSinks;
  /** Aborting it cancels the run: the promise rejects with a RunCancelledError. */
  readonly signal?: AbortSignal;
}

/** What a tool run did. `status` is the command line's exit status; the engine returns for every status. */
export interface ToolRunResult {
  /** The exit status of the command line's run of the same words. */
  readonly status: number;
  /** What the tool wrote to its standard output; empty when `ToolSinks.onStdout` took it. */
  readonly stdout: Blob;
  /** What the tool wrote to its standard error; empty when `ToolSinks.onStderr` took it. */
  readonly stderr: Blob;
  /**
   * The files that the run wrote: those it sent as it wrote them (`windowmasker`'s `-out`, `-logfile`) in the order it
   * opened them, then those it made or changed in its directory (`makeblastdb`'s) in the byte order of their names; those
   * that `ToolSinks.onFile` took are not here.
   */
  readonly files: readonly ToolFile[];
}

/** The verdict of `validateTool`: the command line's error text (its standard error, the usage and the error line) for an argv that it refuses. */
export type ToolValidation = { readonly ok: true } | { readonly ok: false; readonly message: string };

/** The files that a check of a tool's words may read (`ToolRunOptions.directory`). */
export interface ToolValidateOptions {
  readonly directory?: readonly DirectoryFile[];
}

/** What an engine offers for the command line's tools besides the searches (`windowmasker`, `makeblastdb`). */
export interface ToolEngine {
  /**
   * Reads a tool's options, defaults, value lists, ranges and exclusions by running `<tool> -xmlhelp`
   * and parsing NCBI's PrintUsageXml output (`parseToolDescription`).
   */
  describeTool(tool: string): Promise<ToolDescription>;
  /**
   * Runs a tool as the command line runs it in the run's directory (`argv[0]` is the tool's name, `windowmasker` or
   * `makeblastdb`, then the command line's words), with its standard input empty. Streams the tool's standard output, standard error and files to
   * `options.sinks` as they are written (a file's bytes in frames of up to 1 MiB: a small file may arrive only when the run
   * ends), then resolves with the exit status (the command line's, for every status, a refusal
   * of an argument included: its text is on stderr). Rejects for a name that is not a tool that the engine runs, for a
   * cancel and for a failure of the engine or of a sink. One run at a time, a search's included.
   */
  runTool(argv: readonly string[], options?: ToolRunOptions): Promise<ToolRunResult>;
  /**
   * Checks the words of a tool as the command line checks them before the tool opens a file. It reads the run's registry
   * files (`-conffile`, `.ncbirc`, `windowmasker.ini`) as the command line does, so give it the same `directory` as
   * `runTool`: it opens them on the instance that the Engine worker has loaded (a threaded runtime that the searches use
   * stays as it is), and takes its turn there like a run (it waits behind a run in progress). It does not make the engine
   * busy: a run may start while a check is pending, and a check may be asked while a run is in progress. A check whose
   * instance stops rejects and ends the worker as a run's stop does: a run posted behind it on that worker fails (the same
   * error as a worker that ended), and the next run starts a new worker.
   */
  validateTool(argv: readonly string[], options?: ToolValidateOptions): Promise<ToolValidation>;
}
