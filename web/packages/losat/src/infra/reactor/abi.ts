// The ABI v2 binding (docs/web/abi_v2.md) over one reactor instance: memory for the
// arguments (§6), the error text of a failed call (§4), and the `losat_host.emit` stream
// receiver (§5). The Engine worker (searches) and the Data worker (index scan, input check,
// describe, validate, and one instance for each inflation) both use it.
//
// Pointers and lengths are unsigned 32-bit values (ABI v2 §6). JavaScript receives a Wasm
// i32 as a signed number, so a pointer of 2 GiB or more arrives negative: every pointer
// and length is converted with `>>> 0` before it is used as an offset into the memory.
import type { FastaReader, RecordKey, ShownIds } from '../../domain/dataset.js';
import { isPrintedKind } from '../../domain/run-files.js';
import { parseToolFrame, type ToolFrame } from '../../domain/tool-frames.js';
import { FileReadError, type StartResult, type ToolValidation, type ValidationResult } from '../../ports/engine.js';
import type { Compression } from '../../ports/inflate.js';
import { DIAGNOSTICS_STREAM } from '../../ports/run-output.js';

/** The version of ABI v2 that this host speaks (`losat_web2_abi_version`, ABI v2 §10). */
export const ABI_VERSION = 15;
/** Stream 2: the JSON response of describe, validate, register, register_inputs, scan_end, inflate_begin, start and file_end. */
export const RESPONSE_STREAM = 2;
/** Stream 5: the inflated bytes of inflate_chunk and inflate_end. */
export const INFLATED_STREAM = 5;
/**
 * Stream 8: the files of a run besides its outputs (the file of `-logfile`), a ZIP archive of stored
 * entries that the engine sends after the run, also when the run fails (ABI v2 §5; line RG E1w-3b).
 */
export const RUN_FILES_STREAM = 8;
/**
 * Stream 9: a run's signal that an output stream failed, whose received chunks the host discards,
 * and where its error is (ABI v2 §5, §7; line RG T7-W): every output stream that failed alone, or
 * whose search failed in a run of two searches, since ABI version 4 (line WQA's sweep, D7,
 * 2026-11-19 and 2026-11-20). Each `emit` is one signal; `discardSignal` reads it.
 */
export const DISCARD_STREAM = 9;
/**
 * Stream 10: a run's signal of where its messages go, the messages that the CLI writes to its standard error without
 * -logfile (ABI v2 §5, version 4; line WQA's sweep, D7, 2026-11-20): one signal per run that its start did not end,
 * before its outputs. `messagesSignal` reads it.
 */
export const MESSAGES_STREAM = 10;
/**
 * Stream 11: the text that the command line prints to its standard output in place of a run (its help, the XML of
 * `-xmlhelp`, or a version), for an argv whose start `prints` (ABI v2 §5, version 10; line WQA's D14, 2026-11-20).
 */
export const PRINTED_STREAM = 11;
/**
 * Stream 12: the files that a tool run writes (`losat_web2_tool`, ABI v2 §5, version 13; line NP WW3a, D-WW3-6, 2026-11-22),
 * each as it is written: every `emit` is one frame (domain/tool-frames.ts, `parseToolFrame`).
 */
export const TOOL_FILES_STREAM = 12;
/** The place "nowhere" in the signals of streams 9 and 10 (ABI v2 §5): TBLASTX's `-logfile ""` or `-logfile /dev/null`. */
export const PLACE_NOWHERE = 0xffff_ffff;
/**
 * The place "the error that ended the run" in the signals of stream 9 (ABI v2 §5, version 7; line WQA's sweep, D15,
 * 2026-11-20): an output that a run which failed does not keep, because the CLI's run of the argv does not write it.
 */
export const PLACE_RUN = 0xffff_fffe;

/**
 * A signal of stream 9 (ABI v2 §5): the output stream that failed, and the stream that holds its
 * error (3, the diagnostics; 8, the run's files: the file of -logfile), or `undefined` where the
 * run wrote it nowhere (`PLACE_NOWHERE`); `ended` where it is the error that ended the run (`PLACE_RUN`), which no stream
 * holds.
 */
export interface DiscardSignal {
  readonly stream: number;
  readonly errorStream: number | undefined;
  readonly ended: boolean;
}

/** Reads a signal of stream 9: its 8 bytes, the stream and the place of its error, each an unsigned 32-bit little-endian number. */
export function discardSignal(bytes: Uint8Array): DiscardSignal {
  if (bytes.length !== 8) throw new Error(`a signal of stream ${DISCARD_STREAM} has ${bytes.length} bytes, not 8`);
  const view = new DataView(bytes.buffer, bytes.byteOffset, 8);
  const place = view.getUint32(4, true);
  const ended = place === PLACE_RUN;
  return { stream: view.getUint32(0, true), errorStream: place === PLACE_NOWHERE || ended ? undefined : place, ended };
}

/** The output stream that a signal of stream 9 discards (`discardSignal`). */
export function discardedStream(bytes: Uint8Array): number {
  return discardSignal(bytes).stream;
}

/**
 * Reads the signal of stream 10: its 4 bytes, the place of the run's messages as an unsigned 32-bit little-endian
 * number: the stream that holds them (3, the diagnostics; 8, the run's files: the file of -logfile), or `undefined`
 * where the run writes them nowhere (`PLACE_NOWHERE`).
 */
export function messagesSignal(bytes: Uint8Array): number | undefined {
  if (bytes.length !== 4) throw new Error(`a signal of stream ${MESSAGES_STREAM} has ${bytes.length} bytes, not 4`);
  const place = new DataView(bytes.buffer, bytes.byteOffset, 4).getUint32(0, true);
  return place === PLACE_NOWHERE ? undefined : place;
}
export const ROLE_QUERY = 0;
export const ROLE_SUBJECT = 1;
/** The kind of a named input file (`file_begin`, ABI v2 §4): a WindowMasker unit counts file, BLASTN's `-window_masker_db`. */
export const FILE_WINDOW_MASKER_COUNTS = 1;
/** The kind of a named input file (`file_begin`, ABI v2 §4): a registry file of NCBI's application start, the file of `-conffile`. */
export const FILE_REGISTRY = 2;
/** Bytes of a named input file given to the engine at a time, through one buffer of engine memory (`loadFile`). */
export const FILE_CHUNK_BYTES = 8 * 1024 * 1024;

type Export = (...args: number[]) => number;

export interface AbiExports {
  readonly memory: WebAssembly.Memory;
  readonly losat_web2_abi_version: Export;
  readonly losat_web2_alloc: Export;
  readonly losat_web2_dealloc: Export;
  readonly losat_web2_last_error_ptr: Export;
  readonly losat_web2_last_error_len: Export;
  readonly losat_web2_describe: Export;
  readonly losat_web2_start: Export;
  readonly losat_web2_validate: Export;
  readonly losat_web2_register: Export;
  readonly losat_web2_register_inputs: Export;
  readonly losat_web2_release: Export;
  readonly losat_web2_scan_begin: Export;
  readonly losat_web2_scan_chunk: Export;
  readonly losat_web2_scan_end: Export;
  readonly losat_web2_inflate_begin: Export;
  readonly losat_web2_inflate_chunk: Export;
  readonly losat_web2_inflate_end: Export;
  readonly losat_web2_run: Export;
  readonly losat_web2_tool: Export;
  readonly losat_web2_file_begin: Export;
  readonly losat_web2_file_chunk: Export;
  readonly losat_web2_file_end: Export;
  readonly losat_web2_file_release: Export;
}

/** Receives the output streams of a run; `bytes` is valid only during the call. */
export type StreamSink = (stream: number, bytes: Uint8Array) => void;

/**
 * Receives the streams of a tool run as the engine writes them (ABI v2 §4 `losat_web2_tool`): the standard output
 * (stream 11), the standard error (stream 3) and the frames of the files (stream 12). `bytes` of a chunk or of a frame of
 * bytes is valid only during the call: the receiver copies it (`ownedFrame`) or posts it on.
 */
export interface ToolSink {
  stdout(bytes: Uint8Array): void;
  stderr(bytes: Uint8Array): void;
  file(frame: ToolFrame): void;
}

/** An export returned -1; the message is the engine's error text (ABI v2 §4). */
export class EngineCallError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'EngineCallError';
  }
}

/**
 * The host could not place an input in the instance's memory (`losat_web2_alloc` returned
 * 0). The engine did not read the input, so this is not the engine's verdict on it.
 */
export class EngineAllocationError extends EngineCallError {
  constructor(message: string) {
    super(message);
    this.name = 'EngineAllocationError';
  }
}

/**
 * A call into the instance threw (a trap, such as a panic, which aborts, or an out of
 * memory condition). The engine's state may be inconsistent, so the instance is not used
 * again.
 */
export class EngineStoppedError extends Error {
  constructor(message: string, options?: { cause?: unknown }) {
    super(message, options);
    this.name = 'EngineStoppedError';
  }
}

export interface RegisteredInput {
  readonly handle: number;
  readonly records: readonly RecordKey[];
  /**
   * The reader's error that ends a query's records (*register* `reader_error`, the CLI's text): a run
   * of the handle reads the records before it and ends with it in its query batch (line WQA's D10,
   * 2026-11-19).
   */
  readonly readerError?: string;
}

/**
 * The inputs of one run that `losat_web2_register_inputs` read together for its argv (a run with
 * `-parse_deflines`; ABI v2 §4, §9 *inputs*), and the IDs that its outputs show for their records.
 */
export interface RegisteredInputs {
  /** None for a run of `-in_pssm` or `-in_msa`, whose query is read from `/work` (*inputs* `query` null). */
  readonly query: RegisteredInput | undefined;
  /** None for a run of `-db`, whose subjects are the database's (*inputs* `subject` null). */
  readonly subject: RegisteredInput | undefined;
  readonly shownIds: ShownIds;
}

type RecordsResponse = {
  handle: number;
  records: Array<{ id: string; length: number; shown_id?: string | null }>;
  reader_error?: string;
};

/** A registered input of a *register* or *inputs* response. */
function registeredInput(input: RecordsResponse): RegisteredInput {
  return {
    handle: input.handle,
    records: Object.freeze(input.records.map(({ id, length }) => ({ id, length }))),
    ...(typeof input.reader_error === 'string' ? { readerError: input.reader_error } : {}),
  };
}

/**
 * The *file* JSON of `file_end` (ABI v2 §9). The contents of the file do not fail `file_end`: a
 * counts file that BLASTN cannot read has `stat_type` other than `obinary` and `binary`, or an
 * `error`, and the run that names it reports the CLI's own message; the start of `validate` and
 * `run` reads a registry file (kind 2) and reports it where the CLI does.
 */
export type FileInfo = CountsFileInfo | RegistryFileInfo;

interface FileInfoOfAnyKind {
  readonly handle: number;
  readonly name: string;
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes that the engine received. */
  readonly sha256: string;
}

/** The *file* JSON of a WindowMasker unit counts file (kind 1). */
export interface CountsFileInfo extends FileInfoOfAnyKind {
  readonly kind: 'windowmasker_counts';
  readonly stat_type: 'obinary' | 'binary' | 'ascii' | 'oascii' | 'unknown' | null;
  readonly unit_size: number | null;
  readonly error: string | null;
}

/** The *file* JSON of a registry file (kind 2). */
export interface RegistryFileInfo extends FileInfoOfAnyKind {
  readonly kind: 'registry';
}

/** Memory of the instance for the chunks of a file, taken once and used for every chunk (`chunkBuffer`). */
export interface ChunkBuffer {
  readonly ptr: number;
  readonly capacity: number;
}

/** An inflation that `inflate_begin` started, with its *inflate* response (ABI v2 §9). */
export interface BegunInflation {
  readonly handle: number;
  readonly compression: Compression;
  readonly plainName: string;
}

const encoder = new TextEncoder();

/** Bytes copied into engine memory that the caller holds (`ChunkBuffer`). */
interface Placed {
  readonly ptr: number;
  readonly len: number;
  readonly bytes: Uint8Array;
}

export class ReactorAbi {
  private sink: StreamSink | undefined;
  private responses: Uint8Array[] = [];
  private stopped: EngineStoppedError | undefined;

  /**
   * `exports.memory` is read again at every access: the buffer of a memory changes when
   * it grows (for the threaded module, in any thread). `output` returns what the engine
   * wrote to its standard error, for the message when the instance stops.
   */
  constructor(
    private readonly exports: AbiExports,
    private readonly output: () => string = () => '',
  ) {}

  /** The `losat_host.emit` import of the instance. */
  readonly emit = (stream: number, ptr: number, len: number): void => {
    const view = new Uint8Array(this.exports.memory.buffer, ptr >>> 0, len >>> 0);
    if (stream === RESPONSE_STREAM) this.responses.push(view.slice());
    else this.sink?.(stream, view);
  };

  /** The error that stopped the instance, if it stopped. */
  get stoppedBy(): EngineStoppedError | undefined {
    return this.stopped;
  }

  abiVersion(): number {
    return this.call(this.exports.losat_web2_abi_version, []);
  }

  /** Bytes of linear memory now. */
  memoryBytes(): number {
    return this.exports.memory.buffer.byteLength;
  }

  /** The *describe* JSON of a program (ABI v2 §9). */
  describe(program: string): unknown {
    this.check(this.call(this.exports.losat_web2_describe, [encoder.encode(program)]));
    return this.response();
  }

  /**
   * The *start* JSON of an argv (ABI v2 §4 `losat_web2_start`, §9): the registry file of its
   * `-conffile`, the file of its `-logfile`, the file of its `-export_search_strategy` (an engine
   * before line RG item F1w gives none), whether it is a dry run, whether its inputs are read
   * together (`registerInputs`), what it prints in place of a search (`prints`) and whether it searches the database of
   * `-db` (`database`); or the CLI's message for an argv that the engine cannot start (a program that ABI v2 does not run).
   */
  start(argv: readonly string[]): StartResult {
    const status = this.call(this.exports.losat_web2_start, [encoder.encode(argv.join('\0'))]);
    if (status < 0) return { ok: false, message: this.lastError() };
    const json = this.response() as {
      conffile: string | null;
      log_file: string | null;
      log_file_name?: string | null;
      strategy_file?: string | null;
      dry_run: boolean;
      register_inputs: boolean;
      prints: string | null;
      database: boolean;
      query_file?: string | null;
      out_pssm?: string | null;
      out_ascii_pssm?: string | null;
      save_each_pssm?: boolean;
    };
    if (json.prints !== null && !isPrintedKind(json.prints)) throw new EngineCallError(`the engine's start says that a run prints ${JSON.stringify(json.prints)}`);
    return {
      ok: true,
      start: {
        conffile: json.conffile,
        logFile: json.log_file,
        logFileName: json.log_file_name ?? null,
        strategyFile: json.strategy_file ?? null,
        dryRun: json.dry_run,
        registerInputs: json.register_inputs,
        prints: json.prints,
        database: json.database,
        queryFile: json.query_file ?? null,
        outPssm: json.out_pssm ?? null,
        outAsciiPssm: json.out_ascii_pssm ?? null,
        saveEachPssm: json.save_each_pssm === true,
      },
    };
  }

  /**
   * The verdict on an argv: the CLI's error text, or whether its run reads the inputs with NCBI's
   * data loaders on (the *validation* JSON, ABI v2 §9, version 8), which `register` and `scanBegin`
   * then take.
   */
  validate(argv: readonly string[]): ValidationResult {
    const status = this.call(this.exports.losat_web2_validate, [encoder.encode(argv.join('\0'))]);
    // The *validation* comes before the option checks, so a refused argv may have it too (line WQA's D13 round 3, L1):
    // it is read before the error, whose call drops the responses.
    const known = status < 0 && this.responses.length === 0 ? undefined : (this.response() as { data_loaders: boolean }).data_loaders;
    if (status < 0) return { ok: false, message: this.lastError(), ...(known === undefined ? {} : { dataLoaders: known }) };
    return { ok: true, dataLoaders: known! };
  }

  /**
   * The verdict on the argv of a tool run (`losat_web2_validate` with a tool's name as `argv[0]`, ABI v2 §4, version 13): the
   * command line's error text, or valid. The call emits no response for it, and none is read.
   */
  validateTool(argv: readonly string[]): ToolValidation {
    const status = this.call(this.exports.losat_web2_validate, [encoder.encode(argv.join('\0'))]);
    return status < 0 ? { ok: false, message: this.lastError() } : { ok: true };
  }

  /** Reads `bytes` with the engine's reader of `program`'s input of `role`, NCBI's data loaders on or off (ABI v2 §4, version 8). */
  register(program: string, role: number, dataLoaders: boolean, bytes: Uint8Array): RegisteredInput {
    const handle = this.check(this.call(this.exports.losat_web2_register, [encoder.encode(program), role, dataLoaders ? 1 : 0, bytes]));
    try {
      return { ...registeredInput(this.response() as RecordsResponse), handle };
    } catch (error) {
      // A handle whose response cannot be read is not used: give its copy of the input back.
      if (this.stopped === undefined) this.release(handle);
      throw error;
    }
  }

  /**
   * Reads the query and the subject of the run of `argv` together, as the CLI reads its `-query` and
   * `-subject` for that argv (ABI v2 §4 `losat_web2_register_inputs`): the inputs of a run whose
   * *start* JSON says `register_inputs` (`-parse_deflines`, line WQA's sweep, D7, 2026-11-19). Only
   * the run of `argv` takes the two handles; the error is the CLI's message for these inputs. A run of
   * `-db` takes no subject: its `subject` is empty, and the query is read alone (line WQA D12a).
   */
  registerInputs(argv: readonly string[], query: Uint8Array, subject: Uint8Array): RegisteredInputs {
    this.check(this.call(this.exports.losat_web2_register_inputs, [encoder.encode(argv.join('\0')), query, subject]));
    const response = this.response() as { query: RecordsResponse | null; subject: RecordsResponse | null };
    const handles = [response?.query?.handle, response?.subject?.handle].filter((handle) => typeof handle === 'number');
    try {
      // A record past the reader's error that ends the query's records has no shown ID (`null`): no
      // run reads it, and its name is its ID (line WQA's D10, 2026-11-19).
      const shown = (input: RecordsResponse) =>
        Object.freeze(
          input.records.map(({ id, shown_id }) => {
            if (shown_id === null && input.reader_error !== undefined) return id;
            if (typeof shown_id !== 'string') throw new Error('a record of register_inputs has no shown_id');
            return shown_id;
          }),
        );
      const queries = response.query;
      const subjects = response.subject;
      return {
        query: queries === null ? undefined : registeredInput(queries),
        subject: subjects === null ? undefined : registeredInput(subjects),
        shownIds: Object.freeze({
          query: queries === null ? Object.freeze([]) : shown(queries),
          subject: subjects === null ? Object.freeze([]) : shown(subjects),
        }),
      };
    } catch (error) {
      // Handles whose response cannot be read are not used: give their copies of the inputs back.
      if (this.stopped === undefined) for (const handle of handles) this.release(handle);
      throw error;
    }
  }

  release(handle: number): void {
    this.check(this.call(this.exports.losat_web2_release, [handle]));
  }

  /** Starts an index scan with the engine's reader `reader` (ABI v2 §4, version 8). */
  scanBegin(reader: FastaReader): number {
    return this.check(
      this.call(this.exports.losat_web2_scan_begin, [
        encoder.encode(reader.program),
        reader.role === 'query' ? ROLE_QUERY : ROLE_SUBJECT,
        reader.dataLoaders ? 1 : 0,
      ]),
    );
  }

  scanChunk(scanner: number, bytes: Uint8Array): void {
    this.check(this.call(this.exports.losat_web2_scan_chunk, [scanner, bytes]));
  }

  /** Ends a scan and returns its *scan* JSON; rejects with the parser's error. */
  scanEnd(scanner: number): unknown {
    this.check(this.call(this.exports.losat_web2_scan_end, [scanner]));
    return this.response();
  }

  /** Starts a named input file of `size` bytes; returns its handle. `name` is the value the argv gives its option. */
  fileBegin(kind: number, name: string, size: number): number {
    // The engine's size is a 32-bit unsigned number: a larger one would wrap.
    if (!Number.isSafeInteger(size) || size < 0 || size > 0xffff_ffff) {
      throw new EngineCallError(`the file ${name} (${size} bytes) is larger than the engine can hold`);
    }
    return this.check(this.call(this.exports.losat_web2_file_begin, [kind, encoder.encode(name), size]));
  }

  /**
   * Memory for the chunks of files, taken once: every chunk is copied into it and handed to
   * the engine, so a file of many chunks allocates once. Give it back with `freeChunkBuffer`.
   */
  chunkBuffer(capacity: number): ChunkBuffer {
    const ptr = this.trapped(() => this.exports.losat_web2_alloc(capacity)) >>> 0;
    if (ptr === 0) throw new EngineAllocationError(`the engine could not allocate ${capacity} bytes for the chunks of a file`);
    return { ptr, capacity };
  }

  freeChunkBuffer(buffer: ChunkBuffer): void {
    if (this.stopped === undefined) this.trapped(() => this.exports.losat_web2_dealloc(buffer.ptr, buffer.capacity));
  }

  /** Calls the engine outside `call`: a trap stops the instance, as one in `call` does. */
  private trapped<T>(body: () => T): T {
    if (this.stopped !== undefined) throw this.stopped;
    try {
      return body();
    } catch (error) {
      this.stopped = new EngineStoppedError(stopMessage(error, this.output()), { cause: error });
      throw this.stopped;
    }
  }

  /** The next bytes of a file, in order, at most `buffer.capacity` of them. */
  fileChunk(file: number, bytes: Uint8Array, buffer: ChunkBuffer): void {
    if (bytes.length > buffer.capacity) throw new Error(`a chunk of ${bytes.length} bytes does not fit the chunk buffer of ${buffer.capacity}`);
    this.check(this.call(this.exports.losat_web2_file_chunk, [file, { ptr: buffer.ptr, len: bytes.length, bytes }]));
  }

  /** Ends a file; returns its *file* JSON. A short file fails with the engine's message. */
  fileEnd(file: number): FileInfo {
    this.check(this.call(this.exports.losat_web2_file_end, [file]));
    return this.response() as FileInfo;
  }

  /** Forgets a file, registered or loading, and frees its memory. */
  fileRelease(file: number): void {
    this.check(this.call(this.exports.losat_web2_file_release, [file]));
  }

  /**
   * Starts the inflation of a file named `name` (ABI v2 §4): the engine's name rule gives the
   * compression and the plain name; undefined when the name is plain.
   */
  inflateBegin(name: string): BegunInflation | undefined {
    const handle = this.check(this.call(this.exports.losat_web2_inflate_begin, [encoder.encode(name)]));
    if (handle === 0) return undefined;
    try {
      const response = this.response() as { compression: Compression; plain_name: string };
      return { handle, compression: response.compression, plainName: response.plain_name };
    } catch (error) {
      // An inflation whose response cannot be read is not used: free it.
      if (this.stopped === undefined) {
        try {
          this.inflateEnd(handle, () => undefined);
        } catch {
          // The response's error is the one to report.
        }
      }
      throw error;
    }
  }

  /** Inflates the next bytes of a file; the inflated bytes go to `sink` during the call. */
  inflateChunk(inflation: number, bytes: Uint8Array, sink: (bytes: Uint8Array) => void): void {
    this.inflating(sink, () => this.check(this.call(this.exports.losat_web2_inflate_chunk, [inflation, bytes])));
  }

  /** Ends a file: the rest of its inflated bytes go to `sink`; throws the engine's sentence for data cut short. */
  inflateEnd(inflation: number, sink: (bytes: Uint8Array) => void): void {
    this.inflating(sink, () => this.check(this.call(this.exports.losat_web2_inflate_end, [inflation])));
  }

  /**
   * Runs one search (ABI v2 §4); the output streams go to `sink` during the call. A run whose start `prints` reads no
   * input: its handles are 0.
   */
  run(argv: readonly string[], query: number, subject: number, sink: StreamSink): void {
    this.sink = sink;
    try {
      this.check(this.call(this.exports.losat_web2_run, [encoder.encode(argv.join('\0')), query, subject]));
    } finally {
      this.sink = undefined;
    }
  }

  /**
   * Runs a tool (ABI v2 §4 `losat_web2_tool`): its standard output, standard error and the frames of its files go to `sink`
   * during the call, as the tool writes them, and nothing of them is kept here; returns the exit status that the run's
   * *status* response gives (the command line's, for any status). A name that is not a tool that the engine runs throws
   * the engine's message.
   */
  tool(argv: readonly string[], sink: ToolSink): number {
    this.sink = (stream, bytes) => {
      if (stream === PRINTED_STREAM) sink.stdout(bytes);
      else if (stream === DIAGNOSTICS_STREAM) sink.stderr(bytes);
      else if (stream === TOOL_FILES_STREAM) sink.file(parseToolFrame(bytes));
    };
    try {
      this.check(this.call(this.exports.losat_web2_tool, [encoder.encode(argv.join('\0'))]));
    } finally {
      this.sink = undefined;
    }
    const { status } = this.response() as { status?: unknown };
    if (typeof status !== 'number' || !Number.isInteger(status)) throw new EngineCallError('the engine sent no exit status of the tool run');
    return status;
  }

  private inflating(sink: (bytes: Uint8Array) => void, body: () => void): void {
    this.sink = (stream, bytes) => {
      if (stream === INFLATED_STREAM) sink(bytes);
    };
    try {
      body();
    } finally {
      this.sink = undefined;
    }
  }

  private check(status: number): number {
    if (status < 0) throw new EngineCallError(this.lastError());
    return status;
  }

  /**
   * The stream 2 response of the last call. The adapter sends every stream in chunks of
   * at most 1 MiB (ABI v2 §5), so a large response (the scan of thousands of records)
   * arrives in several `emit` calls.
   */
  private response(): unknown {
    const chunks = this.responses;
    this.responses = [];
    if (chunks.length === 0) throw new EngineCallError('the engine sent no response');
    const bytes = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
    let offset = 0;
    for (const chunk of chunks) {
      bytes.set(chunk, offset);
      offset += chunk.length;
    }
    return JSON.parse(new TextDecoder().decode(bytes)) as unknown;
  }

  private lastError(): string {
    const ptr = this.call(this.exports.losat_web2_last_error_ptr, []) >>> 0;
    const len = this.call(this.exports.losat_web2_last_error_len, []) >>> 0;
    // A view of shared memory cannot be decoded directly; decode a copy.
    return new TextDecoder().decode(new Uint8Array(this.exports.memory.buffer, ptr, len).slice());
  }

  /**
   * Calls an export; each byte array becomes a (ptr, len) pair in engine memory, and a `Placed`
   * argument (bytes copied into memory the caller holds) a (ptr, len) pair as it is.
   */
  private call(fn: Export, args: ReadonlyArray<number | Uint8Array | Placed>): number {
    if (this.stopped !== undefined) throw this.stopped;
    // A response belongs to the call that sends it.
    this.responses = [];
    const allocations: Array<readonly [number, number]> = [];
    try {
      const flat: number[] = [];
      for (const arg of args) {
        if (typeof arg === 'number') {
          flat.push(arg);
          continue;
        }
        if (!(arg instanceof Uint8Array)) {
          // The memory view is taken again here: the buffer of a memory changes when it grows.
          new Uint8Array(this.exports.memory.buffer, arg.ptr, arg.len).set(arg.bytes);
          flat.push(arg.ptr, arg.len);
          continue;
        }
        const ptr = this.exports.losat_web2_alloc(arg.length) >>> 0;
        if (ptr === 0) throw new EngineAllocationError(`the engine could not allocate ${arg.length} bytes for an input`);
        allocations.push([ptr, arg.length]);
        new Uint8Array(this.exports.memory.buffer, ptr, arg.length).set(arg);
        flat.push(ptr, arg.length);
      }
      const result = fn(...flat);
      for (const [ptr, len] of allocations.splice(0)) this.exports.losat_web2_dealloc(ptr, len);
      return result;
    } catch (error) {
      if (error instanceof EngineCallError) {
        for (const [ptr, len] of allocations) this.exports.losat_web2_dealloc(ptr, len);
        throw error;
      }
      this.stopped = new EngineStoppedError(stopMessage(error, this.output()), { cause: error });
      throw this.stopped;
    }
  }
}

/** The calls of `loadFile` into an instance. */
export type FileLoader = Pick<ReactorAbi, 'stoppedBy' | 'fileBegin' | 'chunkBuffer' | 'fileChunk' | 'fileEnd' | 'fileRelease' | 'freeChunkBuffer'>;

/**
 * Gives the engine a named input file (ABI v2 §4: a WindowMasker counts file, kind 1, or a registry
 * file, kind 2) and returns its *file* JSON: `file_begin` with the file's size, its bytes in chunks of
 * at most `FILE_CHUNK_BYTES` read from the Blob one at a time and given through one buffer of engine
 * memory taken once, `file_end`. A file that the engine refuses (a chunk or its end), or that cannot
 * be read (`FileReadError`, with the browser's reason), is released, unless the instance stopped, so
 * that its name is free again; the buffer is given back in every case. `call` wraps each call into
 * the instance (the Engine worker counts them as uses of the threaded instance). The reads wait:
 * other calls of the instance may run between the chunks, so a file that must meet no other call
 * goes to an instance of its own.
 */
export async function loadFile(
  abi: FileLoader,
  kind: number,
  name: string,
  file: Blob,
  call: <T>(body: () => T) => T = (body) => body(),
): Promise<FileInfo> {
  const size = file.size;
  const handle = call(() => abi.fileBegin(kind, name, size));
  let buffer: ChunkBuffer | undefined;
  try {
    buffer = call(() => abi.chunkBuffer(Math.max(1, Math.min(FILE_CHUNK_BYTES, size))));
    for (let offset = 0; offset < size; offset += FILE_CHUNK_BYTES) {
      let bytes: Uint8Array;
      try {
        bytes = new Uint8Array(await file.slice(offset, Math.min(size, offset + FILE_CHUNK_BYTES)).arrayBuffer());
      } catch (error) {
        throw new FileReadError(error instanceof Error ? error.message : String(error), { cause: error });
      }
      call(() => abi.fileChunk(handle, bytes, buffer!));
    }
    return call(() => abi.fileEnd(handle));
  } catch (error) {
    // The engine drops a file whose chunk or end it refused; one that it still holds (the file could
    // not be read) is released. Not after a stop, nor when the caller's instance is ending.
    if (abi.stoppedBy === undefined) {
      try {
        call(() => abi.fileRelease(handle));
      } catch {
        // Already dropped by the engine, or the instance is ending.
      }
    }
    throw error;
  } finally {
    if (buffer !== undefined) {
      const taken = buffer;
      try {
        call(() => abi.freeChunkBuffer(taken));
      } catch {
        // The instance is ending: its memory goes with it.
      }
    }
  }
}

/** The message of an instance that stopped: the trap and the last line the engine wrote. */
export function stopMessage(error: unknown, output: string): string {
  const trap = error instanceof Error ? error.message : String(error);
  const lines = output.split('\n').map((line) => line.trim()).filter((line) => line !== '');
  const last = lines[lines.length - 1];
  const memory = /memory allocation of \d+ bytes failed/.test(output) ? ' (it ran out of memory)' : '';
  return `The engine stopped${memory}: ${last === undefined ? trap : `${last} (${trap})`}`;
}
