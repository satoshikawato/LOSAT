// The ABI v2 binding (docs/web/abi_v2.md) over one reactor instance: memory for the
// arguments (§6), the error text of a failed call (§4), and the `losat_host.emit` stream
// receiver (§5). The Engine worker (searches) and the Data worker (index scan, describe,
// validate) both use it.
//
// Pointers and lengths are unsigned 32-bit values (ABI v2 §6). JavaScript receives a Wasm
// i32 as a signed number, so a pointer of 2 GiB or more arrives negative: every pointer
// and length is converted with `>>> 0` before it is used as an offset into the memory.
import type { RecordKey } from '../../domain/dataset';

/** Stream 2: the JSON response of describe, register and scan_end. */
export const RESPONSE_STREAM = 2;
export const ROLE_QUERY = 0;
export const ROLE_SUBJECT = 1;
/** The kind of a named input file (`file_begin`, ABI v2 §4): a WindowMasker unit counts file, BLASTN's `-window_masker_db`. */
export const FILE_WINDOW_MASKER_COUNTS = 1;

type Export = (...args: number[]) => number;

export interface AbiExports {
  readonly memory: WebAssembly.Memory;
  readonly losat_web2_abi_version: Export;
  readonly losat_web2_alloc: Export;
  readonly losat_web2_dealloc: Export;
  readonly losat_web2_last_error_ptr: Export;
  readonly losat_web2_last_error_len: Export;
  readonly losat_web2_describe: Export;
  readonly losat_web2_validate: Export;
  readonly losat_web2_register: Export;
  readonly losat_web2_release: Export;
  readonly losat_web2_scan_begin: Export;
  readonly losat_web2_scan_chunk: Export;
  readonly losat_web2_scan_end: Export;
  readonly losat_web2_run: Export;
  readonly losat_web2_file_begin: Export;
  readonly losat_web2_file_chunk: Export;
  readonly losat_web2_file_end: Export;
  readonly losat_web2_file_release: Export;
}

/** Receives the output streams of a run; `bytes` is valid only during the call. */
export type StreamSink = (stream: number, bytes: Uint8Array) => void;

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
}

/**
 * The *file* JSON of `file_end` (ABI v2 §9). The contents of the file do not fail `file_end`: a
 * counts file that BLASTN cannot read has `stat_type` other than `obinary` and `binary`, or an
 * `error`, and the run that names it reports the CLI's own message.
 */
export interface FileInfo {
  readonly handle: number;
  readonly name: string;
  readonly kind: 'windowmasker_counts';
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes that the engine received. */
  readonly sha256: string;
  readonly stat_type: 'obinary' | 'binary' | 'ascii' | 'oascii' | 'unknown' | null;
  readonly unit_size: number | null;
  readonly error: string | null;
}

/** Memory of the instance for the chunks of a file, taken once and used for every chunk (`chunkBuffer`). */
export interface ChunkBuffer {
  readonly ptr: number;
  readonly capacity: number;
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

  /** Undefined if the argv is valid, else the CLI's error text. */
  validate(argv: readonly string[]): string | undefined {
    const status = this.call(this.exports.losat_web2_validate, [encoder.encode(argv.join('\0'))]);
    return status < 0 ? this.lastError() : undefined;
  }

  register(program: string, role: number, bytes: Uint8Array): RegisteredInput {
    const handle = this.check(this.call(this.exports.losat_web2_register, [encoder.encode(program), role, bytes]));
    try {
      const response = this.response() as { handle: number; records: Array<{ id: string; length: number }> };
      return { handle, records: Object.freeze(response.records.map(({ id, length }) => ({ id, length }))) };
    } catch (error) {
      // A handle whose response cannot be read is not used: give its copy of the input back.
      if (this.stopped === undefined) this.release(handle);
      throw error;
    }
  }

  release(handle: number): void {
    this.check(this.call(this.exports.losat_web2_release, [handle]));
  }

  scanBegin(parser: number): number {
    return this.check(this.call(this.exports.losat_web2_scan_begin, [parser]));
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

  /** Runs one search (ABI v2 §4); the output streams go to `sink` during the call. */
  run(argv: readonly string[], query: number, subject: number, sink: StreamSink): void {
    this.sink = sink;
    try {
      this.check(this.call(this.exports.losat_web2_run, [encoder.encode(argv.join('\0')), query, subject]));
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

/** The message of an instance that stopped: the trap and the last line the engine wrote. */
export function stopMessage(error: unknown, output: string): string {
  const trap = error instanceof Error ? error.message : String(error);
  const lines = output.split('\n').map((line) => line.trim()).filter((line) => line !== '');
  const last = lines[lines.length - 1];
  const memory = /memory allocation of \d+ bytes failed/.test(output) ? ' (it ran out of memory)' : '';
  return `The engine stopped${memory}: ${last === undefined ? trap : `${last} (${trap})`}`;
}
