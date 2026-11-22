// Where the output of a tool run goes on the main thread (D-WW3-6; line NP WW3a, 2026-11-22): the Engine worker posts the
// run's standard output, standard error and file frames as the tool writes them, and the Gateway hands each to this router,
// which gives it to the caller's sink, or, for a stream without one, keeps it as a Blob. Nothing is joined into one array:
// a stream without a sink is a list of chunks that become a Blob (a browser may keep a large Blob out of the page's memory),
// and a stream with one keeps nothing.
//
// The calls come in order and one at a time: every chunk and frame is queued behind the one before, and a sink's call is
// awaited before the next call of the router. The engine does not wait for the router: the worker stops it while more
// than the credit's limit of the bytes that it posted wait for the sinks (tool-credit.ts, D-WW3-11), where the page can
// share memory with the worker; without that, the events queue here until the sink has taken them.
//
// A file has POSIX semantics (D-WW3-9): a frame of bytes carries the offset at which they are written, and a truncate
// frame sets the file's length. The default collection keeps the file's exact bytes for writes at any offset.
import type { ToolFrame } from '../../domain/tool-frames.js';
import type { ToolFile, ToolFileSink, ToolSinks } from '../../ports/engine.js';
import { postedBytes } from './tool-credit.js';

/** Bytes that a collected stream holds as chunks before they are folded into its Blob. */
export const BLOB_FOLD_BYTES = 32 * 1024 * 1024;
/** The most zero bytes that a hole of a collected file puts in one piece. */
const ZEROS_PIECE = 1024 * 1024;
const ZEROS = new Uint8Array(ZEROS_PIECE);

/** What the router tells the page's side of the run's credit (tool-credit.ts, `ToolCreditPage`). */
export interface ToolCredit {
  /** `bytes` that the worker posted have been handled. */
  take(bytes: number): void;
}

/**
 * The bytes of a collected stream or file as a Blob. Appended chunks wait in a list and are folded into the Blob at
 * `BLOB_FOLD_BYTES` (`new Blob([blob, ...chunks])` shares the earlier Blob's data), so that a large output is not held as
 * thousands of arrays and is never joined into one array. A write before the end, a truncate and a write past the end
 * (zero bytes between) rebuild the Blob from slices of itself (`Blob.slice`), without reading its bytes into the page.
 */
class CollectedBlob {
  private blob = new Blob([]);
  private chunks: BlobPart[] = [];
  private pending = 0;
  private length = 0;

  /** The bytes at the end. */
  append(chunk: Uint8Array): void {
    this.chunks.push(chunk as BlobPart);
    this.pending += chunk.length;
    this.length += chunk.length;
    if (this.pending >= BLOB_FOLD_BYTES) this.fold();
  }

  /** `chunk` at `offset`: an append at the end, a hole of zero bytes before a write past the end, an overwrite before it. */
  write(chunk: Uint8Array, offset: number): void {
    if (offset > this.length) this.zeros(offset - this.length);
    if (offset >= this.length) {
      this.append(chunk);
      return;
    }
    this.fold();
    const end = offset + chunk.length;
    const parts: BlobPart[] = [this.blob.slice(0, offset), chunk as BlobPart];
    if (end < this.length) parts.push(this.blob.slice(end));
    this.blob = new Blob(parts);
    this.length = Math.max(this.length, end);
  }

  /** The length becomes `length`: cut, or with zero bytes added. */
  truncate(length: number): void {
    if (length > this.length) {
      this.zeros(length - this.length);
      return;
    }
    if (length === this.length) return;
    this.fold();
    this.blob = this.blob.slice(0, length);
    this.length = length;
  }

  finish(): Blob {
    this.fold();
    return this.blob;
  }

  /** `count` zero bytes at the end, in pieces of a megabyte that share one array. */
  private zeros(count: number): void {
    for (let left = count; left > 0; left -= ZEROS_PIECE) this.append(ZEROS.subarray(0, Math.min(left, ZEROS_PIECE)));
  }

  private fold(): void {
    if (this.chunks.length === 0) return;
    this.blob = new Blob([this.blob, ...this.chunks]);
    this.chunks = [];
    this.pending = 0;
  }
}

interface OpenFile {
  readonly name: string;
  readonly sink: ToolFileSink | undefined;
  readonly collected: CollectedBlob | undefined;
  ended: boolean;
}

/** What a finished run gave to no sink. */
export interface ToolCollected {
  readonly stdout: Blob;
  readonly stderr: Blob;
  readonly files: readonly ToolFile[];
}

export class ToolOutputRouter {
  private tail: Promise<void> = Promise.resolve();
  private failure: { readonly error: unknown } | undefined;
  private readonly stdoutBlob = new CollectedBlob();
  private readonly stderrBlob = new CollectedBlob();
  private readonly files = new Map<number, OpenFile>();

  /**
   * `failed` is told once, with the first error of a sink or of the frames (the Gateway cancels the run). `credit` takes the
   * bytes of each chunk and frame once the router has handled it: its sink call settled, it went into a Blob, or it was
   * dropped after a failure.
   */
  constructor(
    private readonly sinks: ToolSinks | undefined,
    private readonly failed: (error: unknown) => void,
    private readonly credit?: ToolCredit,
  ) {}

  /** A chunk of the tool's standard output; the router owns it from here. */
  stdout(chunk: Uint8Array): void {
    const sink = this.sinks?.onStdout;
    this.queue(() => (sink === undefined ? this.stdoutBlob.append(chunk) : sink(chunk)), chunk.length);
  }

  /** A chunk of the tool's standard error. */
  stderr(chunk: Uint8Array): void {
    const sink = this.sinks?.onStderr;
    this.queue(() => (sink === undefined ? this.stderrBlob.append(chunk) : sink(chunk)), chunk.length);
  }

  /** A frame of stream 12: the beginning, bytes, truncation or end of a file. */
  frame(frame: ToolFrame): void {
    this.queue(() => this.apply(frame), postedBytes(frame));
  }

  /**
   * Waits until every chunk and frame has been taken, and returns what went to no sink. Rejects with the first error of a sink
   * or of the frames, or when a file is still open (a run that ended in its status closed every file, ABI v2 §5).
   */
  async finish(): Promise<ToolCollected> {
    await this.tail;
    if (this.failure !== undefined) throw this.failure.error;
    const open = [...this.files.entries()].find(([, file]) => !file.ended);
    if (open !== undefined) throw new Error(`the tool run ended with its file ${open[0]} (${open[1].name}) still open`);
    const collected = [...this.files.entries()]
      .sort(([a], [b]) => a - b)
      .flatMap(([, file]) => (file.collected === undefined ? [] : [{ name: file.name, file: file.collected.finish() }]));
    return { stdout: this.stdoutBlob.finish(), stderr: this.stderrBlob.finish(), files: collected };
  }

  /**
   * The run ended without its status (a failure, a cancel): drops what is queued behind the call in progress, and tells each
   * sink of an open file (`ToolFileSink.abort`). A sink's own error here is not reported: the run's error is.
   */
  abort(reason: unknown): Promise<void> {
    this.failure ??= { error: reason };
    this.tail = this.tail.then(() => this.abortOpen(reason));
    return this.tail;
  }

  /** Queues `work` for `size` posted bytes, which the credit takes when the work is done or dropped. */
  private queue(work: () => void | Promise<void>, size: number): void {
    this.tail = this.tail.then(async () => {
      try {
        if (this.failure !== undefined) return;
        try {
          await work();
        } catch (error) {
          this.failure = { error };
          this.failed(error);
          await this.abortOpen(error);
        }
      } finally {
        this.credit?.take(size);
      }
    });
  }

  /** After a failure: the files that a sink holds open are told, once. */
  private async abortOpen(reason: unknown): Promise<void> {
    for (const file of this.files.values()) {
      if (file.ended || file.sink?.abort === undefined) continue;
      file.ended = true;
      try {
        await file.sink.abort(reason);
      } catch {
        // The first error is the one to report.
      }
    }
  }

  private async apply(frame: ToolFrame): Promise<void> {
    if (frame.kind === 'begin') {
      if (this.files.has(frame.number)) throw new Error(`the tool run began its file ${frame.number} twice`);
      const sink = await this.sinks?.onFile?.(frame.name);
      this.files.set(frame.number, {
        name: frame.name,
        sink,
        collected: sink === undefined ? new CollectedBlob() : undefined,
        ended: false,
      });
      return;
    }
    const file = this.files.get(frame.number);
    if (file === undefined) throw new Error(`the tool run sent a frame of its file ${frame.number}, which it did not begin`);
    if (file.ended) throw new Error(`the tool run sent a frame of its file ${frame.number} (${file.name}) after its end`);
    if (frame.kind === 'bytes') {
      if (file.sink === undefined) file.collected!.write(frame.bytes, frame.offset);
      else await file.sink.write(frame.bytes, frame.offset);
      return;
    }
    if (frame.kind === 'truncate') {
      if (file.sink === undefined) file.collected!.truncate(frame.length);
      else await file.sink.truncate(frame.length);
      return;
    }
    file.ended = true;
    await file.sink?.close();
  }
}
