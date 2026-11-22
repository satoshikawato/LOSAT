// Where the output of a tool run goes on the main thread (D-WW3-6; line NP WW3a, 2026-11-22): the Engine worker posts the
// run's standard output, standard error and file frames as the tool writes them, and the Gateway hands each to this router,
// which gives it to the caller's sink, or, for a stream without one, keeps it as a Blob. Nothing is joined into one array:
// a stream without a sink is a list of chunks that become a Blob (a browser may keep a large Blob out of the page's memory),
// and a stream with one keeps nothing.
//
// The calls come in order and one at a time: every chunk and frame is queued behind the one before, and a sink's call is
// awaited before the next call of the router. The worker does not wait for a slow sink (the engine's call is synchronous),
// so the events queue here until the sink has taken them.
import type { ToolFrame } from '../../domain/tool-frames.js';
import type { ToolFile, ToolFileSink, ToolSinks } from '../../ports/engine.js';

/** Bytes that a collected stream holds as chunks before they are folded into its Blob. */
export const BLOB_FOLD_BYTES = 32 * 1024 * 1024;

/**
 * The chunks of a collected stream as a Blob. Chunks wait in a list and are folded into the Blob at `BLOB_FOLD_BYTES`
 * (`new Blob([blob, ...chunks])` shares the earlier Blob's data), so that a large output is not held as thousands of arrays.
 */
class BlobBuilder {
  private blob = new Blob([]);
  private chunks: Uint8Array[] = [];
  private pending = 0;

  add(chunk: Uint8Array): void {
    this.chunks.push(chunk);
    this.pending += chunk.length;
    if (this.pending >= BLOB_FOLD_BYTES) this.fold();
  }

  finish(): Blob {
    this.fold();
    return this.blob;
  }

  private fold(): void {
    if (this.chunks.length === 0) return;
    this.blob = new Blob([this.blob, ...(this.chunks as BlobPart[])]);
    this.chunks = [];
    this.pending = 0;
  }
}

interface OpenFile {
  readonly name: string;
  readonly sink: ToolFileSink | undefined;
  readonly collected: BlobBuilder | undefined;
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
  private readonly stdoutBlob = new BlobBuilder();
  private readonly stderrBlob = new BlobBuilder();
  private readonly files = new Map<number, OpenFile>();

  /** `failed` is told once, with the first error of a sink or of the frames (the Gateway cancels the run). */
  constructor(
    private readonly sinks: ToolSinks | undefined,
    private readonly failed: (error: unknown) => void,
  ) {}

  /** A chunk of the tool's standard output; the router owns it from here. */
  stdout(chunk: Uint8Array): void {
    const sink = this.sinks?.onStdout;
    this.queue(() => (sink === undefined ? this.stdoutBlob.add(chunk) : sink(chunk)));
  }

  /** A chunk of the tool's standard error. */
  stderr(chunk: Uint8Array): void {
    const sink = this.sinks?.onStderr;
    this.queue(() => (sink === undefined ? this.stderrBlob.add(chunk) : sink(chunk)));
  }

  /** A frame of stream 12: the beginning, bytes or end of a file. */
  frame(frame: ToolFrame): void {
    this.queue(() => this.apply(frame));
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

  private queue(work: () => void | Promise<void>): void {
    this.tail = this.tail.then(async () => {
      if (this.failure !== undefined) return;
      try {
        await work();
      } catch (error) {
        this.failure = { error };
        this.failed(error);
        await this.abortOpen(error);
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
        collected: sink === undefined ? new BlobBuilder() : undefined,
        ended: false,
      });
      return;
    }
    const file = this.files.get(frame.number);
    if (file === undefined) throw new Error(`the tool run sent a frame of its file ${frame.number}, which it did not begin`);
    if (file.ended) throw new Error(`the tool run sent a frame of its file ${frame.number} (${file.name}) after its end`);
    if (frame.kind === 'bytes') {
      if (file.sink === undefined) file.collected!.add(frame.bytes);
      else await file.sink.write(frame.bytes);
      return;
    }
    file.ended = true;
    await file.sink?.close();
  }
}
