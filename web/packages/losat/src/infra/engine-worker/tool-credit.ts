// Backpressure of a tool run (D-WW3-11; line NP WW3a, 2026-11-23): the Engine worker posts the tool's standard output,
// standard error and file frames as the engine emits them, and the engine's call is synchronous, so a sink that takes
// the bytes slower than the engine writes them would make the page hold what has not been taken (a genome masked as
// FASTA is 3.1 GB). Where the page can share memory with the worker (a `SharedArrayBuffer`: a cross-origin isolated
// page, Node), the two count the bytes that were posted and not yet taken in a shared cell; the worker adds each
// chunk's bytes before it posts the chunk, and, while the cell is above the limit, waits (`Atomics.wait`, in slices) in
// the engine's call until the page has taken enough: the router (tool-output.ts) takes a chunk's bytes when its sink
// call settled, or when it went into a Blob. A run that ends, fails or is cancelled closes the cell, so that a worker
// that waits goes on (its instance ends with it). Where there is no `SharedArrayBuffer`, nothing waits: the chunks
// queue on the page (ports/engine.ts, `ToolSinks`).
import type { ToolFrame } from '../../domain/tool-frames.js';

/** Bytes that the worker has posted and the page has not taken, above which the worker waits: 64 MiB. */
export const TOOL_CREDIT_LIMIT = 64 * 1024 * 1024;
/** The longest that the worker waits at a time before it looks at the counters again, in milliseconds. */
export const TOOL_CREDIT_SLICE_MS = 100;

/** The cells of the shared buffer (an Int32Array): the bytes posted and not taken, and the flag that the run is over. */
const BACKLOG = 0;
const CLOSED = 1;
/** The bytes of the shared buffer. */
const CELL_BYTES = 8;

/**
 * The shared buffer for one tool run, or undefined where the page cannot share memory with its workers (no
 * `SharedArrayBuffer`, or a page that is not cross-origin isolated): nothing waits then.
 */
export function createToolCredit(): SharedArrayBuffer | undefined {
  if (typeof SharedArrayBuffer !== 'function' || globalThis.crossOriginIsolated === false) return undefined;
  return new SharedArrayBuffer(CELL_BYTES);
}

/** The bytes of a frame that count toward the backlog: those of a frame of bytes (the others are a few bytes). */
export function postedBytes(frame: ToolFrame): number {
  return frame.kind === 'bytes' ? frame.bytes.length : 0;
}

/** The page's side: it takes bytes as the sinks settle, and closes the cell when the run is over. */
export class ToolCreditPage {
  private readonly cells: Int32Array;

  constructor(buffer: SharedArrayBuffer) {
    this.cells = new Int32Array(buffer);
  }

  /** `bytes` that were posted have been handled (a sink took them, or they went into a Blob). */
  take(bytes: number): void {
    if (bytes === 0) return;
    Atomics.sub(this.cells, BACKLOG, bytes);
    Atomics.notify(this.cells, BACKLOG);
  }

  /** The run is over: a worker that waits goes on, and none waits again. */
  close(): void {
    Atomics.store(this.cells, CLOSED, 1);
    Atomics.notify(this.cells, BACKLOG);
  }
}

/** The worker's side: it counts the bytes that it posts, and waits in the engine's call while the backlog is above the limit. */
export class ToolCreditWorker {
  private readonly cells: Int32Array;
  private readonly limit: number;
  private readonly sliceMs: number;

  constructor(buffer: SharedArrayBuffer, limit: number = TOOL_CREDIT_LIMIT, sliceMs: number = TOOL_CREDIT_SLICE_MS) {
    this.cells = new Int32Array(buffer);
    this.limit = limit;
    this.sliceMs = sliceMs;
  }

  /** Counts `bytes` that are about to be posted (before the post, so that the page never takes what was not counted). */
  add(bytes: number): void {
    if (bytes !== 0) Atomics.add(this.cells, BACKLOG, bytes);
  }

  /** Waits until the backlog is at most the limit, or the run is over. */
  wait(): void {
    for (;;) {
      const backlog = Atomics.load(this.cells, BACKLOG);
      if (backlog <= this.limit || Atomics.load(this.cells, CLOSED) !== 0) return;
      // Wakes at once when the page takes bytes (`notify`), and in any case after a slice, to look at the flag again.
      Atomics.wait(this.cells, BACKLOG, backlog, this.sliceMs);
    }
  }
}
