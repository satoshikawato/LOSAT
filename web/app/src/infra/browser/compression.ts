// The Compression port with the browser's CompressionStream and DecompressionStream ('gzip';
// Chromium 80, Firefox 113, Safari 16.4), which Node also has for the unit tests.
import type { Compression } from '../../ports/compression';
import type { ExportSink } from '../../ports/download';
import { concatBytes } from '../bytes';
import { nextTask } from './next-task';

/**
 * Compressed bytes are handed on in blocks of about this size: the stream's own chunks are small
 * (tens of kB), and each block that a download's sink takes lets the page paint once.
 */
const OUT_BYTES = 1024 * 1024;
/**
 * Bytes handed to the compressor at once. It compresses on the page's thread, a whole write at a
 * time: a write of 8 MiB held the page for 100-450 ms (fix round 2); a piece of 256 KiB takes a few
 * milliseconds, and the page may draw between two pieces.
 */
const IN_BYTES = 256 * 1024;

export const browserCompression: Compression = {
  gzip(out) {
    return new GzipSink(out);
  },
  async *gunzip(file) {
    const reader = file.stream().pipeThrough(new DecompressionStream('gzip')).getReader();
    let finished = false;
    try {
      for (;;) {
        const { done, value } = await reader.read();
        if (done) {
          finished = true;
          return;
        }
        yield value;
      }
    } finally {
      // An iteration that stopped early (a refused file) stops reading it.
      if (!finished) await reader.cancel().catch(() => undefined);
    }
  },
};

class GzipSink implements ExportSink {
  private readonly writer: WritableStreamDefaultWriter<BufferSource>;
  /** Reads the compressed blocks and hands them to `out`, in order, while blocks are written. */
  private readonly pumping: Promise<void>;
  private failure: unknown;
  private state: 'open' | 'closed' | 'aborted' = 'open';

  constructor(out: (bytes: Uint8Array) => Promise<void>) {
    const stream = new CompressionStream('gzip');
    this.writer = stream.writable.getWriter();
    const reader = stream.readable.getReader();
    this.pumping = (async () => {
      let parts: Uint8Array[] = [];
      let held = 0;
      for (;;) {
        const { done, value } = await reader.read();
        if (!done) {
          parts.push(value);
          held += value.length;
        }
        if ((done || held >= OUT_BYTES) && held > 0) {
          const block = concatBytes(parts);
          parts = [];
          held = 0;
          await out(block);
        }
        if (done) return;
      }
    })();
    // A failure of `out` stops the writes that wait for the reader (backpressure).
    this.pumping.catch((error: unknown) => {
      this.failure = error;
      void this.writer.abort(error).catch(() => undefined);
      void reader.cancel(error).catch(() => undefined);
    });
  }

  async write(bytes: Uint8Array): Promise<void> {
    if (this.state !== 'open') throw new Error('the compressed file is no longer open');
    if (this.failure !== undefined) throw this.failure;
    for (let at = 0; at < bytes.length; at += IN_BYTES) {
      if (at > 0) await nextTask();
      try {
        // A copy: a stream may keep a chunk after its write settles (Node's does), and the caller
        // may reuse `bytes` then (the Writer reuses its block).
        await this.writer.write(bytes.slice(at, at + IN_BYTES) as BufferSource);
      } catch (error) {
        throw this.failure ?? error;
      }
    }
  }

  async close(): Promise<void> {
    if (this.state !== 'open') throw new Error('the compressed file is no longer open');
    this.state = 'closed';
    try {
      await this.writer.close();
    } catch (error) {
      throw this.failure ?? error;
    }
    await this.pumping;
  }

  abort(): void {
    if (this.state !== 'open') return;
    this.state = 'aborted';
    void this.writer.abort(new Error('the compressed file was aborted')).catch(() => undefined);
    this.pumping.catch(() => undefined);
  }
}
