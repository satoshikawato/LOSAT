// Engine side of the run output channel (ports/run-output.ts). The Engine worker (S09) and
// the FakeEngine write through it.
import type { OutputStream, RunOutputMessage } from '../../ports/run-output';

export class RunOutputWriter {
  private chunks = 0;
  private bytes = 0;
  private ended = false;

  constructor(private readonly port: MessagePort) {}

  /**
   * Sends a copy of `chunk` and transfers the copy. The caller may reuse its buffer when
   * `write` returns, so a view into Wasm memory that is valid only during an ABI `emit`
   * call can be passed as is. Empty chunks are not sent.
   */
  write(stream: OutputStream, chunk: Uint8Array): void {
    if (this.ended) throw new Error('the output of this run has already ended');
    if (chunk.length === 0) return;
    const bytes = chunk.slice();
    const message: RunOutputMessage = { type: 'chunk', stream, bytes };
    // Count before the transfer detaches the copy.
    this.chunks++;
    this.bytes += bytes.length;
    this.port.postMessage(message, [bytes.buffer]);
  }

  /** Sends `end`. Call it once, after the last chunk and before the run reports success. */
  end(): void {
    if (this.ended) throw new Error('the output of this run has already ended');
    this.ended = true;
    const message: RunOutputMessage = { type: 'end', chunks: this.chunks, bytes: this.bytes };
    this.port.postMessage(message);
  }
}
