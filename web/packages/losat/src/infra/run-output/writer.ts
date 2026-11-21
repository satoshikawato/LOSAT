// Engine side of the run output channel (ports/run-output.ts). The Engine worker (S09) and
// the FakeEngine write through it.
import type { RunOutput } from '../../domain/output-format.js';
import type { OutputPlace } from '../../domain/run-files.js';
import type { OutputStream, RunOutputMessage } from '../../ports/run-output.js';

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

  /**
   * Tells the receiver that the output `stream` failed, where its error is (`error`), and to drop
   * what it received of it, which then stays empty: the engine's signal of stream 9
   * (docs/web/abi_v2.md §5), for each output that failed alone or with its search in a run of two
   * searches, also where none of its chunks was sent. The chunks sent still count in the totals of
   * `end`.
   */
  discard(stream: RunOutput, error: OutputPlace): void {
    if (this.ended) throw new Error('the output of this run has already ended');
    const message: RunOutputMessage = { type: 'discard', stream, error };
    this.port.postMessage(message);
  }

  /** Sends `end`. Call it once, after the last chunk and before the run reports success. */
  end(): void {
    this.finish('end');
  }

  /**
   * Sends `fail` in place of `end`: the engine failed the run after it began to write. Call it once,
   * before the run's rejection, which then says so (`withRunOutputFailed`, ports/engine.ts).
   */
  fail(): void {
    this.finish('fail');
  }

  private finish(type: 'end' | 'fail'): void {
    if (this.ended) throw new Error('the output of this run has already ended');
    this.ended = true;
    const message: RunOutputMessage = { type, chunks: this.chunks, bytes: this.bytes };
    this.port.postMessage(message);
  }
}
