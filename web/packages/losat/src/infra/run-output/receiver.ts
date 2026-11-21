// Data worker side of the run output channel (ports/run-output.ts).
import type { RunOutput } from '../../domain/output-format.js';
import { isOutputPlace, type OutputPlace } from '../../domain/run-files.js';
import { isRunOutput, OUTPUT_STREAMS, type OutputStream, type RunOutputMessage } from '../../ports/run-output.js';

export type RunOutputResult =
  /** `end` arrived and its totals agree with what arrived. */
  | { readonly state: 'ended' }
  /** `fail` arrived and its totals agree with what arrived: the engine failed the run. */
  | { readonly state: 'failed' }
  /** `end` or `fail` arrived but the totals disagree, or a message was not part of the protocol. */
  | { readonly state: 'broken'; readonly detail: string }
  /** The receiver was closed before `end`. */
  | { readonly state: 'closed' };

export class RunOutputReceiver {
  /** Settles once: on `end`, on `fail`, on a broken message, or on `close`. */
  readonly finished: Promise<RunOutputResult>;
  private settle!: (result: RunOutputResult) => void;
  private done = false;
  private chunks = 0;
  private bytes = 0;

  /**
   * `onChunk` receives each chunk in order until the output finishes, and `onDiscard` each
   * `discard` among them (drop what arrived of that output, whose error is at `error`); neither may
   * throw.
   */
  constructor(
    private readonly port: MessagePort,
    onChunk: (stream: OutputStream, bytes: Uint8Array) => void,
    onDiscard: (stream: RunOutput, error: OutputPlace) => void,
  ) {
    this.finished = new Promise((resolve) => {
      this.settle = resolve;
    });
    port.onmessage = (event: MessageEvent<unknown>) => this.receive(event.data, onChunk, onDiscard);
  }

  /** Stops receiving; later messages are dropped. */
  close(): void {
    this.finish({ state: 'closed' });
  }

  private receive(
    message: unknown,
    onChunk: (stream: OutputStream, bytes: Uint8Array) => void,
    onDiscard: (stream: RunOutput, error: OutputPlace) => void,
  ): void {
    if (this.done) return;
    const data = message as Partial<RunOutputMessage> | null;
    if (data?.type === 'chunk' && OUTPUT_STREAMS.includes(data.stream as OutputStream) && data.bytes instanceof Uint8Array) {
      this.chunks++;
      this.bytes += data.bytes.length;
      onChunk(data.stream as OutputStream, data.bytes);
      return;
    }
    if (data?.type === 'discard' && isRunOutput(data.stream) && isOutputPlace(data.error)) {
      onDiscard(data.stream, data.error);
      return;
    }
    if ((data?.type === 'end' || data?.type === 'fail') && typeof data.chunks === 'number' && typeof data.bytes === 'number') {
      const complete = data.chunks === this.chunks && data.bytes === this.bytes;
      this.finish(
        complete
          ? { state: data.type === 'end' ? 'ended' : 'failed' }
          : {
              state: 'broken',
              detail:
                `the engine sent ${data.chunks} chunks (${data.bytes} bytes), ` +
                `but ${this.chunks} chunks (${this.bytes} bytes) arrived`,
            },
      );
      return;
    }
    this.finish({ state: 'broken', detail: 'a message that is not part of the run output protocol arrived' });
  }

  private finish(result: RunOutputResult): void {
    if (this.done) return;
    this.done = true;
    this.port.onmessage = null;
    this.port.close();
    this.settle(result);
  }
}
