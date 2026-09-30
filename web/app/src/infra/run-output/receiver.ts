// Data worker side of the run output channel (ports/run-output.ts).
import { OUTPUT_STREAMS, type OutputStream, type RunOutputMessage } from '../../ports/run-output';

export type RunOutputResult =
  /** `end` arrived and its totals agree with what arrived. */
  | { readonly state: 'ended' }
  /** `end` arrived but the totals disagree, or a message was not part of the protocol. */
  | { readonly state: 'broken'; readonly detail: string }
  /** The receiver was closed before `end`. */
  | { readonly state: 'closed' };

export class RunOutputReceiver {
  /** Settles once: on `end`, on a broken message, or on `close`. */
  readonly finished: Promise<RunOutputResult>;
  private settle!: (result: RunOutputResult) => void;
  private done = false;
  private chunks = 0;
  private bytes = 0;

  /** `onChunk` receives each chunk in order until the output finishes; it must not throw. */
  constructor(
    private readonly port: MessagePort,
    onChunk: (stream: OutputStream, bytes: Uint8Array) => void,
  ) {
    this.finished = new Promise((resolve) => {
      this.settle = resolve;
    });
    port.onmessage = (event: MessageEvent<unknown>) => this.receive(event.data, onChunk);
  }

  /** Stops receiving; later messages are dropped. */
  close(): void {
    this.finish({ state: 'closed' });
  }

  private receive(message: unknown, onChunk: (stream: OutputStream, bytes: Uint8Array) => void): void {
    if (this.done) return;
    const data = message as Partial<RunOutputMessage> | null;
    if (data?.type === 'chunk' && OUTPUT_STREAMS.includes(data.stream as OutputStream) && data.bytes instanceof Uint8Array) {
      this.chunks++;
      this.bytes += data.bytes.length;
      onChunk(data.stream as OutputStream, data.bytes);
      return;
    }
    if (data?.type === 'end' && typeof data.chunks === 'number' && typeof data.bytes === 'number') {
      const complete = data.chunks === this.chunks && data.bytes === this.bytes;
      this.finish(
        complete
          ? { state: 'ended' }
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
