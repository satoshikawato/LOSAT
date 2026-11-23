// Streams a tool run's stdout, stderr and created files over the dedicated MessagePort returned
// by DataGateway.openToolRun so bytes go straight into session OPFS (tmp/) with transferred
// ArrayBuffers (zero accumulation on the page JS heap) and awaited acks for backpressure.
import type { ToolFileSink, ToolSinks } from 'losat/ports';
import type { ToolOutputCommand, ToolOutputReply, ToolOutputsSummary } from '../ports/data';

interface Pending {
  readonly resolve: (reply: Extract<ToolOutputReply, { ok: true }>) => void;
  readonly reject: (error: Error) => void;
}

export class ToolStorageSink {
  private nextId = 1;
  private nextFile = 0;
  private readonly pending = new Map<number, Pending>();
  private closed = false;
  private broken: Error | undefined;

  readonly sinks: ToolSinks;

  constructor(private readonly port: MessagePort) {
    port.onmessage = (event: MessageEvent<ToolOutputReply>) => {
      const reply = event.data;
      const entry = this.pending.get(reply.id);
      if (entry === undefined) return;
      this.pending.delete(reply.id);
      if (reply.ok) {
        entry.resolve(reply);
      } else {
        const error = new Error(reply.error.message);
        error.name = reply.error.name;
        entry.reject(error);
      }
    };
    port.onmessageerror = () => {
      this.fail(new Error('The output channel to temporary storage stopped'));
    };

    const sendBytes = async (
      chunk: Uint8Array,
      make: (bytes: Uint8Array) => ToolOutputCommand,
    ): Promise<void> => {
      if (chunk.length === 0) return;
      const { bytes, transfer } = transferable(chunk);
      await this.send(make(bytes), transfer);
    };

    this.sinks = {
      onStdout: (chunk) => sendBytes(chunk, (bytes) => ({ id: 0, op: 'stdout', bytes })),
      onStderr: (chunk) => sendBytes(chunk, (bytes) => ({ id: 0, op: 'stderr', bytes })),
      onFile: async (name) => {
        const key = `file:${this.nextFile++}`;
        await this.send({ id: 0, op: 'open', key, name }, []);
        const fileSink: ToolFileSink = {
          write: (chunk, offset) => sendBytes(chunk, (bytes) => ({ id: 0, op: 'write', key, offset, bytes })),
          truncate: async (length) => {
            await this.send({ id: 0, op: 'truncate', key, length }, []);
          },
          close: async () => {
            await this.send({ id: 0, op: 'close', key }, []);
          },
          abort: () => undefined,
        };
        return fileSink;
      },
    };
  }

  /** Seals all staged streams/files and closes the port. */
  async commit(): Promise<ToolOutputsSummary> {
    try {
      const reply = await this.send({ id: 0, op: 'commit' }, []);
      if (reply.summary === undefined) throw new Error('temporary storage returned no summary');
      return reply.summary;
    } finally {
      this.close();
    }
  }

  /** Drops staged blocks and closes the port. Safe to call after a failure or commit. */
  async discard(): Promise<void> {
    if (this.closed) return;
    try {
      if (this.broken === undefined) await this.send({ id: 0, op: 'discard' }, []);
    } catch {
      // Discard best-effort; caller also cleans up via deleteRun.
    } finally {
      this.close();
    }
  }

  /** Fails every pending ack immediately (used on cancel or Data worker end). */
  fail(reason: unknown): void {
    const error = reason instanceof Error ? reason : new Error(String(reason));
    if (this.broken === undefined) this.broken = error;
    for (const entry of this.pending.values()) entry.reject(error);
    this.pending.clear();
    this.close();
  }

  private close(): void {
    if (this.closed) return;
    this.closed = true;
    this.port.onmessage = null;
    this.port.onmessageerror = null;
    this.port.close();
  }

  private send(command: ToolOutputCommand, transfer: Transferable[]): Promise<Extract<ToolOutputReply, { ok: true }>> {
    if (this.broken !== undefined) return Promise.reject(this.broken);
    if (this.closed) return Promise.reject(new Error('tool output channel is closed'));
    const id = this.nextId++;
    const payload = { ...command, id } as ToolOutputCommand;
    return new Promise((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      try {
        this.port.postMessage(payload, transfer);
      } catch (error) {
        this.pending.delete(id);
        reject(error instanceof Error ? error : new Error(String(error)));
      }
    });
  }
}

/** Transfers the chunk's buffer when it owns the whole ArrayBuffer, otherwise copies to a standalone buffer. */
function transferable(chunk: Uint8Array): { readonly bytes: Uint8Array; readonly transfer: Transferable[] } {
  if (chunk.buffer instanceof ArrayBuffer && chunk.byteOffset === 0 && chunk.byteLength === chunk.buffer.byteLength) {
    return { bytes: chunk, transfer: [chunk.buffer] };
  }
  const owned = chunk.slice();
  return { bytes: owned, transfer: [owned.buffer] };
}
