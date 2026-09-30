// Contract harness page for tests/e2e/contracts.spec.ts. Playwright builds it in memory and
// serves it under /__harness/; it is never part of the application build. It runs the
// contract suites in real browser workers:
// - the BlockStore contract against OPFS and memory in a dedicated worker;
// - the run output contract with the writer in a separate worker (the Engine worker's
//   place) and the real Data worker of the application as the receiver.
import { startDataWorker } from '../../../src/infra/data-worker/gateway';
import type { StorageBackend } from '../../../src/ports/data';
import type { OutputStream } from '../../../src/ports/run-output';
import { runCases, type CaseResult } from '../../contract/contract';
import { RUN_OUTPUT_CASES, type RemoteWriter } from '../../contract/run-output.contract';
import type { BlockStorePageMessage, BlockStoreWorkerMessage, EngineDoubleCommand, EngineDoubleReply } from './protocol';

export interface Harness {
  blockStore(backend: 'opfs' | 'memory'): Promise<CaseResult[]>;
  runOutput(): Promise<{ readonly backend: StorageBackend; readonly results: CaseResult[] }>;
}

declare global {
  interface Window {
    losatHarness?: Harness;
    /** Exposed by Playwright: sets the origin's storage quota in bytes, or resets it (null). */
    setStorageQuota?: (bytes: number | null) => Promise<void>;
  }
}

/**
 * Makes the origin's storage refuse more data. The quota is set to zero rather than to
 * the current usage, because `navigator.storage.estimate()` can lag behind the usage that
 * Chromium checks, which would leave room for more data.
 */
async function exhaust(): Promise<void> {
  await window.setStorageQuota!(0);
}

function logResult(label: string, result: CaseResult): void {
  console.log(`${result.ok ? 'pass' : 'FAIL'} [${label}] ${result.name}${result.error ? `: ${result.error}` : ''}`);
}

async function restore(): Promise<void> {
  await window.setStorageQuota!(null);
}

function blockStore(backend: 'opfs' | 'memory'): Promise<CaseResult[]> {
  const worker = new Worker(new URL('./block-store-worker.ts', import.meta.url), { type: 'module' });
  const post = (message: BlockStorePageMessage) => worker.postMessage(message);
  return new Promise((resolve, reject) => {
    worker.onerror = (event) => reject(new Error(`the block store worker failed: ${event.message}`));
    worker.onmessage = (event: MessageEvent<BlockStoreWorkerMessage>) => {
      const message = event.data;
      if (message.type === 'results') {
        worker.terminate();
        resolve(message.results);
        return;
      }
      if (message.type === 'progress') {
        logResult(backend, message.result);
        return;
      }
      void (message.type === 'exhaust' ? exhaust() : restore()).then(() => post({ type: 'done', id: message.id }));
    };
    post({ type: 'start', backend });
  });
}

/** The page side of the engine double worker. */
class EngineDouble {
  private readonly worker = new Worker(new URL('./engine-double-worker.ts', import.meta.url), { type: 'module' });
  private readonly pending = new Map<number, { resolve(value: number | undefined): void; reject(error: Error): void }>();
  private nextId = 0;

  constructor() {
    this.worker.onmessage = (event: MessageEvent<EngineDoubleReply>) => {
      const reply = event.data;
      const call = this.pending.get(reply.id);
      this.pending.delete(reply.id);
      if (reply.ok) {
        call?.resolve(reply.value);
      } else {
        const error = new Error(reply.message);
        error.name = reply.name;
        call?.reject(error);
      }
    };
  }

  async writer(port: MessagePort): Promise<RemoteWriter> {
    const writer = (await this.call({ type: 'open', port }, [port]))!;
    return {
      write: async (stream: OutputStream, bytes: Uint8Array) => {
        await this.call({ type: 'write', writer, stream, bytes });
      },
      end: async () => {
        await this.call({ type: 'end', writer });
      },
      post: async (message: unknown) => {
        await this.call({ type: 'post', writer, message });
      },
    };
  }

  terminate(): void {
    this.worker.terminate();
  }

  private call(command: DistributiveOmit<EngineDoubleCommand, 'id'>, transfer: Transferable[] = []) {
    const id = ++this.nextId;
    return new Promise<number | undefined>((resolve, reject) => {
      this.pending.set(id, { resolve, reject });
      this.worker.postMessage({ ...command, id }, transfer);
    });
  }
}

type DistributiveOmit<T, K extends PropertyKey> = T extends unknown ? Omit<T, K> : never;

async function runOutput(): Promise<{ backend: StorageBackend; results: CaseResult[] }> {
  // One Data worker serves every case; the cases use distinct run ids.
  const data = startDataWorker();
  const engine = new EngineDouble();
  try {
    const backend = (await data.storageInfo()).backend;
    const results = await runCases(
      RUN_OUTPUT_CASES,
      () => ({ data, writer: (port: MessagePort) => engine.writer(port), exhaust, restore }),
      { onResult: (result) => logResult('run output', result) },
    );
    return { backend, results };
  } finally {
    engine.terminate();
  }
}

window.losatHarness = { blockStore, runOutput };
