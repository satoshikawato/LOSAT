// Runs the BlockStore contract in a dedicated worker, where OPFS synchronous access handles
// exist. The OPFS store works in a scratch directory outside tmp/, so the application's
// start-up cleanup never touches it.
import { MemoryBlockStore } from '../../../src/infra/data/memory-block-store';
import { OpfsBlockStore } from '../../../src/infra/data/opfs-block-store';
import { BLOCK_STORE_CASES, type BlockStoreEnv } from '../../contract/block-store.contract';
import { runCases, type CaseResult } from '../../contract/contract';
import type { BlockStorePageMessage, BlockStoreWorkerMessage } from './protocol';

const scope = self as unknown as {
  postMessage(message: BlockStoreWorkerMessage): void;
  onmessage: ((event: MessageEvent<BlockStorePageMessage>) => void) | null;
};

let nextId = 0;

/** Reports each case as it finishes, so a case that hangs can be found. */
function progress(result: CaseResult): void {
  scope.postMessage({ type: 'progress', result });
}

const waiting = new Map<number, () => void>();

/** Asks the page to change the origin's quota (Playwright does it through CDP). */
function ask(type: 'exhaust' | 'restore'): Promise<void> {
  const id = ++nextId;
  return new Promise((resolve) => {
    waiting.set(id, resolve);
    scope.postMessage({ type, id });
  });
}

async function runOpfs(): Promise<CaseResult[]> {
  const root = await navigator.storage.getDirectory();
  const scratch = await root.getDirectoryHandle(`contract-${crypto.randomUUID()}`, { create: true });
  let count = 0;
  try {
    return await runCases(
      BLOCK_STORE_CASES,
      async (): Promise<BlockStoreEnv> => ({
        store: new OpfsBlockStore(await scratch.getDirectoryHandle(`case-${++count}`, { create: true })),
        exhaust: () => ask('exhaust'),
        restore: () => ask('restore'),
      }),
      { onResult: progress },
    );
  } finally {
    // A failed case can leave a block open, which keeps its directory in use.
    await root.removeEntry(scratch.name, { recursive: true }).catch(() => undefined);
  }
}

async function runMemory(): Promise<CaseResult[]> {
  return runCases(
    BLOCK_STORE_CASES,
    (): BlockStoreEnv => {
      const store = new MemoryBlockStore();
      return {
        store,
        exhaust: async () => store.setCapacity(store.usage()),
        restore: async () => store.setCapacity(Number.POSITIVE_INFINITY),
      };
    },
    { onResult: progress },
  );
}

scope.onmessage = (event) => {
  const message = event.data;
  if (message.type === 'done') {
    waiting.get(message.id)?.();
    waiting.delete(message.id);
    return;
  }
  void (message.backend === 'opfs' ? runOpfs() : runMemory()).then(
    (results) => scope.postMessage({ type: 'results', results }),
    (error: unknown) =>
      scope.postMessage({ type: 'results', results: [{ name: 'the harness', ok: false, error: String(error) }] }),
  );
};
