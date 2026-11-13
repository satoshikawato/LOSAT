// The worker of the memory keeper (memory-keeper.ts).
import { MemoryKeeper } from './memory-keeper.js';
import type { KeeperCommand } from './protocol.js';

const scope = self as unknown as {
  onmessage: ((event: MessageEvent<KeeperCommand>) => void) | null;
  postMessage(message: unknown): void;
};
const keeper = new MemoryKeeper();
scope.onmessage = (event) => scope.postMessage(keeper.receive(event.data));
