// The keeper of the threaded instance's shared memory (thread-host.ts makes one for each
// memory): a worker of the Engine worker that holds the memory and runs nothing else, so that
// it can grow the memory while the Engine worker and the thread workers run a search. Before
// the Gateway ends the Engine worker while its threaded instance is in use (a cancel, a thread
// that reports a fault; instance-use.ts), it asks the keeper on the keeper's link (link.ts) to grow the
// memory to its maximum (gateway.ts, retire; growToMaximum) and to end the waits of its threads,
// which WebKit would otherwise keep alive with the memory (endWaits).
//
// The main thread does not hold the memory itself. A thread that is sent a shared memory keeps
// it until that thread collects the message. In WebKit 26.5 (build 2311 of Playwright 1.61.1, on Linux), a
// page held at most 15 shared memories of 4 GiB, and a new memory did not make the main thread
// collect: the memories that the main thread had been sent, one for each runtime, kept the
// threaded engine from starting after some runtimes. The keeper ends with the thread workers
// of its memory (ThreadHost.terminate), and with it its reference.
import { endWaits, growToMaximum } from '../reactor/shared-memory.js';
import { openLink } from './link.js';
import type { KeeperCommand, KeeperGrown, KeeperReady } from './protocol.js';

export class MemoryKeeper {
  private port: MessagePort | undefined;

  /**
   * Keeps the memory and listens on its link; resolves with the answer that the worker posts back once the link carries
   * messages (link.ts).
   */
  async receive(command: KeeperCommand): Promise<KeeperReady> {
    const shared = { memory: command.memory, maximumPages: command.maximumPages };
    this.port?.close();
    const { port } = command;
    this.port = port;
    await openLink(port, (message) => {
      if (message?.type !== 'grow') return;
      try {
        growToMaximum(shared);
      } catch {
        // The Gateway ends a runtime whose threaded instance is in use only on the answer
        // (gateway.ts, retire), so it comes also then; the memory stays as it is, as when the
        // browser cannot give its maximum.
      }
      try {
        // The threads that wait in the threaded module wait in slices from now on, so that the
        // browser can end their workers (WebKit, endWaits).
        endWaits(shared.memory);
      } catch {
        // As for the grow: the answer comes all the same.
      }
      const grown: KeeperGrown = { type: 'grown' };
      port.postMessage(grown);
    });
    return { type: 'ready' };
  }

  /** Stops listening (the worker is ended instead; tests close it). */
  close(): void {
    this.port?.close();
    this.port = undefined;
  }
}
