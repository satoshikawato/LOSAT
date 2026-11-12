// The keeper of the threaded instance's shared memory (thread-host.ts makes one for each
// memory): a worker of the Engine worker that holds the memory and runs nothing else, so that
// it can grow the memory while the Engine worker and the thread workers run a search. Before
// the Gateway ends the Engine worker while its threaded instance is in use (a cancel, a thread
// that reports a fault; instance-use.ts), it asks the keeper on a BroadcastChannel to grow the
// memory to its maximum (gateway.ts, retire; growToMaximum).
//
// The main thread does not hold the memory itself. A thread that is sent a shared memory keeps
// it until that thread collects the message. In WebKit 26.5 (Playwright's build on Linux), a
// page held at most 15 shared memories of 4 GiB, and a new memory did not make the main thread
// collect: the memories that the main thread had been sent, one for each runtime, kept the
// threaded engine from starting after some runtimes. The keeper ends with the thread workers
// of its memory (ThreadHost.terminate), and with it its reference.
//
// The Gateway and the keeper each open the channel by its name. A MessagePort would have to
// be moved to them by the Engine worker, which then runs a search synchronously. In
// Playwright's Firefox 151, a port that a worker moved carried no message while that worker
// ran synchronously, also after a setTimeout of 0 ms before it; after one of 5 ms, it did.
import { growToMaximum } from '../reactor/shared-memory';
import type { KeeperCommand, KeeperGrow, KeeperGrown, KeeperReady } from './protocol';

export class MemoryKeeper {
  private channel: BroadcastChannel | undefined;

  /** Keeps the memory and listens on the channel; returns the answer that the worker posts back. */
  receive(command: KeeperCommand): KeeperReady {
    const shared = { memory: command.memory, maximumPages: command.maximumPages };
    this.channel?.close();
    const channel = new BroadcastChannel(command.channel);
    this.channel = channel;
    channel.onmessage = (event: MessageEvent<KeeperGrow>) => {
      if (event.data?.type !== 'grow') return;
      try {
        growToMaximum(shared);
      } catch {
        // The Gateway ends a runtime whose threaded instance is in use only on the answer
        // (gateway.ts, retire), so it comes also then; the memory stays as it is, as when the
        // browser cannot give its maximum.
      }
      const grown: KeeperGrown = { type: 'grown' };
      channel.postMessage(grown);
    };
    return { type: 'ready' };
  }

  /** Stops listening (the worker is ended instead; tests close it). */
  close(): void {
    this.channel?.close();
    this.channel = undefined;
  }
}
