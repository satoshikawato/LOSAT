// The links between the main thread (the Gateway) and the thread workers and the keeper of the shared memory, which the
// Engine worker starts: a thread worker reports on its link a thread that stopped, and the Gateway asks the keeper on its
// link to grow the memory before it ends a runtime in use (gateway.ts, retire). Both must work while the Engine worker
// runs a search, whose call into the engine is synchronous: the Engine worker takes no message then.
//
// A link is a MessageChannel that the Engine worker makes (MainLinks.link): it posts one end to the main thread (a 'link'
// event) and gives the other to the thread worker or the keeper. The link carries messages only once the worker's end
// said hello and the main thread's answer came back (openLink, answerLink). In Firefox 159 (Playwright 1.65), a port
// that a worker moved, to a nested worker that already ran, carried no message while the moving worker then ran
// synchronously for 3 s, and all of them once it returned; after a message had gone both ways first, it carried them
// during that time, in Chromium 158, Firefox 159 and WebKit 27.2 (losat-baselines/web-qa-20261119/fix-d8/probe/
// probe-started*.log). An earlier note (memory-keeper.ts, Playwright's Firefox 151) had seen the same also after a
// setTimeout of 0 ms before the synchronous run, and not after one of 5 ms. The ThreadHost therefore starts no thread
// before every thread worker and the keeper have had their answer (thread-host.ts, prepare).
//
// The links were BroadcastChannels, opened by name on each side. Firefox refuses `new BroadcastChannel()` with
// SecurityError ("The operation is insecure.") where the site's data is blocked (cookies blocked), on the page and in
// every worker, and so every search failed with those words (line WQA's sweep, D8, 2026-11-20). A MessageChannel is
// not the site's data: no browser refuses it.
import type { LinkHello, LinkMessage, LinkRole, ThreadFault } from './protocol.js';

/** How the ThreadHost reaches the main thread, also while the Engine worker runs a search synchronously. */
export interface MainLinks {
  /** A new link for a thread worker or the keeper: the main thread has the other end. */
  link(role: LinkRole): MessagePort;
  /** Tells the main thread of a fault of the ThreadHost itself (a thread that did not start in time). */
  fault(fault: ThreadFault): void;
}

const HELLO: LinkHello = { type: 'hello' };

/**
 * The end of a link that a thread worker or the keeper holds: says hello and resolves once the main thread's answer came
 * back, from when on the link carries messages also while the Engine worker runs synchronously. `receive` gets every
 * other message.
 */
export function openLink(port: MessagePort, receive: (message: LinkMessage) => void = () => undefined): Promise<void> {
  return new Promise<void>((resolve) => {
    port.onmessage = (event: MessageEvent<LinkMessage>) => {
      if (event.data?.type === 'hello') resolve();
      else receive(event.data);
    };
    port.postMessage(HELLO);
  });
}

/** The main thread's end of a link: answers its hello, and gives `receive` every other message. */
export function answerLink(port: MessagePort, receive: (message: LinkMessage) => void): void {
  port.onmessage = (event: MessageEvent<LinkMessage>) => {
    if (event.data?.type === 'hello') port.postMessage(HELLO);
    else receive(event.data);
  };
}
