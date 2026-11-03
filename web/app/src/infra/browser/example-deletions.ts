// Deletions of saved Example copies between the tabs of the site (ports/example.ts), over a
// BroadcastChannel. The message carries only the set's id from the table. A browser without
// BroadcastChannel tells no other tab.
import type { ExampleDeletions } from '../../ports/example';

const CHANNEL = 'losat-web:examples';

export function browserExampleDeletions(): ExampleDeletions {
  if (typeof BroadcastChannel !== 'function') return { announce: () => undefined, subscribe: () => () => undefined };
  const channel = new BroadcastChannel(CHANNEL);
  return {
    announce: (id) => channel.postMessage({ deleted: id }),
    subscribe(listener) {
      const onMessage = (event: MessageEvent) => {
        const id = (event.data as { deleted?: unknown } | null)?.deleted;
        if (typeof id === 'string') listener(id);
      };
      channel.addEventListener('message', onMessage);
      return () => channel.removeEventListener('message', onMessage);
    },
  };
}
