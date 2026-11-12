// Deletions of saved WindowMasker counts between the tabs of the site (ports/counts.ts), over a
// BroadcastChannel. The message carries only the taxonomy ID. A browser without BroadcastChannel
// tells no other tab.
import type { CountsDeletions } from '../../ports/counts';

const CHANNEL = 'losat-web:windowmasker';

export function browserCountsDeletions(): CountsDeletions {
  if (typeof BroadcastChannel !== 'function') return { announce: () => undefined, subscribe: () => () => undefined };
  const channel = new BroadcastChannel(CHANNEL);
  return {
    announce: (taxid) => channel.postMessage({ deleted: taxid }),
    subscribe(listener) {
      const onMessage = (event: MessageEvent) => {
        const taxid = (event.data as { deleted?: unknown } | null)?.deleted;
        if (typeof taxid === 'number' && Number.isSafeInteger(taxid)) listener(taxid);
      };
      channel.addEventListener('message', onMessage);
      return () => channel.removeEventListener('message', onMessage);
    },
  };
}
