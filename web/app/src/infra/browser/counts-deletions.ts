// Deletions of saved WindowMasker counts between the tabs of the site (ports/counts.ts), over a
// BroadcastChannel. The message carries only the taxonomy ID. A browser without BroadcastChannel,
// or one that refuses it (site-channel.ts), tells no other tab.
import type { CountsDeletions } from '../../ports/counts';
import { siteChannel } from './site-channel';

const CHANNEL = 'losat-web:windowmasker';

export function browserCountsDeletions(): CountsDeletions {
  const channel = siteChannel(CHANNEL);
  if (channel === undefined) return { announce: () => undefined, subscribe: () => () => undefined };
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
