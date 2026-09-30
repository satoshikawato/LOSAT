// The Data worker (plan §3.1): it owns the sources, the record tables, the run registry
// and the temporary storage of one working session, and it outlives Engine workers.
// Only this worker uses OPFS synchronous access handles. This module is the worker's
// composition root.
import { sha256Hex } from '../browser/platform';
import { DataService } from '../data/data-service';
import { openOpfsSession } from '../data/opfs-block-store';
import { startDataSession, type SessionLocks } from '../data/session';
import { FakeScanner } from '../fake/fake-fasta';
import { DATA_GATEWAY_METHODS } from './methods';
import { serveRpc, type RpcEndpoint } from './rpc';

const token = crypto.randomUUID();
const locks = 'locks' in navigator ? (navigator.locks as unknown as SessionLocks) : undefined;

const service = startDataSession({ token, locks, openOpfs: openOpfsSession }).then(
  (session) =>
    new DataService({
      store: session.store,
      ...(session.fallbackReason === undefined ? {} : { fallbackReason: session.fallbackReason }),
      cleanup: session.cleanup,
      // S09 replaces the FakeScanner with the ABI v2 scan of the adapter's serial reactor.
      scanner: new FakeScanner(),
      digest: sha256Hex,
      newToken: () => crypto.randomUUID(),
      estimate: estimateStorage,
    }),
);

serveRpc(self as unknown as RpcEndpoint, service, DATA_GATEWAY_METHODS);

async function estimateStorage(): Promise<{ usage: number; quota: number } | undefined> {
  if (typeof navigator.storage?.estimate !== 'function') return undefined;
  const { usage, quota } = await navigator.storage.estimate();
  return usage === undefined || quota === undefined ? undefined : { usage, quota };
}
