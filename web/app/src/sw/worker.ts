// The site's service worker, `/sw.js` (plan §5.9, design §14), built by
// build/service-worker.ts and registered by the boot page once a release has started. It is
// site-level: one worker serves the pages of every release, so it speaks the versioned
// messages of protocol.ts.
// - It stores a release only when a page asks (`prepare`), after checking every file
//   (prepare.ts); it answers a release's files from a complete cache and the start page from
//   the network first, else from the stored shell (serve.ts); it removes a release's cache only
//   when no tab holds that release's lock (cleanup.ts). Nothing else is answered or stored.
// - Lifecycle: nothing is stored on install, and there is no unconditional skipWaiting. A newer
//   worker waits until no page uses this one, unless a page whose user chose to switch releases
//   tells it to take over (`skip-waiting`). On activate it claims the pages that have no worker
//   (a first visit): it answers only release files and the shell, so nothing changes for them.
import { parseReleaseManifest } from '../domain/release-manifest';
import { parseSiteIndex, RELEASE_ID_PATTERN } from '../domain/site-index';
import { releaseCachesToRemove } from './cleanup';
import { prepareRelease, type PrepareDeps, type PrepareOutcome } from './prepare';
import { PROTOCOL_VERSION, RELEASE_CACHE_PREFIX, type PrepareProgress, type WorkerMessage } from './protocol';
import { routeOf, serveRelease, serveShell, type ServeDeps } from './serve';

declare const self: ServiceWorkerGlobalScope;

const origin = self.location.origin;
const message = (error: unknown): string => (error instanceof Error ? error.message : String(error));

async function sha256(bytes: ArrayBuffer): Promise<string> {
  const digest = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
  return Array.from(digest, (byte) => byte.toString(16).padStart(2, '0')).join('');
}

const prepareDeps: PrepareDeps = {
  origin,
  fetch: (url, init) => fetch(url, init),
  caches,
  sha256,
  // An index without a shell lists no start page: the preparation then fails and says so.
  readShell: (index) => parseSiteIndex(index).shell ?? [],
  readManifest: (manifest, release) => parseReleaseManifest(manifest, release).assets,
};
const serveDeps: ServeDeps = { origin, caches, fetch: (request) => fetch(request) };

function post(port: MessagePort, reply: WorkerMessage): void {
  port.postMessage(reply);
}

/** A preparation, queued or running, and the pages that wait for it. */
interface Job {
  readonly ports: Set<MessagePort>;
  last: PrepareProgress | undefined;
  done: Promise<PrepareOutcome>;
}

/** By release; a release being prepared is never removed by a cleanup. */
const jobs = new Map<string, Job>();
/** One preparation at a time: two would both write the shell. */
let queue: Promise<unknown> = Promise.resolve();

function startJob(release: string): Job {
  const job: Job = { ports: new Set(), last: undefined, done: Promise.resolve({ ok: false, reason: 'not started' }) };
  const run = async (): Promise<PrepareOutcome> => {
    try {
      return await prepareRelease(prepareDeps, release, (progress) => {
        job.last = progress;
        for (const port of job.ports) post(port, { type: 'progress', progress });
      });
    } catch (error) {
      return { ok: false, reason: message(error) };
    } finally {
      jobs.delete(release);
    }
  };
  job.done = queue.then(run, run);
  queue = job.done;
  jobs.set(release, job);
  return job;
}

async function cleanup(keep: readonly string[]): Promise<string[]> {
  const locks = 'locks' in navigator ? await navigator.locks.query().catch(() => undefined) : undefined;
  const remove = releaseCachesToRemove(await caches.keys(), locks, keep, jobs.keys());
  for (const name of remove) await caches.delete(name);
  return remove.map((name) => name.slice(RELEASE_CACHE_PREFIX.length));
}

self.addEventListener('activate', (event) => event.waitUntil(self.clients.claim()));

self.addEventListener('fetch', (event) => {
  const route = routeOf(event.request, origin);
  if (route.kind === 'release') event.respondWith(serveRelease(serveDeps, event.request, route.release));
  else if (route.kind === 'shell') event.respondWith(serveShell(serveDeps, event.request, route.path));
});

self.addEventListener('message', (event) => {
  const port = event.ports[0];
  const data = event.data as { v?: unknown; type?: unknown; release?: unknown; keep?: unknown } | null;
  if (data?.v !== PROTOCOL_VERSION) {
    if (port !== undefined) post(port, { type: 'refused', reason: 'the offline worker of this site speaks another version' });
    return;
  }
  if (data.type === 'skip-waiting') {
    event.waitUntil(self.skipWaiting());
    return;
  }
  if (port === undefined) return;
  if (data.type === 'prepare') {
    const release = data.release;
    if (typeof release !== 'string' || !RELEASE_ID_PATTERN.test(release)) {
      post(port, { type: 'refused', reason: 'not a release ID' });
      return;
    }
    const job = jobs.get(release) ?? startJob(release);
    job.ports.add(port);
    if (job.last !== undefined) post(port, { type: 'progress', progress: job.last });
    event.waitUntil(
      job.done.then((outcome) => {
        job.ports.delete(port);
        post(port, outcome.ok ? { type: 'prepared', ok: true } : { type: 'prepared', ok: false, reason: outcome.reason });
      }),
    );
    return;
  }
  if (data.type === 'cleanup') {
    const keep = Array.isArray(data.keep) ? data.keep.filter((id): id is string => typeof id === 'string') : [];
    event.waitUntil(
      cleanup(keep).then(
        (removed) => post(port, { type: 'cleaned', removed }),
        (error: unknown) => post(port, { type: 'refused', reason: message(error) }),
      ),
    );
    return;
  }
  post(port, { type: 'refused', reason: 'unknown message' });
});
