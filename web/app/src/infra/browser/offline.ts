// The browser side of offline use and the explicit update (ports/offline.ts): the client of
// the site's service worker (src/sw/), the release lock, the pinned release in localStorage
// (the key that the boot page reads), the BroadcastChannel of the site's tabs, and the
// page's network state. Storage that throws or is missing never stops the app.
import { RELEASE_ID_PATTERN } from '../../domain/site-index';
import type { OfflinePort, PinWrite, PrepareResult, ReleaseChannel, ReleasePage, ReleasePin } from '../../ports/offline';
import {
  MANIFEST_NAME,
  PROTOCOL_VERSION,
  RELEASE_LOCK_PREFIX,
  SW_SCOPE,
  SW_URL,
  releaseCacheName,
  releasePath,
  type PageMessage,
  type PrepareProgress,
  type WorkerMessage,
} from '../../sw/protocol';

/** The pinned release (DECISIONS.md 3); the boot page reads the same key. */
export const PIN_KEY = 'losat-web:release';
/** The channel on which a tab tells the others that it switched the pinned release. */
export const CHANNEL_NAME = 'losat-web:release';

/** How long the worker may take to become active after it is registered. */
const CONNECT_TIMEOUT_MS = 30_000;
/**
 * A preparation whose worker says nothing for this long has stopped (a browser may end a busy
 * worker). Longer than the worker's own limit for a download that stalls (src/sw/prepare.ts),
 * so that its reason arrives first.
 */
const SILENCE_MS = 150_000;
/** How long a cleanup may take before the page stops waiting for it. */
const CLEANUP_TIMEOUT_MS = 30_000;
/** How long a waiting worker may take to take over after it is told to. */
const ACTIVATE_TIMEOUT_MS = 3_000;
/** After the page's load event: a pause, then an idle period (or this long without one). */
const SETTLE_MS = 2_000;
const IDLE_TIMEOUT_MS = 5_000;

const message = (error: unknown): string => (error instanceof Error ? error.message : String(error));

/**
 * Whether a rejected `register()` means that the browser does not let this site use service
 * workers (cookies or site data blocked for the site). Then the browser refuses every call of
 * the container, `getRegistration()` too: Chromium with a NotSupportedError ("The user denied
 * permission to use Service Worker."), Firefox with a SecurityError ("The operation is
 * insecure."). A SecurityError also comes from a `sw.js` that is not JavaScript (a proxy's or a
 * captive portal's HTML page, in all three browsers; Firefox with the same message) and, in
 * Chromium, from a certificate error; then `getRegistration()` still answers, with or without a
 * registration, and the rejection stays a failure to try again (code review 3 L2; probed in fix
 * rounds 2 and 4 of S16). A page that a worker controls is never refused.
 */
async function refusesServiceWorkers(container: ServiceWorkerContainer, error: unknown): Promise<boolean> {
  if (!(error instanceof DOMException) || (error.name !== 'SecurityError' && error.name !== 'NotSupportedError')) return false;
  if ((container.controller ?? null) !== null) return false;
  try {
    await container.getRegistration(SW_SCOPE);
    return false;
  } catch {
    return true;
  }
}

function isProgress(value: unknown): value is PrepareProgress {
  if (typeof value !== 'object' || value === null) return false;
  const p = value as Record<string, unknown>;
  return ['files', 'filesTotal', 'bytes', 'bytesTotal'].every((key) => typeof p[key] === 'number');
}

/** A reply of the worker, or undefined for anything else (a worker of another version may say more). */
function workerMessage(value: unknown): WorkerMessage | undefined {
  if (typeof value !== 'object' || value === null) return undefined;
  const m = value as Record<string, unknown>;
  switch (m['type']) {
    case 'progress':
      return isProgress(m['progress']) ? { type: 'progress', progress: m['progress'] } : undefined;
    case 'prepared':
      return m['ok'] === true ? { type: 'prepared', ok: true } : { type: 'prepared', ok: false, reason: String(m['reason'] ?? 'unknown') };
    case 'cleaned':
      return { type: 'cleaned', removed: Array.isArray(m['removed']) ? m['removed'].map(String) : [] };
    case 'refused':
      return { type: 'refused', reason: String(m['reason'] ?? 'unknown') };
    default:
      return undefined;
  }
}

/** The service-worker client: one per page. */
export function browserOffline(): OfflinePort {
  const container = typeof navigator !== 'undefined' && 'serviceWorker' in navigator ? navigator.serviceWorker : undefined;
  let registration: ServiceWorkerRegistration | undefined;
  let blocked = false;

  const active = async (): Promise<ServiceWorker | undefined> => {
    if (container === undefined) return undefined;
    registration ??= await container.getRegistration(SW_SCOPE);
    return registration?.active ?? undefined;
  };

  /** Sends a message with a port, and resolves with what `onReply` returns once it returns something. */
  const ask = <T>(
    worker: ServiceWorker,
    request: PageMessage,
    silenceMs: number,
    onReply: (reply: WorkerMessage) => T | undefined,
    onSilence: () => T,
  ): Promise<T> =>
    new Promise((resolve) => {
      const channel = new MessageChannel();
      let timer = setTimeout(() => finish(onSilence()), silenceMs);
      const finish = (value: T) => {
        clearTimeout(timer);
        channel.port1.close();
        resolve(value);
      };
      channel.port1.onmessage = (event: MessageEvent<unknown>) => {
        const reply = workerMessage(event.data);
        if (reply === undefined) return;
        clearTimeout(timer);
        timer = setTimeout(() => finish(onSilence()), silenceMs);
        const value = onReply(reply);
        if (value !== undefined) finish(value);
      };
      worker.postMessage(request, [channel.port2]);
    });

  return {
    // False without service workers, and once a registration showed that they are blocked.
    get supported(): boolean {
      return container !== undefined && !blocked;
    },

    async connect(): Promise<void> {
      if (container === undefined) throw new Error('this browser has no service workers');
      let registered: ServiceWorkerRegistration | undefined;
      try {
        registered = (await container.register(SW_URL, { scope: SW_SCOPE, updateViaCache: 'none' })) as ServiceWorkerRegistration | undefined;
      } catch (error) {
        // A browser that does not let this site use service workers (cookies or site data blocked)
        // refuses every registration: that is no offline use, as without service workers, never a
        // failure to try again on each load (code review 2 L1). Any other rejection stays a failure,
        // so that an update is still downloaded and checked before the switch (code review 3 L2).
        if (await refusesServiceWorkers(container, error)) {
          blocked = true;
          throw new Error(`this browser does not let this site use service workers: ${message(error)}`, { cause: error });
        }
        throw error;
      }
      // Something that blocks service workers by replacing register() (Playwright's
      // `serviceWorkers: 'block'`, a privacy tool) may resolve with nothing instead of refusing:
      // say so, not what reading an undefined registration throws.
      if (registered === undefined) {
        blocked = true;
        throw new Error('service workers are blocked in this browser');
      }
      registration = registered;
      if (registration.active !== null) return;
      let timer: ReturnType<typeof setTimeout> | undefined;
      try {
        await Promise.race([
          container.ready,
          new Promise((_, reject) => {
            timer = setTimeout(() => reject(new Error('the offline worker did not start')), CONNECT_TIMEOUT_MS);
          }),
        ]);
      } finally {
        clearTimeout(timer);
      }
    },

    async prepare(release, onProgress): Promise<PrepareResult> {
      const worker = await active();
      if (worker === undefined) return { ok: false, reason: 'the offline worker is not running' };
      return ask<PrepareResult>(
        worker,
        { v: PROTOCOL_VERSION, type: 'prepare', release },
        SILENCE_MS,
        (reply) => {
          if (reply.type === 'progress') {
            onProgress(reply.progress);
            return undefined;
          }
          if (reply.type === 'prepared') return reply.ok ? { ok: true } : { ok: false, reason: reply.reason };
          if (reply.type === 'refused') return { ok: false, reason: reply.reason };
          return undefined;
        },
        () => ({ ok: false, reason: 'the offline worker stopped answering' }),
      );
    },

    async cleanup(keep): Promise<void> {
      const worker = await active();
      if (worker === undefined) return;
      await ask<true>(
        worker,
        { v: PROTOCOL_VERSION, type: 'cleanup', keep: [...keep] },
        CLEANUP_TIMEOUT_MS,
        (reply) => (reply.type === 'cleaned' || reply.type === 'refused' ? true : undefined),
        () => true,
      );
    },

    async activateWaiting(): Promise<boolean> {
      if (container === undefined) return false;
      registration ??= await container.getRegistration(SW_SCOPE);
      const waiting = registration?.waiting ?? null;
      if (waiting === null) return false;
      await new Promise<void>((resolve) => {
        const timer = setTimeout(resolve, ACTIVATE_TIMEOUT_MS);
        waiting.addEventListener('statechange', () => {
          if (waiting.state !== 'activated' && waiting.state !== 'redundant') return;
          clearTimeout(timer);
          resolve();
        });
        waiting.postMessage({ v: PROTOCOL_VERSION, type: 'skip-waiting' } satisfies PageMessage);
      });
      return true;
    },

    holdRelease(release): () => void {
      if (typeof navigator === 'undefined' || !('locks' in navigator)) return () => undefined;
      let letGo = (): void => undefined;
      const held = new Promise<void>((resolve) => (letGo = resolve));
      // Shared: every tab of the release holds it at once. The page's end releases it too.
      navigator.locks.request(RELEASE_LOCK_PREFIX + release, { mode: 'shared' }, () => held).catch(() => undefined);
      return () => letGo();
    },

    async stored(release): Promise<boolean> {
      if (typeof caches === 'undefined') return false;
      try {
        // caches.match with a cacheName never creates the cache, unlike caches.open.
        const manifest = await caches.match(location.origin + releasePath(release, MANIFEST_NAME), { cacheName: releaseCacheName(release) });
        return manifest !== undefined;
      } catch {
        return false;
      }
    },
  };
}

/** The pin in a Storage (localStorage in the page); every access may throw. */
export function releasePin(storage: () => Storage | undefined): ReleasePin {
  return {
    read(): string | undefined {
      try {
        const value = storage()?.getItem(PIN_KEY) ?? null;
        return value !== null && RELEASE_ID_PATTERN.test(value) ? value : undefined;
      } catch {
        return undefined;
      }
    },
    write(release): PinWrite {
      try {
        const store = storage();
        if (store === undefined) return { ok: false, reason: 'this browser gives the page no storage' };
        store.setItem(PIN_KEY, release);
        return store.getItem(PIN_KEY) === release ? { ok: true } : { ok: false, reason: 'the browser did not keep it' };
      } catch (error) {
        return { ok: false, reason: message(error) };
      }
    },
  };
}

export const browserReleasePin: ReleasePin = releasePin(() => globalThis.localStorage);

/** The site's tabs, through a BroadcastChannel; without one, a switch is not told (each tab still keeps its release). */
export function browserReleaseChannel(): ReleaseChannel {
  let channel: BroadcastChannel | undefined;
  try {
    channel = new BroadcastChannel(CHANNEL_NAME);
  } catch {
    channel = undefined;
  }
  return {
    announce(release): void {
      try {
        channel?.postMessage({ type: 'switched', release });
      } catch {
        // Not told: the other tabs see the pin at their next check.
      }
    },
    onSwitch(listener): () => void {
      if (channel === undefined) return () => undefined;
      const target = channel;
      const handler = (event: MessageEvent<unknown>) => {
        const data = event.data as { type?: unknown; release?: unknown } | null;
        if (data?.type === 'switched' && typeof data.release === 'string' && RELEASE_ID_PATTERN.test(data.release)) listener(data.release);
      };
      target.addEventListener('message', handler);
      return () => target.removeEventListener('message', handler);
    },
  };
}

export const browserReleasePage: ReleasePage = {
  online: () => navigator.onLine !== false,

  onOnline(listener): () => void {
    window.addEventListener('online', listener);
    return () => window.removeEventListener('online', listener);
  },

  whenIdle(): Promise<void> {
    return new Promise((resolve) => {
      const idle = () => {
        // WebKit has no requestIdleCallback: there the pause alone decides.
        if (typeof requestIdleCallback === 'function') requestIdleCallback(() => resolve(), { timeout: IDLE_TIMEOUT_MS });
        else resolve();
      };
      const settle = () => setTimeout(idle, SETTLE_MS);
      if (document.readyState === 'complete') settle();
      else window.addEventListener('load', settle, { once: true });
    });
  },

  reload(): void {
    location.reload();
  },
};
