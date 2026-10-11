// What the service worker answers (plan §5.9). A release's files come from its cache once
// the cache is complete, else from the network without being stored. The start page and the
// other shell files come from the network first and from the stored shell when the network
// fails or takes too long, so that a start without a network still works. Every other request
// (`/releases.json`, `/api/*`, `/sw.js`, any method but GET, other origins, `blob:` and
// `data:` URLs) is not answered by the worker at all: the browser handles it as if there were
// no worker, and nothing of it is stored.
import { RELEASE_ID_PATTERN } from '../domain/site-index';
import { isComplete, type CacheStorageLike } from './prepare';
import { SHELL_CACHE, SW_URL, releaseCacheName } from './protocol';

export type Route =
  | { readonly kind: 'release'; readonly release: string }
  | { readonly kind: 'shell'; readonly path: string }
  | { readonly kind: 'none' };

/** A file at the top of the site: the boot page's files (`/boot.js`, `/boot.css`, `/404.html`, ...). */
const SHELL_FILE = /^\/[0-9A-Za-z][0-9A-Za-z._-]*\.(?:html|js|css)$/;
const RELEASE_FILE = /^\/r\/([^/]+)\/./;

/**
 * Which answer a request gets. The decision is made from the request alone, at once: the worker
 * must call `respondWith` in the event's own task, before it could read a cache.
 */
export function routeOf(request: { readonly method: string; readonly url: string }, origin: string): Route {
  if (request.method !== 'GET') return { kind: 'none' };
  let url: URL;
  try {
    url = new URL(request.url);
  } catch {
    return { kind: 'none' };
  }
  // A `blob:` URL has the origin of the page that made it: the scheme rules it out first.
  if ((url.protocol !== 'http:' && url.protocol !== 'https:') || url.origin !== origin) return { kind: 'none' };
  const release = RELEASE_FILE.exec(url.pathname)?.[1];
  if (release !== undefined) return RELEASE_ID_PATTERN.test(release) ? { kind: 'release', release } : { kind: 'none' };
  if (url.pathname === '/' || (SHELL_FILE.test(url.pathname) && url.pathname !== SW_URL)) return { kind: 'shell', path: url.pathname };
  return { kind: 'none' };
}

export interface ServeDeps {
  readonly origin: string;
  readonly caches: CacheStorageLike;
  readonly fetch: (request: Request) => Promise<Response>;
  /** How long the network may take for a shell file before the stored one is used. */
  readonly shellTimeoutMs?: number;
}

const DEFAULT_SHELL_TIMEOUT_MS = 3_000;

/** A file of a release: from its complete cache, else from the network (not stored). */
export async function serveRelease(deps: ServeDeps, request: Request, release: string): Promise<Response> {
  if (await isComplete(deps.caches, deps.origin, release)) {
    const stored = await (await deps.caches.open(releaseCacheName(release))).match(request.url);
    if (stored !== undefined) return stored;
  }
  return deps.fetch(request);
}

/**
 * A shell file: the network's answer, unless the network fails, answers a server error or
 * takes longer than the timeout and the shell is stored; then the stored file.
 */
export async function serveShell(deps: ServeDeps, request: Request, path: string): Promise<Response> {
  const stored = (await deps.caches.has(SHELL_CACHE))
    ? await (await deps.caches.open(SHELL_CACHE)).match(deps.origin + path)
    : undefined;
  if (stored === undefined) return deps.fetch(request);
  const network = deps.fetch(request);
  let timer: ReturnType<typeof setTimeout> | undefined;
  const timeout = new Promise<undefined>((resolve) => {
    timer = setTimeout(() => resolve(undefined), deps.shellTimeoutMs ?? DEFAULT_SHELL_TIMEOUT_MS);
  });
  try {
    const response = await Promise.race([network.catch(() => undefined), timeout]);
    return response !== undefined && response.status < 500 ? response : stored;
  } finally {
    clearTimeout(timer);
  }
}
