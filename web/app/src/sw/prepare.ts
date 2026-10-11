// The explicit preparation of a release for offline use (plan §5.9, design §14), run by the
// service worker when a page asks for it. It reads the site index's shell (the site-level
// files that start the app) and the release's ReleaseManifest, downloads every file they
// list, and checks each one's size and SHA-256. Only when every file matched does it store
// them: the release's files in `losat-web-release-<id>`, the ReleaseManifest last (it marks
// the cache complete), then the shell in `losat-web-shell`. A failure before that stores
// nothing, and a failure while storing puts the caches back as they were, so a release is
// either stored whole or not at all. Responses are stored as the server sent them, status
// and headers included: COOP, COEP, CORP and the CSP must survive for `crossOriginIsolated`.
// Nothing but these listed files ever enters a cache; research data is never a network
// response (S15 hand-over).
import {
  MANIFEST_NAME,
  SHELL_CACHE,
  SITE_INDEX_URL,
  releaseCacheName,
  releasePath,
  type PrepareProgress,
} from './protocol';

/** A file that a verified list names: the site index's shell, or a ReleaseManifest's assets. */
export interface ListedFile {
  /** A site path for the shell (`/`, `/boot.js`); a path relative to the release directory for an asset. */
  readonly path: string;
  readonly size: number;
  /** Lower-case hex SHA-256 of the bytes. */
  readonly sha256: string;
}

/** The part of a Cache that the worker uses; the browser's Cache is one. */
export interface CacheLike {
  match(key: string): Promise<Response | undefined>;
  put(key: string, response: Response): Promise<void>;
  delete(key: string): Promise<boolean>;
  keys(): Promise<ReadonlyArray<{ readonly url: string }>>;
}

/** The part of CacheStorage that the worker uses; the browser's `caches` is one. */
export interface CacheStorageLike {
  has(name: string): Promise<boolean>;
  open(name: string): Promise<CacheLike>;
  delete(name: string): Promise<boolean>;
  keys(): Promise<readonly string[]>;
}

export interface PrepareDeps {
  /** The site's origin, for the cache keys. */
  readonly origin: string;
  readonly fetch: (url: string, init: RequestInit) => Promise<Response>;
  readonly caches: CacheStorageLike;
  /** Lower-case hex SHA-256 (crypto.subtle in the worker). */
  readonly sha256: (bytes: ArrayBuffer) => Promise<string>;
  /** The shell of a parsed `/releases.json`; throws when the index is not valid. */
  readonly readShell: (index: unknown) => readonly ListedFile[];
  /** The assets of the parsed ReleaseManifest of `release`; throws when it is not valid or names another release. */
  readonly readManifest: (manifest: unknown, release: string) => readonly ListedFile[];
  /** Downloads at a time. */
  readonly concurrency?: number;
  /** The preparation fails when no file arrives for this long. */
  readonly stallMs?: number;
}

export type PrepareOutcome =
  | { readonly ok: true; readonly progress: PrepareProgress }
  | { readonly ok: false; readonly reason: string };

const DEFAULT_CONCURRENCY = 4;
const DEFAULT_STALL_MS = 120_000;

/** A downloaded file whose bytes matched, waiting to be stored. */
interface Verified {
  readonly key: string;
  readonly response: Response;
}

const message = (error: unknown): string => (error instanceof Error ? error.message : String(error));

/** Whether the cache of `release` holds its ReleaseManifest, which is stored last. */
export async function isComplete(caches: CacheStorageLike, origin: string, release: string): Promise<boolean> {
  const name = releaseCacheName(release);
  if (!(await caches.has(name))) return false;
  return (await (await caches.open(name)).match(origin + releasePath(release, MANIFEST_NAME))) !== undefined;
}

async function hasShell(caches: CacheStorageLike, origin: string): Promise<boolean> {
  if (!(await caches.has(SHELL_CACHE))) return false;
  return (await (await caches.open(SHELL_CACHE)).match(origin + '/')) !== undefined;
}

// The reasons below reach the screen after "could not be downloaded and checked:" or "could not
// be prepared:" (UpdateNotice.vue, OfflineNotice.vue): each says in plain words what failed and
// why, and ends with the file's path where one file failed (screen review 1 L5).

/** The site index, as the screen names it. */
const SITE_INDEX = "the site's list of releases";
/** A ReleaseManifest, as the screen names it. */
const FILE_LIST = "the release's file list";

/** A response that may be stored: 200 from this origin, not redirected. `what` names the file ("a file (<path>)"). */
function check(response: Response, what: string): void {
  if (response.redirected) throw new Error(`the site redirected the request for ${what}`);
  if (response.status !== 200) throw new Error(`the site answered ${response.status} for ${what}`);
  if (response.type !== 'basic' && response.type !== 'default') throw new Error(`${what} did not come from this site`);
}

/** Reads a JSON list and checks it with `read`, which throws when it is not valid. */
async function fetchJson<T>(
  deps: PrepareDeps,
  url: string,
  init: RequestInit,
  what: string,
  read: (value: unknown) => T,
): Promise<{ value: T; response: Response }> {
  let response: Response;
  try {
    response = await deps.fetch(url, { credentials: 'same-origin', redirect: 'error', ...init });
  } catch (error) {
    throw new Error(`${what} could not be read (${message(error)})`);
  }
  check(response, what);
  const kept = response.clone();
  try {
    return { value: read(JSON.parse(await response.text()) as unknown), response: kept };
  } catch (error) {
    throw new Error(`${what} is not valid (${message(error)})`);
  }
}

/**
 * Prepares `release`: stores its files and the shell after checking every one, or stores
 * nothing. A release already stored whole is not downloaded again; only the shell is read
 * again (a new deployment may have changed it), and when that fails the stored shell stays.
 */
export async function prepareRelease(
  deps: PrepareDeps,
  release: string,
  onProgress: (progress: PrepareProgress) => void,
): Promise<PrepareOutcome> {
  const { origin, caches } = deps;
  const complete = await isComplete(caches, origin, release);
  const keepStored = async (reason: string): Promise<PrepareOutcome> =>
    complete && (await hasShell(caches, origin))
      ? { ok: true, progress: { files: 0, filesTotal: 0, bytes: 0, bytesTotal: 0 } }
      : { ok: false, reason };

  let shell: readonly ListedFile[];
  try {
    shell = (await fetchJson(deps, SITE_INDEX_URL, { cache: 'no-store' }, SITE_INDEX, deps.readShell)).value;
  } catch (error) {
    return keepStored(message(error));
  }
  if (!shell.some((file) => file.path === '/')) return keepStored(`${SITE_INDEX} names no start page`);

  let assets: readonly ListedFile[] = [];
  let manifest: Response | undefined;
  if (!complete) {
    try {
      const url = origin + releasePath(release, MANIFEST_NAME);
      // A 404 under /r/ may be the HTTP cache's (the rules make it immutable there): a failed
      // read is tried once more past the HTTP cache before the release counts as not served.
      const readManifest = (value: unknown) => deps.readManifest(value, release);
      const read = await fetchJson(deps, url, { cache: 'no-cache' }, FILE_LIST, readManifest).catch(() =>
        fetchJson(deps, url, { cache: 'reload' }, FILE_LIST, readManifest),
      );
      assets = read.value;
      manifest = read.response;
    } catch (error) {
      return { ok: false, reason: message(error) };
    }
  }

  const files = [
    ...shell.map((file) => ({ file, key: origin + file.path, shell: true })),
    ...assets.map((file) => ({ file, key: origin + releasePath(release, file.path), shell: false })),
  ];
  const progress = { files: 0, filesTotal: files.length, bytes: 0, bytesTotal: files.reduce((sum, { file }) => sum + file.size, 0) };
  onProgress({ ...progress });

  const downloaded = await downloadAll(deps, files, (size) => {
    progress.files++;
    progress.bytes += size;
    onProgress({ ...progress });
  });
  if (!downloaded.ok) return complete ? keepStored(downloaded.reason) : { ok: false, reason: downloaded.reason };

  const shellFiles = downloaded.verified.filter((_, n) => files[n]?.shell === true);
  const releaseFiles = downloaded.verified.filter((_, n) => files[n]?.shell === false);
  const name = releaseCacheName(release);
  if (!complete && manifest !== undefined) {
    try {
      // A cache without its manifest is what an interrupted preparation left: start it again.
      if (await caches.has(name)) await caches.delete(name);
      const cache = await caches.open(name);
      for (const { key, response } of releaseFiles) await cache.put(key, response);
      await cache.put(origin + releasePath(release, MANIFEST_NAME), manifest);
    } catch (error) {
      await caches.delete(name).catch(() => false);
      return { ok: false, reason: `the files could not be stored (${message(error)})` };
    }
  }
  try {
    await storeShell(caches, shellFiles);
  } catch (error) {
    if (!complete) await caches.delete(name).catch(() => false);
    const reason = `the start-up files could not be stored (${message(error)})`;
    return complete ? keepStored(reason) : { ok: false, reason };
  }
  return { ok: true, progress: { ...progress } };
}

/** Replaces the shell; on a failure puts back what was there. Entries no longer listed are removed. */
async function storeShell(caches: CacheStorageLike, files: readonly Verified[]): Promise<void> {
  const cache = await caches.open(SHELL_CACHE);
  const previous = new Map<string, Response | undefined>();
  for (const { key } of files) previous.set(key, await cache.match(key));
  try {
    for (const { key, response } of files) await cache.put(key, response);
  } catch (error) {
    for (const [key, response] of previous) {
      await (response === undefined ? cache.delete(key) : cache.put(key, response)).catch(() => undefined);
    }
    throw error;
  }
  for (const request of await cache.keys()) {
    if (!previous.has(request.url)) await cache.delete(request.url).catch(() => false);
  }
}

type Downloaded = { readonly ok: true; readonly verified: readonly Verified[] } | { readonly ok: false; readonly reason: string };

/** Downloads and checks every file, a few at a time; stops at the first failure. */
async function downloadAll(
  deps: PrepareDeps,
  files: ReadonlyArray<{ readonly file: ListedFile; readonly key: string; readonly shell: boolean }>,
  done: (size: number) => void,
): Promise<Downloaded> {
  const verified: Verified[] = new Array<Verified>(files.length);
  const controller = new AbortController();
  let failure: string | undefined;
  const stallMs = deps.stallMs ?? DEFAULT_STALL_MS;
  let stall = setTimeout(() => controller.abort(new Error('no file arrived for a while')), stallMs);
  let next = 0;
  const runner = async (): Promise<void> => {
    while (failure === undefined && next < files.length) {
      const n = next++;
      const { file, key, shell } = files[n]!;
      try {
        verified[n] = { key, response: await downloadVerified(deps, key, file, shell, controller.signal) };
      } catch (error) {
        if (failure === undefined) failure = message(error);
        controller.abort();
        return;
      }
      clearTimeout(stall);
      stall = setTimeout(() => controller.abort(new Error('no file arrived for a while')), stallMs);
      done(file.size);
    }
  };
  try {
    await Promise.all(Array.from({ length: Math.min(deps.concurrency ?? DEFAULT_CONCURRENCY, files.length) }, runner));
  } finally {
    clearTimeout(stall);
  }
  return failure === undefined ? { ok: true, verified } : { ok: false, reason: failure };
}

/**
 * Downloads one file and checks it. A mismatch, or an answer that may not be stored (a 404
 * that the HTTP cache kept, say), is tried once more past the HTTP cache (a stale entry there
 * is not the server's file); a second one fails.
 */
async function downloadVerified(deps: PrepareDeps, url: string, file: ListedFile, shell: boolean, signal: AbortSignal): Promise<Response> {
  const what = `a file (${file.path})`;
  const list = shell ? SITE_INDEX : FILE_LIST;
  let mismatch = '';
  for (const cache of [shell ? 'no-cache' : 'default', 'reload'] as const) {
    let response: Response;
    let bytes: ArrayBuffer;
    try {
      response = await deps.fetch(url, { cache, credentials: 'same-origin', redirect: 'error', signal });
    } catch (error) {
      throw new Error(`a file could not be downloaded (${file.path}: ${message(signal.reason ?? error)})`);
    }
    if (cache !== 'reload') {
      try {
        check(response, what);
      } catch (error) {
        mismatch = message(error);
        continue;
      }
    } else check(response, what);
    const kept = response.clone();
    try {
      bytes = await response.arrayBuffer();
    } catch (error) {
      throw new Error(`a file could not be downloaded (${file.path}: ${message(signal.aborted ? signal.reason : error)})`);
    }
    if (bytes.byteLength !== file.size) {
      mismatch = `a file had ${bytes.byteLength} bytes, not the ${file.size} that ${list} gives (${file.path})`;
    } else if ((await deps.sha256(bytes)) !== file.sha256) {
      mismatch = `a file did not match its checksum in ${list} (${file.path})`;
    } else {
      return kept;
    }
  }
  throw new Error(mismatch);
}
