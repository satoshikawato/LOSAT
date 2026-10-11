// The names and messages that the site's service worker (src/sw/) and the pages of every
// release share (plan §5.9, design §14). The worker is site-level and the pages are
// release-level, so a page of an older release talks to the worker of a newer site and the
// other way round: these names are a contract between releases. Renaming a cache or a lock
// would make a new worker delete what an older tab still uses, or keep what nobody uses.

/** The worker's script and scope: site-level, not versioned. */
export const SW_URL = '/sw.js';
export const SW_SCOPE = '/';

/** The site index (src/domain/site-index.ts): read from the network only, never stored. */
export const SITE_INDEX_URL = '/releases.json';

/** The cache of one release: `losat-web-release-<release-id>`, holding exactly the files of its ReleaseManifest. */
export const RELEASE_CACHE_PREFIX = 'losat-web-release-';
/** The cache of the site-level files that start the app (the site index's `shell`). */
export const SHELL_CACHE = 'losat-web-shell';
/** The ReleaseManifest's name in a release directory; stored last, it marks a complete release cache. */
export const MANIFEST_NAME = 'release-manifest.json';
/**
 * Every page holds the shared Web Lock `losat-web:release:<release-id>` of the release it runs,
 * for its whole life; the worker removes a release's cache only when nobody holds that lock.
 */
export const RELEASE_LOCK_PREFIX = 'losat-web:release:';

/** The version of the messages below; a worker answers a message of another version with `refused`. */
export const PROTOCOL_VERSION = 1;

/** Files and bytes of a preparation: done of total. */
export interface PrepareProgress {
  readonly files: number;
  readonly filesTotal: number;
  readonly bytes: number;
  readonly bytesTotal: number;
}

/** A page's request; the reply comes on the MessagePort sent with it (none for `skip-waiting`). */
export type PageMessage =
  | { readonly v: typeof PROTOCOL_VERSION; readonly type: 'prepare'; readonly release: string }
  | { readonly v: typeof PROTOCOL_VERSION; readonly type: 'cleanup'; readonly keep: readonly string[] }
  | { readonly v: typeof PROTOCOL_VERSION; readonly type: 'skip-waiting' };

/** The worker's replies. `progress` may come several times before `prepared`. */
export type WorkerMessage =
  | { readonly type: 'progress'; readonly progress: PrepareProgress }
  | { readonly type: 'prepared'; readonly ok: true }
  | { readonly type: 'prepared'; readonly ok: false; readonly reason: string }
  | { readonly type: 'cleaned'; readonly removed: readonly string[] }
  | { readonly type: 'refused'; readonly reason: string };

/** The cache name of a release. */
export function releaseCacheName(release: string): string {
  return RELEASE_CACHE_PREFIX + release;
}

/** The path of a file of a release directory. */
export function releasePath(release: string, path: string): string {
  return `/r/${release}/${path}`;
}
