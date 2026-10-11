// Which release caches the service worker may remove (plan §5.9: an older release's assets
// stay while a tab uses them). Every page holds the shared Web Lock of its release
// (protocol.ts RELEASE_LOCK_PREFIX) for its life, the way a Data worker holds the lock of its
// tmp/ directory (design §7.2, S10 decision 8). A release's cache goes only when nobody holds
// or waits for its lock, it is not a release the asking page names (the pinned one and its
// own), and no preparation of it is running. Without the Web Locks API nothing is removed:
// time stamps are never used.
import { RELEASE_CACHE_PREFIX, RELEASE_LOCK_PREFIX } from './protocol';

/** What `navigator.locks.query()` returns, as far as the rule reads it. */
export interface LockSnapshot {
  readonly held?: ReadonlyArray<{ readonly name?: string }>;
  readonly pending?: ReadonlyArray<{ readonly name?: string }>;
}

/** The names of the release caches to remove; none when `locks` is undefined (no Web Locks). */
export function releaseCachesToRemove(
  cacheNames: readonly string[],
  locks: LockSnapshot | undefined,
  keep: Iterable<string>,
  preparing: Iterable<string>,
): string[] {
  if (locks === undefined) return [];
  const used = new Set<string>([...keep, ...preparing]);
  for (const lock of [...(locks.held ?? []), ...(locks.pending ?? [])]) {
    if (lock.name?.startsWith(RELEASE_LOCK_PREFIX)) used.add(lock.name.slice(RELEASE_LOCK_PREFIX.length));
  }
  return cacheNames.filter((name) => name.startsWith(RELEASE_CACHE_PREFIX) && !used.has(name.slice(RELEASE_CACHE_PREFIX.length)));
}
