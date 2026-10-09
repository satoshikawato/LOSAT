// Ownership of temporary storage (plan §5.6, design §7.2). Every working session (one
// tab's Data worker) keeps its data under tmp/<session-token>/ and holds the Web Lock
// `losat-web:tmp:<session-token>` for as long as the worker lives. It takes the lock before
// it creates the directory. At start-up, a session removes only the directories under
// tmp/ whose lock nobody holds. Running, stopped (frozen or in the background) and
// older-release tabs keep their locks, so their data stays. Time stamps are never used,
// and without Web Locks nothing is removed.
import type { CleanupState } from '../../ports/data';
import { isStorageFull, MEMORY_FULL_MESSAGE, MEMORY_RESULTS_CAPACITY_BYTES, type BlockStore } from './block-store';
import { MemoryBlockStore } from './memory-block-store';

/** The OPFS directory of the working sessions. Shared by every release. */
export const TMP_DIRECTORY = 'tmp';
/** Lock name prefix; shared by every release so that an older tab keeps its data. */
export const SESSION_LOCK_PREFIX = 'losat-web:tmp:';

/** The part of the Web Locks API (`navigator.locks`) that sessions use. */
export interface SessionLocks {
  request(
    name: string,
    options: { readonly ifAvailable?: boolean },
    callback: (lock: object | null) => Promise<void>,
  ): Promise<unknown>;
}

/** The directories of the working sessions (tmp/ in OPFS). */
export interface SessionSpace {
  list(): Promise<readonly string[]>;
  remove(token: string): Promise<void>;
}

/** Takes the session's lock and keeps it until the worker ends. Resolves once it is held. */
export function holdSessionLock(locks: SessionLocks, token: string): Promise<void> {
  return new Promise((resolve, reject) => {
    locks
      .request(SESSION_LOCK_PREFIX + token, {}, () => {
        resolve();
        return new Promise<void>(() => {});
      })
      .catch(reject);
  });
}

/**
 * Removes the directory of every other session whose lock is free, holding that lock
 * while it removes so that two tabs never remove the same directory. Returns how many
 * directories it removed. A directory that cannot be removed (still in use) is left.
 */
export async function reclaimAbandonedSessions(
  space: SessionSpace,
  locks: SessionLocks,
  ownToken: string,
): Promise<number> {
  let removed = 0;
  for (const token of await space.list()) {
    if (token === ownToken) continue;
    await locks.request(SESSION_LOCK_PREFIX + token, { ifAvailable: true }, async (lock) => {
      if (lock === null) return;
      try {
        await space.remove(token);
        removed++;
      } catch {
        // In use or already gone: leave it for a later start-up.
      }
    });
  }
  return removed;
}

/** OPFS as a session uses it; implemented by opfs-block-store.ts. */
export interface OpfsAccess {
  /** The session directories under tmp/. Rejects when this browser has no usable OPFS. */
  space(): Promise<SessionSpace>;
  /** Opens tmp/<token>/ as a block store after checking it; rejects when it cannot be used. */
  open(token: string): Promise<BlockStore>;
}

export interface DataSessionEnv {
  readonly token: string;
  /** `navigator.locks`, or undefined where the browser has no Web Locks. */
  readonly locks: SessionLocks | undefined;
  readonly opfs: OpfsAccess;
}

export interface DataSession {
  readonly store: BlockStore;
  /** Why the session uses memory instead of OPFS. */
  readonly fallbackReason?: string;
  /** Settles when the removal of abandoned sessions has finished. */
  readonly cleanup: Promise<CleanupState>;
}

/**
 * Starts a working session. It holds its lock, then chooses the storage by trying OPFS
 * (never by browser name), and removes abandoned sessions in the background. A tab keeps
 * its data in OPFS only while it holds its lock; without one (no Web Locks API, or a
 * request that failed) other tabs could not tell that the data is in use, so the tab
 * keeps it in memory. If OPFS cannot be opened, for example because abandoned sessions
 * fill the quota, the session removes them first and tries once more.
 */
export async function startDataSession(env: DataSessionEnv): Promise<DataSession> {
  const { locks, token } = env;
  if (locks === undefined) {
    return memorySession(
      'this browser has no Web Locks API, which keeps the temporary files of a tab safe from other tabs',
      { state: 'unavailable', reason: 'this browser has no Web Locks API' },
    );
  }
  try {
    await holdSessionLock(locks, token);
  } catch (error) {
    return memorySession(`the Web Lock of this tab could not be taken (${messageOf(error)})`, {
      state: 'unavailable',
      reason: 'this tab holds no Web Lock',
    });
  }
  let space: SessionSpace;
  try {
    space = await env.opfs.space();
  } catch (error) {
    return memorySession(messageOf(error), { state: 'done', removedSessions: 0 });
  }
  const reclaim = () =>
    reclaimAbandonedSessions(space, locks, token).then(
      (removedSessions): CleanupState => ({ state: 'done', removedSessions }),
      (error: unknown): CleanupState => ({ state: 'unavailable', reason: messageOf(error) }),
    );
  try {
    return { store: await env.opfs.open(token), cleanup: reclaim() };
  } catch {
    const cleanup = await reclaim();
    try {
      return { store: await env.opfs.open(token), cleanup: Promise.resolve(cleanup) };
    } catch (error) {
      return memorySession(openFailure(error), cleanup);
    }
  }
}

function memorySession(fallbackReason: string, cleanup: CleanupState): DataSession {
  const store = new MemoryBlockStore({ capacityBytes: MEMORY_RESULTS_CAPACITY_BYTES, fullMessage: MEMORY_FULL_MESSAGE });
  return { store, fallbackReason, cleanup: Promise.resolve(cleanup) };
}

function openFailure(error: unknown): string {
  return isStorageFull(error) ? 'the browser storage of this site is full' : messageOf(error);
}

function messageOf(error: unknown): string {
  return error instanceof Error ? error.message : String(error);
}
