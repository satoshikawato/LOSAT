// Ownership of temporary storage (plan §5.6, design §7.2). Every working session (one
// tab's Data worker) keeps its data under tmp/<session-token>/ and holds the Web Lock
// `losat-web:tmp:<session-token>` for as long as the worker lives. It takes the lock before
// it creates the directory. At start-up, a session removes only the directories under
// tmp/ whose lock nobody holds. Running, stopped (frozen or in the background) and
// older-release tabs keep their locks, so their data stays. Time stamps are never used,
// and without Web Locks nothing is removed.
import type { CleanupState } from '../../ports/data';
import type { BlockStore } from './block-store';
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

export interface DataSessionEnv {
  readonly token: string;
  /** `navigator.locks`, or undefined where the browser has no Web Locks. */
  readonly locks: SessionLocks | undefined;
  /** Opens OPFS for this session; rejects when OPFS cannot be used here. */
  readonly openOpfs: (token: string) => Promise<{ readonly store: BlockStore; readonly space: SessionSpace }>;
}

export interface DataSession {
  readonly store: BlockStore;
  /** Why the session uses memory instead of OPFS. */
  readonly fallbackReason?: string;
  /** Settles when the removal of abandoned sessions has finished. */
  readonly cleanup: Promise<CleanupState>;
}

/**
 * Starts a working session: holds its lock, then chooses the storage by trying OPFS
 * (never by browser name) and falls back to memory, then removes abandoned sessions in
 * the background.
 */
export async function startDataSession(env: DataSessionEnv): Promise<DataSession> {
  if (env.locks !== undefined) await holdSessionLock(env.locks, env.token);
  let opfs: { readonly store: BlockStore; readonly space: SessionSpace } | undefined;
  let fallbackReason: string | undefined;
  try {
    opfs = await env.openOpfs(env.token);
  } catch (error) {
    fallbackReason = error instanceof Error ? error.message : String(error);
  }
  const locks = env.locks;
  let cleanup: Promise<CleanupState>;
  if (locks === undefined) {
    cleanup = Promise.resolve({ state: 'unavailable', reason: 'this browser has no Web Locks API' });
  } else if (opfs === undefined) {
    cleanup = Promise.resolve({ state: 'done', removedSessions: 0 });
  } else {
    cleanup = reclaimAbandonedSessions(opfs.space, locks, env.token).then(
      (removedSessions): CleanupState => ({ state: 'done', removedSessions }),
      (error: unknown): CleanupState => ({
        state: 'unavailable',
        reason: error instanceof Error ? error.message : String(error),
      }),
    );
  }
  return {
    store: opfs?.store ?? new MemoryBlockStore(),
    ...(fallbackReason === undefined ? {} : { fallbackReason }),
    cleanup,
  };
}
