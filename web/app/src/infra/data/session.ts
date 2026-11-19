// Ownership of temporary storage (plan §5.6, design §7.2). Every working session (one
// tab's Data worker) keeps its data under tmp/<session-token>/ and holds the Web Lock
// `losat-web:tmp:<session-token>` for as long as the worker lives. It takes the lock before
// it creates the directory. At start-up, a session removes only the directories under
// tmp/ whose lock nobody holds. Running, stopped (frozen or in the background) and
// older-release tabs keep their locks, so their data stays. Time stamps are never used,
// and without Web Locks nothing is removed.
import type { CleanupState } from '../../ports/data';
import { isStorageFull, MEMORY_FULL_MESSAGE, MEMORY_RESULTS_CAPACITY_BYTES, StorageFullError, type BlockStore } from './block-store';
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
  /**
   * The session directories under tmp/. Rejects when this browser has no usable OPFS: with
   * StorageUnsupportedError when it lacks the API, with StorageRefusedError when the browser
   * refuses the site's storage.
   */
  space(): Promise<SessionSpace>;
  /**
   * Opens tmp/<token>/ as a block store after checking it; rejects when it cannot be used
   * (StorageUnsupportedError when this browser cannot write the files from a worker).
   */
  open(token: string): Promise<BlockStore>;
}

/**
 * The browser refused the site's storage from the start: the first `getDirectory()` of a session
 * (`space()`) rejected with SecurityError (Firefox's private windows, a site whose data is blocked),
 * or with UnknownError twice (WebKit's ephemeral sessions, as Safari's private windows). A failure
 * after the storage worked is never one. The browser's error is the cause.
 */
export class StorageRefusedError extends Error {
  constructor(cause: unknown) {
    super('the browser refused the storage of this site', { cause });
    this.name = 'StorageRefusedError';
  }
}

/** This browser lacks a part of OPFS or Web Locks that LOSAT stores files with. */
export class StorageUnsupportedError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'StorageUnsupportedError';
  }
}

/** The session directories under tmp/ could not be listed (the cleanup); the browser's error is the cause. */
export class SessionsUnreadableError extends Error {
  constructor(cause: unknown) {
    super('the session directories could not be listed', { cause });
    this.name = 'SessionsUnreadableError';
  }
}

/*
 * Why a session keeps its data in memory, in LOSAT's words: the note under the page says
 * "Results are kept in this tab's memory only and are lost when it is closed or reloaded,
 * because <reason>." The browser's own message is never shown: Firefox's private window said
 * "Security error when calling GetDirectory", and WebKit's "The operation failed for an unknown
 * transient reason (e.g. out of memory)", which reads as if memory ran out (line WQA's sweep, D5,
 * 2026-11-19). Each reason says what happened at the step that failed; the browser's error, its
 * name and message, stays in the note's data-error attribute, for developers.
 */
/** The storage was refused from the start (StorageRefusedError), or the tab's lock with SecurityError. */
export const REFUSED_REASON = 'the browser does not let this site store files here, as in a private window';
/** QuotaExceededError (or StorageFullError), at any step. */
export const FULL_REASON = 'the browser storage of this site is full';
/** No OPFS, no writing from a worker, or no Web Locks (without which no other tab could tell the data is in use). */
export const UNSUPPORTED_REASON = 'LOSAT cannot store files in this browser';
/** The tab's lock was not given (other than SecurityError): the storage was never tried. */
export const UNRESERVED_REASON = 'the browser could not reserve storage for this tab';
/** `space()` failed otherwise: the root or tmp/ could not be opened. */
export const UNOPENED_REASON = 'the browser storage of this site could not be opened';
/** `open()` failed after `space()` had worked: making the session's directory, or its check. */
export const BROKEN_REASON = 'the browser storage of this site failed when LOSAT used it';

/* The cleanup's note: "Temporary data of closed tabs is not removed: <reason>." */
/** No Web Locks: whether another tab holds a directory cannot be known. */
export const NO_LOCKS_CLEANUP = 'this browser cannot tell whether other tabs still use it';
/** The tab's lock, or the lock of another directory, could not be requested. */
export const LOCK_FAILED_CLEANUP = 'the browser could not tell whether other tabs still use it';
/** tmp/ could not be listed. */
export const UNREADABLE_CLEANUP = 'the browser storage of this site could not be read';

export interface DataSessionEnv {
  readonly token: string;
  /** `navigator.locks`, or undefined where the browser has no Web Locks. */
  readonly locks: SessionLocks | undefined;
  readonly opfs: OpfsAccess;
}

export interface DataSession {
  readonly store: BlockStore;
  /** Why the session uses memory instead of OPFS, in LOSAT's words (REFUSED_REASON…). */
  readonly fallbackReason?: string;
  /** The error behind fallbackReason, "<name>: <message>" (the browser's, under LOSAT's wrappers), for developers. */
  readonly fallbackError?: string;
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
    const missing = new StorageUnsupportedError('this browser has no Web Locks API');
    return memorySession(
      { reason: UNSUPPORTED_REASON, error: detail(missing) },
      { state: 'unavailable', reason: NO_LOCKS_CLEANUP, error: detail(missing) },
    );
  }
  try {
    await holdSessionLock(locks, token);
  } catch (error) {
    // A SecurityError is the browser's refusal of the site's storage; any other failure leaves the storage untried.
    const reason = errorName(error) === 'SecurityError' ? REFUSED_REASON : UNRESERVED_REASON;
    return memorySession({ reason, error: detail(error) }, { state: 'unavailable', reason: LOCK_FAILED_CLEANUP, error: detail(error) });
  }
  let space: SessionSpace;
  try {
    space = await env.opfs.space();
  } catch (error) {
    return memorySession(fallbackOf(error, UNOPENED_REASON), { state: 'done', removedSessions: 0 });
  }
  const reclaim = () =>
    reclaimAbandonedSessions(readable(space), locks, token).then(
      (removedSessions): CleanupState => ({ state: 'done', removedSessions }),
      (error: unknown): CleanupState => ({
        state: 'unavailable',
        reason: error instanceof SessionsUnreadableError ? UNREADABLE_CLEANUP : LOCK_FAILED_CLEANUP,
        error: detail(error),
      }),
    );
  try {
    return { store: await env.opfs.open(token), cleanup: reclaim() };
  } catch {
    const cleanup = await reclaim();
    try {
      return { store: await env.opfs.open(token), cleanup: Promise.resolve(cleanup) };
    } catch (error) {
      return memorySession(fallbackOf(error, BROKEN_REASON), cleanup);
    }
  }
}

interface Fallback {
  readonly reason: string;
  readonly error?: string;
}

function memorySession(fallback: Fallback, cleanup: CleanupState): DataSession {
  const store = new MemoryBlockStore({ capacityBytes: MEMORY_RESULTS_CAPACITY_BYTES, fullMessage: MEMORY_FULL_MESSAGE });
  return {
    store,
    fallbackReason: fallback.reason,
    ...(fallback.error === undefined ? {} : { fallbackError: fallback.error }),
    cleanup: Promise.resolve(cleanup),
  };
}

/** The reason for a failure of OPFS; `otherwise` is what any other failure at that step means. */
function fallbackOf(error: unknown, otherwise: string): Fallback {
  const reason =
    error instanceof StorageRefusedError
      ? REFUSED_REASON
      : error instanceof StorageUnsupportedError
        ? UNSUPPORTED_REASON
        : isStorageFull(error)
          ? FULL_REASON
          : otherwise;
  return { reason, error: detail(error) };
}

/** The session space whose listing failures say so (the cleanup tells them from lock failures). */
function readable(space: SessionSpace): SessionSpace {
  return {
    list: () =>
      space.list().catch((error: unknown) => {
        throw new SessionsUnreadableError(error);
      }),
    remove: (token) => space.remove(token),
  };
}

/**
 * "<name>: <message>" of the error a developer needs: the browser's, under LOSAT's wrappers
 * (StorageRefusedError, StorageFullError, SessionsUnreadableError keep it as their cause).
 */
export function detail(error: unknown): string {
  let inner = error;
  while (inner instanceof Error && inner.cause !== undefined && isWrapper(inner)) inner = inner.cause;
  if (inner instanceof Error || (typeof inner === 'object' && inner !== null)) {
    const { name, message } = inner as { name?: unknown; message?: unknown };
    const text = typeof message === 'string' ? message : '';
    const label = typeof name === 'string' && name !== '' ? name : 'Error';
    return text === '' ? label : `${label}: ${text}`;
  }
  return String(inner);
}

function isWrapper(error: Error): boolean {
  return error instanceof StorageRefusedError || error instanceof StorageFullError || error instanceof SessionsUnreadableError;
}

function errorName(error: unknown): unknown {
  return typeof error === 'object' && error !== null ? (error as { name?: unknown }).name : undefined;
}
