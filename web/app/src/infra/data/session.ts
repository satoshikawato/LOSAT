// Ownership of temporary storage (plan §5.6, design §7.2). Every working session (one
// tab's Data worker) keeps its data under tmp/<session-token>/ and holds the Web Lock
// `losat-web:tmp:<session-token>` for as long as the worker lives. It takes the lock before
// it creates the directory. At start-up, a session removes only the directories under
// tmp/ whose lock nobody holds. Running, stopped (frozen or in the background) and
// older-release tabs keep their locks, so their data stays. Time stamps are never used,
// and without Web Locks nothing is removed.
import { MEMORY_RESULTS_CAPACITY_BYTES, type CleanupState } from '../../ports/data';
import { isStorageFull, MEMORY_FULL_MESSAGE, StorageFullError, type BlockStore } from './block-store';
import { MemoryBlockStore } from './memory-block-store';

/** The OPFS directory of the working sessions. Shared by every release. */
export const TMP_DIRECTORY = 'tmp';
/** Lock name prefix; shared by every release so that an older tab keeps its data. */
export const SESSION_LOCK_PREFIX = 'losat-web:tmp:';

/** The part of the Web Locks API (`navigator.locks`) that sessions and saved-counts stores use. */
export interface SessionLocks {
  request(
    name: string,
    options: { readonly ifAvailable?: boolean; readonly mode?: 'exclusive' | 'shared' },
    callback: (lock: object | null) => Promise<void>,
  ): Promise<unknown>;
  query?(): Promise<{
    readonly held?: ReadonlyArray<{ readonly name?: string; readonly mode?: 'exclusive' | 'shared' }>;
  }>;
}

/** The directories of the working sessions (tmp/ in OPFS). */
export interface SessionSpace {
  list(): Promise<readonly string[]>;
  remove(token: string): Promise<void>;
}

/**
 * Takes the session's lock and keeps it until the worker ends. Resolves once it is held. A
 * SecurityError is the browser refusing the site's storage, and rejects as StorageRefusedError:
 * where the site's data is blocked, Chromium and Firefox refuse this lock before any storage call.
 */
export function holdSessionLock(locks: SessionLocks, token: string): Promise<void> {
  return new Promise((resolve, reject) => {
    locks
      .request(SESSION_LOCK_PREFIX + token, {}, () => {
        resolve();
        return new Promise<void>(() => {});
      })
      .catch((error: unknown) => reject(errorName(error) === 'SecurityError' ? new StorageRefusedError(error) : error));
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
 * The browser refused the site's storage from the start: the tab's lock was refused with
 * SecurityError (a site whose data is blocked, in Chromium and Firefox), or the first
 * `getDirectory()` of a session (`space()`) rejected with SecurityError (Firefox's private
 * windows, a site whose data is blocked). Only a SecurityError is one: WebKit's UnknownError comes
 * from its ephemeral sessions and from a failing storage alike (opfs-block-store.ts). A failure
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
 * "Results are kept in this tab's memory (at most 512 MB), because <reason>." The browser's own
 * message is never shown: Firefox's private window said
 * "Security error when calling GetDirectory", and WebKit's "The operation failed for an unknown
 * transient reason (e.g. out of memory)", which reads as if memory ran out (line WQA's sweep, D5,
 * 2026-11-19). Each reason says what happened at the step that failed; the browser's error, its
 * name and message, stays in the note's data-error attribute, for developers.
 */
/**
 * The storage was refused from the start (StorageRefusedError: the tab's lock, or the first
 * getDirectory()). It names both causes and so what helps: Chromium's private windows store files,
 * so there it is always blocked site data (review storage-r2 of line WQA's sweep, D5, L1, 2026-11-19).
 */
export const REFUSED_REASON = "the browser does not let this site store files (a private window, or site data blocked in the browser's settings)";
/** QuotaExceededError (or StorageFullError), at any step. */
export const FULL_REASON = 'the browser storage of this site is full';
/** No OPFS, no writing from a worker, or no Web Locks (without which no other tab could tell the data is in use). */
export const UNSUPPORTED_REASON = 'LOSAT cannot store files in this browser';
/** The tab's lock was not given (other than SecurityError): the storage was never tried, and no other tab was asked. */
export const UNRESERVED_REASON = 'the browser could not reserve storage for this tab';
/**
 * `space()` failed otherwise: the root or tmp/ could not be opened. The root's lasting UnknownError
 * comes here: WebKit gives it for its ephemeral sessions (Safari's private windows) and for a
 * storage that fails in a normal window, and Firefox for the latter, so the sentence is true of both
 * and says only that a private window may have none (review storage-r2 of line WQA's sweep, D5, M1,
 * 2026-11-19: a normal window had been told it was a private one).
 */
export const UNOPENED_REASON = 'the browser storage of this site could not be opened (a private window may have none)';
/** `open()` failed after `space()` had worked: making the session's directory, or its check. */
export const BROKEN_REASON = 'the browser storage of this site failed when LOSAT used it';

/* The cleanup's note: "Temporary data of closed tabs is not removed: <reason>." */
/** No Web Locks: whether another tab holds a directory cannot be known. */
export const NO_LOCKS_CLEANUP = 'this browser cannot tell whether other tabs still use it';
/** The lock of another tab's directory could not be requested (the cleanup asks for each one). */
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
  /**
   * The name of the session's lock, while the worker holds it (until it ends); undefined without one. The page learns
   * from its release that the worker ended (data-worker/gateway.ts; line WQA's sweep, D2).
   */
  readonly lock?: string;
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
    // Refused (StorageRefusedError), as space() can be; any other failure leaves the storage untried.
    // Either way no other tab was asked, so nothing is said of them: the page said "the browser could
    // not tell whether other tabs still use it" where the browser had refused this site's storage
    // (review of line WQA's sweep, D5, M2, 2026-11-19).
    return memorySession(fallbackOf(error, UNRESERVED_REASON), STORAGE_UNREACHED);
  }
  const lock = SESSION_LOCK_PREFIX + token;
  let space: SessionSpace;
  try {
    space = await env.opfs.space();
  } catch (error) {
    return memorySession(fallbackOf(error, UNOPENED_REASON), STORAGE_UNREACHED, lock);
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
    return { store: await env.opfs.open(token), lock, cleanup: reclaim() };
  } catch {
    const cleanup = await reclaim();
    try {
      return { store: await env.opfs.open(token), lock, cleanup: Promise.resolve(cleanup) };
    } catch (error) {
      return memorySession(fallbackOf(error, BROKEN_REASON), cleanup, lock);
    }
  }
}

interface Fallback {
  readonly reason: string;
  readonly error?: string;
}

/** The cleanup of a session that never reached the storage: it removes nothing, and the memory note says why. */
const STORAGE_UNREACHED: CleanupState = { state: 'done', removedSessions: 0 };

/** A session in memory; `lock` when the worker holds its lock all the same (its release tells the page it ended). */
function memorySession(fallback: Fallback, cleanup: CleanupState, lock?: string): DataSession {
  const store = new MemoryBlockStore({ capacityBytes: MEMORY_RESULTS_CAPACITY_BYTES, fullMessage: MEMORY_FULL_MESSAGE });
  return {
    store,
    fallbackReason: fallback.reason,
    ...(fallback.error === undefined ? {} : { fallbackError: fallback.error }),
    ...(lock === undefined ? {} : { lock }),
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
