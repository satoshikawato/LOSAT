import { afterEach, describe, expect, it, vi } from 'vitest';
import { MEMORY_RESULTS_CAPACITY_BYTES, StorageFullError } from '../../src/infra/data/block-store';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { opfsAccess } from '../../src/infra/data/opfs-block-store';
import {
  BROKEN_REASON,
  detail,
  LOCK_FAILED_CLEANUP,
  NO_LOCKS_CLEANUP,
  FULL_REASON,
  holdSessionLock,
  reclaimAbandonedSessions,
  REFUSED_REASON,
  SESSION_LOCK_PREFIX,
  startDataSession,
  StorageRefusedError,
  StorageUnsupportedError,
  UNOPENED_REASON,
  UNREADABLE_CLEANUP,
  UNRESERVED_REASON,
  UNSUPPORTED_REASON,
  type OpfsAccess,
  type SessionLocks,
  type SessionSpace,
} from '../../src/infra/data/session';

/** Web Locks for one origin: exclusive locks, no queueing (the tests never wait). */
class FakeLocks implements SessionLocks {
  readonly held = new Set<string>();

  async request(name: string, options: { ifAvailable?: boolean }, callback: (lock: object | null) => Promise<void>) {
    if (this.held.has(name)) {
      if (options.ifAvailable) return callback(null);
      throw new Error('FakeLocks does not queue requests');
    }
    this.held.add(name);
    try {
      return await callback({ name });
    } finally {
      this.held.delete(name);
    }
  }
}

class FakeSpace implements SessionSpace {
  readonly removed: string[] = [];
  constructor(
    readonly tokens: string[],
    private readonly onRemove: (token: string) => void = () => undefined,
  ) {}
  async list() {
    return [...this.tokens];
  }
  async remove(token: string) {
    this.onRemove(token);
    this.tokens.splice(this.tokens.indexOf(token), 1);
    this.removed.push(token);
  }
}

describe('working session ownership', () => {
  it('holds the session lock until the worker ends', async () => {
    const locks = new FakeLocks();
    await holdSessionLock(locks, 'own');
    expect(locks.held.has(`${SESSION_LOCK_PREFIX}own`)).toBe(true);
  });

  it('removes only the sessions whose lock is free, holding their lock while it removes', async () => {
    const locks = new FakeLocks();
    await holdSessionLock(locks, 'own');
    await holdSessionLock(locks, 'running-tab');
    await holdSessionLock(locks, 'frozen-tab');
    const space = new FakeSpace(['own', 'running-tab', 'closed-tab', 'frozen-tab', 'crashed-tab'], (token) => {
      expect(locks.held.has(SESSION_LOCK_PREFIX + token)).toBe(true);
    });
    expect(await reclaimAbandonedSessions(space, locks, 'own')).toBe(2);
    expect(space.removed).toEqual(['closed-tab', 'crashed-tab']);
    expect(space.tokens).toEqual(['own', 'running-tab', 'frozen-tab']);
    expect(locks.held.has(`${SESSION_LOCK_PREFIX}closed-tab`)).toBe(false);
  });

  it('leaves a directory that cannot be removed and goes on', async () => {
    const locks = new FakeLocks();
    const space = new FakeSpace(['busy', 'gone'], (token) => {
      if (token === 'busy') throw new Error('NoModificationAllowedError');
    });
    expect(await reclaimAbandonedSessions(space, locks, 'own')).toBe(1);
    expect(space.tokens).toEqual(['busy']);
  });

  it('takes its lock before it opens its directory, then reclaims in the background', async () => {
    const locks = new FakeLocks();
    const space = new FakeSpace(['dead']);
    const store = new MemoryBlockStore();
    const session = await startDataSession({
      token: 'own',
      locks,
      opfs: {
        space: async () => space,
        open: async (token) => {
          expect(locks.held.has(SESSION_LOCK_PREFIX + token)).toBe(true);
          return store;
        },
      },
    });
    expect(session.store).toBe(store);
    expect(session.fallbackReason).toBeUndefined();
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 1 });
  });

  it('reclaims abandoned sessions first when OPFS is too full to open, then opens it', async () => {
    const space = new FakeSpace(['dead']);
    const store = new MemoryBlockStore();
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      opfs: {
        space: async () => space,
        open: async () => {
          if (space.tokens.length > 0) throw new StorageFullError();
          return store;
        },
      },
    });
    expect(session.store).toBe(store);
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 1 });
  });

  it('keeps its data in memory when OPFS stays full after the reclaim', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      opfs: {
        space: async () => new FakeSpace([]),
        open: async () => {
          throw new StorageFullError();
        },
      },
    });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe(FULL_REASON);
    expect(session.fallbackError).toBe('StorageFullError: ' + new StorageFullError().message);
    // Results kept in memory have a budget (S09).
    expect((session.store as MemoryBlockStore).capacityBytes()).toBe(MEMORY_RESULTS_CAPACITY_BYTES);
  });

  it('keeps its data in memory and removes nothing without Web Locks', async () => {
    const space = new FakeSpace(['other']);
    let opened = false;
    const session = await startDataSession({
      token: 'own',
      locks: undefined,
      opfs: {
        space: async () => space,
        open: async () => {
          opened = true;
          return new MemoryBlockStore();
        },
      },
    });
    expect(opened).toBe(false);
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe(UNSUPPORTED_REASON);
    expect(session.fallbackError).toBe('StorageUnsupportedError: this browser has no Web Locks API');
    expect(await session.cleanup).toEqual({
      state: 'unavailable',
      reason: NO_LOCKS_CLEANUP,
      error: 'StorageUnsupportedError: this browser has no Web Locks API',
    });
    expect(space.tokens).toEqual(['other']);
  });

  it('keeps its data in memory when its lock cannot be taken, so no other tab can remove it', async () => {
    let opened = false;
    const session = await startDataSession({
      token: 'own',
      locks: {
        request: async () => {
          throw new DOMException('The request was denied.', 'SecurityError');
        },
      },
      opfs: {
        space: async () => new FakeSpace([]),
        open: async () => {
          opened = true;
          return new MemoryBlockStore();
        },
      },
    });
    expect(opened).toBe(false);
    expect(session.store.backend).toBe('memory');
    // A SecurityError for the tab's lock, the first storage call, is the browser's refusal.
    expect(session.fallbackReason).toBe(REFUSED_REASON);
    expect(session.fallbackError).toBe('SecurityError: The request was denied.');
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: LOCK_FAILED_CLEANUP, error: 'SecurityError: The request was denied.' });
  });

  it('falls back to memory with the reason when OPFS cannot be used', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      opfs: {
        space: async () => {
          throw new StorageUnsupportedError('this browser has no Origin Private File System');
        },
        open: async () => new MemoryBlockStore(),
      },
    });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe(UNSUPPORTED_REASON);
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 0 });
  });
});

// Line WQA's sweep, D5 (2026-11-19): in a private window the note under the page quoted the
// browser's exception ("(The operation failed for an unknown transient reason (e.g. out of
// memory).)" in WebKit). Every reason is LOSAT's own words and true of the step that failed;
// "as in a private window" only for a refusal from the start. The browser's error, name and
// message, is kept apart for developers.
describe('why a session keeps its data in memory', () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  const opfs = (failure: { space?: unknown; open?: unknown; list?: unknown }): OpfsAccess => ({
    space: async () => {
      if (failure.space !== undefined) throw failure.space;
      return {
        list: async () => {
          if (failure.list !== undefined) throw failure.list;
          return [];
        },
        remove: async () => undefined,
      };
    },
    open: async () => {
      if (failure.open !== undefined) throw failure.open;
      return new MemoryBlockStore();
    },
  });

  const WEBKIT = new DOMException('The operation failed for an unknown transient reason (e.g. out of memory).', 'UnknownError');
  const FIREFOX = new DOMException('Security error when calling GetDirectory', 'SecurityError');
  const QUOTA = new DOMException('The quota has been exceeded.', 'QuotaExceededError');
  const at = (error: DOMException) => `${error.name}: ${error.message}`;

  it.each([
    ['refused from the start (WebKit, ephemeral)', { space: new StorageRefusedError(WEBKIT) }, REFUSED_REASON, at(WEBKIT)],
    ['refused from the start (Firefox, private window)', { space: new StorageRefusedError(FIREFOX) }, REFUSED_REASON, at(FIREFOX)],
    ['tmp/ could not be made', { space: new DOMException('Disk I/O failed.', 'UnknownError') }, UNOPENED_REASON, 'UnknownError: Disk I/O failed.'],
    ['full at the start', { space: QUOTA }, FULL_REASON, at(QUOTA)],
    ['no OPFS', { space: new StorageUnsupportedError('this browser has no Origin Private File System') }, UNSUPPORTED_REASON, 'StorageUnsupportedError: this browser has no Origin Private File System'],
    ['no writing from a worker', { open: new StorageUnsupportedError('no createSyncAccessHandle') }, UNSUPPORTED_REASON, 'StorageUnsupportedError: no createSyncAccessHandle'],
    ['full after the reclaim, the browser error kept', { open: new StorageFullError(undefined, QUOTA) }, FULL_REASON, at(QUOTA)],
    ['a failure after the storage worked', { open: WEBKIT }, BROKEN_REASON, at(WEBKIT)],
    ['the check of OPFS failed', { open: new Error('OPFS returned other bytes than were written') }, BROKEN_REASON, 'Error: OPFS returned other bytes than were written'],
  ] as const)('%s', async (_case, failure, reason, error) => {
    const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfs(failure) });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe(reason);
    expect(session.fallbackError).toBe(error);
  });

  it('a lock that the browser fails to give leaves the storage untried', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: {
        request: async () => {
          throw new DOMException('The lock manager failed.', 'AbortError');
        },
      },
      opfs: opfs({}),
    });
    expect(session.fallbackReason).toBe(UNRESERVED_REASON);
    expect(session.fallbackError).toBe('AbortError: The lock manager failed.');
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: LOCK_FAILED_CLEANUP, error: 'AbortError: The lock manager failed.' });
  });

  it('a cleanup tells a listing that failed from a lock that failed', async () => {
    const listed = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfs({ list: WEBKIT }) });
    expect(listed.fallbackReason).toBeUndefined();
    expect(await listed.cleanup).toEqual({ state: 'unavailable', reason: UNREADABLE_CLEANUP, error: at(WEBKIT) });
    const own = new FakeLocks();
    const locks: SessionLocks = {
      request: async (name, options, callback) => {
        if (options.ifAvailable) throw new DOMException('The lock manager failed.', 'AbortError');
        return own.request(name, options, callback);
      },
    };
    const space = new FakeSpace(['other']);
    const locked = await startDataSession({ token: 'own', locks, opfs: { space: async () => space, open: async () => new MemoryBlockStore() } });
    expect(await locked.cleanup).toEqual({ state: 'unavailable', reason: LOCK_FAILED_CLEANUP, error: 'AbortError: The lock manager failed.' });
  });

  /** navigator.storage whose getDirectory() answers in turn: an error to reject with, or a root. */
  function storage(...answers: unknown[]): { calls: () => number } {
    let calls = 0;
    const handle = (): unknown => ({
      kind: 'directory',
      getDirectoryHandle: async () => handle(),
      removeEntry: async () => undefined,
      entries: async function* () {},
    });
    vi.stubGlobal('navigator', {
      storage: {
        getDirectory: async () => {
          const answer = answers[Math.min(calls++, answers.length - 1)];
          if (answer !== undefined) throw answer;
          return handle();
        },
      },
    });
    return { calls: () => calls };
  }

  it('the first getDirectory() refused (SecurityError, or UnknownError twice) is a private window', async () => {
    for (const answers of [[FIREFOX], [WEBKIT, WEBKIT]]) {
      storage(...answers);
      const error = await opfsAccess.space().then(
        () => undefined,
        (rejected: unknown) => rejected,
      );
      expect(error).toBeInstanceOf(StorageRefusedError);
      expect((error as StorageRefusedError).cause).toBe(answers[0]);
      const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfsAccess });
      expect(session.fallbackReason).toBe(REFUSED_REASON);
      expect(session.fallbackError).toBe(at(answers[0]!));
    }
  });

  it('an UnknownError that passes is not a refusal, and a full storage at the root is full', async () => {
    const passing = storage(WEBKIT, undefined);
    await expect(opfsAccess.space()).resolves.toBeDefined();
    expect(passing.calls()).toBe(2);
    storage(QUOTA);
    const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfsAccess });
    expect(session.fallbackReason).toBe(FULL_REASON);
    expect(session.fallbackError).toBe(at(QUOTA));
  });

  it('getDirectory() failing after the storage worked is not a private window (open, reopen)', async () => {
    // space() gets the root; open() then meets UnknownError, or QuotaExceededError, from getDirectory().
    for (const [failure, reason] of [
      [WEBKIT, BROKEN_REASON],
      [QUOTA, FULL_REASON],
    ] as const) {
      storage(undefined, failure);
      const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfsAccess });
      expect(session.fallbackReason).toBe(reason);
      expect(session.fallbackError).toBe(at(failure));
    }
  });

  it('keeps the error a developer needs: the browser error under the wrappers, its message too', () => {
    expect(detail(new StorageRefusedError(FIREFOX))).toBe(at(FIREFOX));
    expect(detail(new StorageFullError(undefined, QUOTA))).toBe(at(QUOTA));
    expect(detail(new Error('OPFS returned other bytes than were written'))).toBe('Error: OPFS returned other bytes than were written');
    expect(detail('a string')).toBe('a string');
  });

  it('no reason quotes the browser or names its mechanisms', () => {
    const reasons = [REFUSED_REASON, FULL_REASON, UNSUPPORTED_REASON, UNRESERVED_REASON, UNOPENED_REASON, BROKEN_REASON];
    const cleanups = [NO_LOCKS_CLEANUP, LOCK_FAILED_CLEANUP, UNREADABLE_CLEANUP];
    for (const reason of [...reasons, ...cleanups]) {
      expect(reason).not.toMatch(/error|exception|transient|GetDirectory|OPFS|Web Lock|lock|API/i);
    }
  });
});
