import { afterEach, describe, expect, it, vi } from 'vitest';
import { MEMORY_RESULTS_CAPACITY_BYTES, StorageFullError } from '../../src/infra/data/block-store';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import { opfsAccess } from '../../src/infra/data/opfs-block-store';
import {
  FAILED_REASON,
  FULL_REASON,
  holdSessionLock,
  reclaimAbandonedSessions,
  REFUSED_REASON,
  SESSION_LOCK_PREFIX,
  startDataSession,
  StorageRefusedError,
  StorageUnsupportedError,
  UNREADABLE_REASON,
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
    expect(session.fallbackErrorName).toBe('StorageFullError');
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
    expect(session.fallbackErrorName).toBeUndefined();
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: 'this browser has no Web Locks API' });
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
    // A SecurityError is the browser's refusal, wherever it comes from.
    expect(session.fallbackReason).toBe(REFUSED_REASON);
    expect(session.fallbackErrorName).toBe('SecurityError');
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: 'this tab holds no Web Lock' });
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
// memory).)" in WebKit). Every reason is LOSAT's own words; the error's name is kept apart.
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

  it.each([
    ['WebKit, ephemeral (a private window)', { space: new StorageRefusedError(WEBKIT) }, REFUSED_REASON, 'UnknownError'],
    ['Firefox, private window', { space: new StorageRefusedError(FIREFOX) }, REFUSED_REASON, 'SecurityError'],
    ['a SecurityError past the root', { space: FIREFOX }, REFUSED_REASON, 'SecurityError'],
    ['tmp/ could not be made', { space: new DOMException('Disk I/O failed.', 'UnknownError') }, FAILED_REASON, 'UnknownError'],
    ['no OPFS', { space: new StorageUnsupportedError('this browser has no Origin Private File System') }, UNSUPPORTED_REASON, undefined],
    ['no writing from a worker', { open: new StorageUnsupportedError('no createSyncAccessHandle') }, UNSUPPORTED_REASON, undefined],
    ['full after the reclaim', { open: new DOMException('Quota exceeded.', 'QuotaExceededError') }, FULL_REASON, 'QuotaExceededError'],
    ['the check of OPFS failed', { open: new Error('OPFS returned other bytes than were written') }, FAILED_REASON, 'Error'],
  ] as const)('%s', async (_case, failure, reason, errorName) => {
    const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfs(failure) });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe(reason);
    expect(session.fallbackErrorName).toBe(errorName);
  });

  it('a lock that the browser fails to give is a failure, not a refusal', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: {
        request: async () => {
          throw new DOMException('The lock manager failed.', 'AbortError');
        },
      },
      opfs: opfs({}),
    });
    expect(session.fallbackReason).toBe(FAILED_REASON);
    expect(session.fallbackErrorName).toBe('AbortError');
  });

  it('a cleanup that cannot list the sessions gives LOSAT words and the name of the error', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      opfs: opfs({ list: new DOMException('The operation failed for an unknown transient reason (e.g. out of memory).', 'UnknownError') }),
    });
    expect(session.fallbackReason).toBeUndefined();
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: UNREADABLE_REASON, errorName: 'UnknownError' });
  });

  it("OPFS's root refused by the browser is a refusal, whatever the error's name", async () => {
    for (const refusal of [FIREFOX, WEBKIT]) {
      vi.stubGlobal('navigator', { storage: { getDirectory: () => Promise.reject(refusal) } });
      const error = await opfsAccess.space().then(
        () => undefined,
        (rejected: unknown) => rejected,
      );
      expect(error).toBeInstanceOf(StorageRefusedError);
      expect((error as StorageRefusedError).cause).toBe(refusal);
      const session = await startDataSession({ token: 'own', locks: new FakeLocks(), opfs: opfsAccess });
      expect(session.fallbackReason).toBe(REFUSED_REASON);
      expect(session.fallbackErrorName).toBe(refusal.name);
    }
  });

  it('no reason quotes the browser', () => {
    for (const reason of [REFUSED_REASON, FULL_REASON, UNSUPPORTED_REASON, FAILED_REASON, UNREADABLE_REASON]) {
      expect(reason).not.toMatch(/error|exception|transient|GetDirectory|OPFS|Web Lock/i);
    }
  });
});
