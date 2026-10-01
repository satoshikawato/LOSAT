import { describe, expect, it } from 'vitest';
import { StorageFullError } from '../../src/infra/data/block-store';
import { MemoryBlockStore } from '../../src/infra/data/memory-block-store';
import {
  holdSessionLock,
  reclaimAbandonedSessions,
  SESSION_LOCK_PREFIX,
  startDataSession,
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
    expect(session.fallbackReason).toBe('the browser storage of this site is full');
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
    expect(session.fallbackReason).toMatch(/^this browser has no Web Locks API/);
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: 'this browser has no Web Locks API' });
    expect(space.tokens).toEqual(['other']);
  });

  it('keeps its data in memory when its lock cannot be taken, so no other tab can remove it', async () => {
    let opened = false;
    const session = await startDataSession({
      token: 'own',
      locks: {
        request: async () => {
          throw new Error('SecurityError');
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
    expect(session.fallbackReason).toBe('the Web Lock of this tab could not be taken (SecurityError)');
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: 'this tab holds no Web Lock' });
  });

  it('falls back to memory with the reason when OPFS cannot be used', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      opfs: {
        space: async () => {
          throw new Error('this browser has no Origin Private File System');
        },
        open: async () => new MemoryBlockStore(),
      },
    });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe('this browser has no Origin Private File System');
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 0 });
  });
});
