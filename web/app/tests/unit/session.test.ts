import { describe, expect, it } from 'vitest';
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
      openOpfs: async (token) => {
        expect(locks.held.has(SESSION_LOCK_PREFIX + token)).toBe(true);
        return { store, space };
      },
    });
    expect(session.store).toBe(store);
    expect(session.fallbackReason).toBeUndefined();
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 1 });
  });

  it('removes nothing without Web Locks', async () => {
    const space = new FakeSpace(['other']);
    const session = await startDataSession({
      token: 'own',
      locks: undefined,
      openOpfs: async () => ({ store: new MemoryBlockStore(), space }),
    });
    expect(await session.cleanup).toEqual({ state: 'unavailable', reason: 'this browser has no Web Locks API' });
    expect(space.tokens).toEqual(['other']);
  });

  it('falls back to memory with the reason when OPFS cannot be used', async () => {
    const session = await startDataSession({
      token: 'own',
      locks: new FakeLocks(),
      openOpfs: async () => {
        throw new Error('this browser cannot write OPFS files from a worker (no createSyncAccessHandle)');
      },
    });
    expect(session.store.backend).toBe('memory');
    expect(session.fallbackReason).toBe('this browser cannot write OPFS files from a worker (no createSyncAccessHandle)');
    expect(await session.cleanup).toEqual({ state: 'done', removedSessions: 0 });
  });
});
