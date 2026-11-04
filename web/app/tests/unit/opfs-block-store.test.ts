// The OPFS block store when the browser removes its directory under it (the site's data was
// cleared while the tab stayed open): a fake directory tree that fails as Chromium and Firefox do.
import { describe, expect, it } from 'vitest';
import {
  RESULTS_REMOVED_MESSAGE,
  StorageFullError,
  StorageRemovedError,
  STORAGE_REMOVED_MESSAGE,
} from '../../src/infra/data/block-store';
import { OpfsBlockStore } from '../../src/infra/data/opfs-block-store';

class FakeFile {
  readonly kind = 'file';
  bytes = new Uint8Array(0);
  constructor(private readonly parent: FakeDirectory) {}
  async createSyncAccessHandle() {
    const failure = this.parent.handleFailures.shift();
    if (failure !== undefined) throw new DOMException('fault', failure);
    return {
      write: (buffer: Uint8Array, options?: { at?: number }) => {
        const at = options?.at ?? 0;
        const next = new Uint8Array(Math.max(this.bytes.length, at + buffer.length));
        next.set(this.bytes);
        next.set(buffer, at);
        this.bytes = next;
        return buffer.length;
      },
      flush: () => undefined,
      close: () => undefined,
    };
  }
  async getFile() {
    this.parent.assertThere('NotFoundError', 'Entry not found');
    return new Blob([this.bytes as BlobPart]);
  }
}

/** A directory that fails when it, or a directory above it, was removed: the way a browser answers for each name. */
class FakeDirectory {
  readonly kind = 'directory';
  removed = false;
  /** Errors for the next openings of a file's access handle. */
  readonly handleFailures: string[] = [];
  readonly entries = new Map<string, FakeDirectory | FakeFile>();
  constructor(
    private readonly parent: FakeDirectory | undefined,
    private readonly createFailure: string,
  ) {}
  private gone(): boolean {
    return this.removed || (this.parent?.gone() ?? false);
  }
  assertThere(name: string, message: string): void {
    if (this.gone()) throw new DOMException(message, name);
  }
  async getDirectoryHandle(name: string, options?: { create?: boolean }): Promise<FakeDirectory> {
    this.assertThere(options?.create ? this.createFailure : 'NotFoundError', 'removed');
    let entry = this.entries.get(name);
    if (entry === undefined) {
      if (options?.create !== true) throw new DOMException('Entry not found', 'NotFoundError');
      entry = new FakeDirectory(this, this.createFailure);
      this.entries.set(name, entry);
    }
    return entry as FakeDirectory;
  }
  async getFileHandle(name: string, options?: { create?: boolean }): Promise<FakeFile> {
    this.assertThere(options?.create ? this.createFailure : 'NotFoundError', 'removed');
    let entry = this.entries.get(name);
    if (entry === undefined) {
      if (options?.create !== true) throw new DOMException('Entry not found', 'NotFoundError');
      entry = new FakeFile(this);
      this.entries.set(name, entry);
    }
    return entry as FakeFile;
  }
  async removeEntry(name: string): Promise<void> {
    this.assertThere('NotFoundError', 'removed');
    this.entries.delete(name);
  }
}

const asHandle = (directory: FakeDirectory) => directory as unknown as FileSystemDirectoryHandle;


for (const [browser, failure] of [
  ['Chromium', 'NotFoundError'],
  ['Firefox', 'UnknownError'],
] as const) {
  describe(`OpfsBlockStore when the browser removes its directory (${browser})`, () => {
    async function setup(reopen: (() => Promise<FileSystemDirectoryHandle>) | 'again' = 'again') {
      const top = new FakeDirectory(undefined, failure);
      const session = () => top.getDirectoryHandle('session', { create: true });
      const first = await session();
      const store = new OpfsBlockStore(asHandle(first), {
        removed: async () => !top.entries.has('session'),
        reopen: reopen === 'again' ? async () => asHandle(await session()) : reopen,
      });
      // What clearing the site's data does: the directory goes, the page and its worker stay.
      const wipe = () => {
        first.removed = true;
        top.entries.delete('session');
      };
      return { store, wipe, first };
    }

    it('makes the directory again for the next block, and says that the block before it is gone', async () => {
      const { store, wipe } = await setup();
      const before = await store.create('runs/r1/out6');
      before.append(new Uint8Array([1, 2, 3]));
      await before.seal();
      wipe();
      const writer = await store.create('runs/r2/out6');
      writer.append(new Uint8Array([4, 5]));
      await writer.seal();
      expect([...(await store.read('runs/r2/out6', 0, 2))]).toEqual([4, 5]);
      await expect(store.read('runs/r1/out6', 0, 3)).rejects.toThrow(RESULTS_REMOVED_MESSAGE);
    });

    it('says that the storage was removed, and how to go on, when the directory cannot be made again', async () => {
      const cause = new DOMException('nope', 'SecurityError');
      const { store, wipe } = await setup(() => Promise.reject(cause));
      wipe();
      const error = await store.create('runs/r1/out6').catch((e: unknown) => e);
      expect(error).toBeInstanceOf(StorageRemovedError);
      expect((error as Error).message).toBe(STORAGE_REMOVED_MESSAGE);
      expect((error as Error).message).toMatch(/Reload the page\.$/);
      expect((error as Error).cause).toBe(cause);
    });

    it('reports a full storage when that is why the directory cannot be made again', async () => {
      const { store, wipe } = await setup(() => Promise.reject(new DOMException('full', 'QuotaExceededError')));
      wipe();
      await expect(store.create('runs/r1/out6')).rejects.toBeInstanceOf(StorageFullError);
    });
  });
}

describe('OpfsBlockStore with a directory that is still there', () => {
  async function setup() {
    const top = new FakeDirectory(undefined, 'UnknownError');
    const session = await top.getDirectoryHandle('session', { create: true });
    let reopened = 0;
    const store = new OpfsBlockStore(asHandle(session), {
      removed: async () => !top.entries.has('session'),
      reopen: async () => {
        reopened++;
        return asHandle(session);
      },
    });
    return { store, session, reopened: () => reopened };
  }

  it('leaves a fault of the browser as it is: no new directory, no message about a removed storage', async () => {
    const { store, session, reopened } = await setup();
    const dir = await session.getDirectoryHandle('runs', { create: true });
    const run = await dir.getDirectoryHandle('r1', { create: true });
    run.handleFailures.push('UnknownError');
    const error = await store.create('runs/r1/out6').catch((e: unknown) => e);
    expect(error).toBeInstanceOf(DOMException);
    expect((error as Error).name).toBe('UnknownError');
    expect(reopened()).toBe(0);
    // The empty file that the attempt made is gone, so that the next attempt makes a block, not "already exists".
    expect([...run.entries.keys()]).toEqual([]);
    const writer = await store.create('runs/r1/out6');
    writer.append(new Uint8Array([9]));
    expect(await writer.seal()).toBe(1);
  });
});

// WA2-SC-05: Chromium also removed them by itself (a private context whose tab was hidden for half an hour), and the
// message blamed a cleared site data that the user had not cleared. The e2e tests read this constant.
describe('RESULTS_REMOVED_MESSAGE', () => {
  it('says that the browser removed the stored results, and names no cause', () => {
    expect(RESULTS_REMOVED_MESSAGE).toBe('The browser removed the stored results of this run. Run it again.');
  });
});
