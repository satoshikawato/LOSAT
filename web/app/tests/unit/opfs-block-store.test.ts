// The OPFS block store when the browser removes its directory under it (the site's data was
// cleared while the tab stayed open): a fake directory tree that fails as Chromium and Firefox do.
import { describe, expect, it } from 'vitest';
import { RESULTS_REMOVED_MESSAGE, StorageRemovedError, STORAGE_REMOVED_MESSAGE } from '../../src/infra/data/block-store';
import { OpfsBlockStore } from '../../src/infra/data/opfs-block-store';

class FakeFile {
  readonly kind = 'file';
  bytes = new Uint8Array(0);
  constructor(private readonly parent: FakeDirectory) {}
  async createSyncAccessHandle() {
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
    function setup(reopen = true) {
      const top = new FakeDirectory(undefined, failure);
      const session = () => top.getDirectoryHandle('session', { create: true });
      let first!: FakeDirectory;
      return session().then((made) => {
        first = made;
        const store = new OpfsBlockStore(asHandle(first), reopen ? async () => asHandle(await session()) : undefined);
        // What clearing the site's data does: the directory goes, the page and its worker stay.
        const wipe = () => {
          first.removed = true;
          top.entries.delete('session');
        };
        return { store, wipe };
      });
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
      const { store, wipe } = await setup(false);
      wipe();
      const error = await store.create('runs/r1/out6').catch((e: unknown) => e);
      expect(error).toBeInstanceOf(StorageRemovedError);
      expect((error as Error).message).toBe(STORAGE_REMOVED_MESSAGE);
      expect((error as Error).message).toMatch(/Reload the page\.$/);
    });
  });
}
