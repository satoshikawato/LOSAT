// OPFS implementation of the BlockStore contract (plan §5.6, design §7.1). It writes with
// synchronous access handles, which exist only in dedicated workers, so it runs only in
// the Data worker. Blocks of one working session live under tmp/<session-token>/; their
// names are internal tokens, never input names.
import {
  checkRange,
  isStorageFull,
  pathNames,
  StorageFullError,
  type BlockStore,
  type BlockWriter,
} from './block-store';
import { TMP_DIRECTORY, type OpfsAccess, type SessionSpace } from './session';

/** The part of FileSystemSyncAccessHandle that the store uses (worker-only API). */
interface SyncAccessHandle {
  write(buffer: BufferSource, options?: { at?: number }): number;
  flush(): void;
  close(): void;
}

type SyncFileHandle = FileSystemFileHandle & { createSyncAccessHandle?: () => Promise<SyncAccessHandle> };

interface OpfsBlock {
  readonly file: FileSystemFileHandle;
  readonly directory: FileSystemDirectoryHandle;
  readonly name: string;
  handle: SyncAccessHandle | undefined;
  length: number;
  state: 'open' | 'sealed' | 'full';
}

export class OpfsBlockStore implements BlockStore {
  readonly backend = 'opfs';
  private readonly blocks = new Map<string, OpfsBlock>();
  private readonly removals = new Map<string, Promise<void>>();
  private used = 0;

  constructor(private readonly root: FileSystemDirectoryHandle) {}

  async create(path: string): Promise<BlockWriter> {
    const names = pathNames(path);
    const name = names[names.length - 1]!;
    if (this.blocks.has(path)) throw new Error(`block ${path} already exists`);
    let directory: FileSystemDirectoryHandle;
    let file: FileSystemFileHandle;
    let handle: SyncAccessHandle;
    try {
      directory = await this.directory(names.slice(0, -1), true);
      if (await hasEntry(directory, name)) throw new Error(`block ${path} already exists`);
      file = await directory.getFileHandle(name, { create: true });
      handle = await openSyncAccessHandle(file);
    } catch (error) {
      throw isStorageFull(error) ? new StorageFullError() : error;
    }
    if (this.blocks.has(path)) {
      handle.close();
      throw new Error(`block ${path} already exists`);
    }
    const block: OpfsBlock = { file, directory, name, handle, length: 0, state: 'open' };
    this.blocks.set(path, block);
    const live = () => this.blocks.get(path) === block;
    return {
      append: (bytes) => {
        if (!live() || block.state !== 'open' || block.handle === undefined) {
          throw new Error(`block ${path} is not open`);
        }
        let written: number;
        try {
          written = block.handle.write(bytes as BufferSource, { at: block.length });
        } catch (error) {
          if (!isStorageFull(error)) throw error;
          block.state = 'full';
          throw new StorageFullError();
        }
        block.length += written;
        this.used += written;
        if (written !== bytes.length) {
          block.state = 'full';
          throw new StorageFullError();
        }
      },
      seal: async () => {
        if (live() && block.state === 'full') throw new StorageFullError();
        if (!live() || block.state !== 'open' || block.handle === undefined) {
          throw new Error(`block ${path} is not open`);
        }
        try {
          block.handle.flush();
        } catch (error) {
          throw isStorageFull(error) ? new StorageFullError() : error;
        }
        block.handle.close();
        block.handle = undefined;
        block.state = 'sealed';
        return block.length;
      },
      discard: async () => {
        if (!live()) return;
        this.forget(path, block);
        await removeEntry(block.directory, block.name, false);
      },
    };
  }

  async read(path: string, offset: number, length: number): Promise<Uint8Array> {
    const block = this.sealed(path);
    checkRange(path, offset, length, block.length);
    const file = await block.file.getFile();
    return new Uint8Array(await file.slice(offset, offset + length).arrayBuffer());
  }

  async size(path: string): Promise<number> {
    return this.sealed(path).length;
  }

  async removeAll(prefix: string): Promise<void> {
    const names = pathNames(prefix);
    for (const [path, block] of [...this.blocks]) {
      if (path.startsWith(`${prefix}/`)) this.forget(path, block);
    }
    // Removals of one directory run one after another, never at the same time.
    const removal = (this.removals.get(prefix) ?? Promise.resolve()).then(async () => {
      let parent: FileSystemDirectoryHandle;
      try {
        parent = await this.directory(names.slice(0, -1), false);
      } catch (error) {
        ignoreNotFound(error);
        return;
      }
      await removeEntry(parent, names[names.length - 1]!, true);
    });
    const settled = removal.catch(() => undefined);
    this.removals.set(prefix, settled);
    try {
      await removal;
    } finally {
      if (this.removals.get(prefix) === settled) this.removals.delete(prefix);
    }
  }

  usage(): number {
    return this.used;
  }

  private sealed(path: string): OpfsBlock {
    const block = this.blocks.get(path);
    if (block === undefined) throw new Error(`block ${path} does not exist`);
    if (block.state !== 'sealed') throw new Error(`block ${path} is not sealed`);
    return block;
  }

  private forget(path: string, block: OpfsBlock): void {
    this.blocks.delete(path);
    this.used -= block.length;
    try {
      block.handle?.close();
    } catch {
      // Already closed.
    }
    block.handle = undefined;
  }

  private async directory(names: readonly string[], create: boolean): Promise<FileSystemDirectoryHandle> {
    let directory = this.root;
    for (const name of names) directory = await directory.getDirectoryHandle(name, { create });
    return directory;
  }
}

/** OPFS for the working sessions (session.ts): tmp/ and the directory of each session. */
export const opfsAccess: OpfsAccess = {
  async space() {
    return opfsSessionSpace(await tmpDirectory());
  },
  /**
   * Opens tmp/<token>/, checked by writing, sealing, reading and removing a small block.
   * Call it only while the session's lock is held. Rejects, after removing what it
   * created, when this browser cannot use OPFS here.
   */
  async open(token) {
    const tmp = await tmpDirectory();
    const store = new OpfsBlockStore(await tmp.getDirectoryHandle(token, { create: true }));
    try {
      await probe(store);
    } catch (error) {
      await tmp.removeEntry(token, { recursive: true }).catch(() => undefined);
      throw error;
    }
    return store;
  },
};

async function tmpDirectory(): Promise<FileSystemDirectoryHandle> {
  const storage = (globalThis.navigator as Navigator | undefined)?.storage;
  if (typeof storage?.getDirectory !== 'function') {
    throw new Error('this browser has no Origin Private File System');
  }
  return (await storage.getDirectory()).getDirectoryHandle(TMP_DIRECTORY, { create: true });
}

/** The session directories under tmp/ (session.ts). */
function opfsSessionSpace(tmp: FileSystemDirectoryHandle): SessionSpace {
  return {
    async list() {
      const names: string[] = [];
      const entries = (tmp as unknown as { entries(): AsyncIterable<[string, FileSystemHandle]> }).entries();
      for await (const [name, handle] of entries) {
        if (handle.kind === 'directory') names.push(name);
      }
      return names;
    },
    async remove(token) {
      await tmp.removeEntry(token, { recursive: true });
    },
  };
}

async function probe(store: OpfsBlockStore): Promise<void> {
  const bytes = new Uint8Array([0x4c, 0x4f, 0x53, 0x41, 0x54]);
  const writer = await store.create('probe/check');
  writer.append(bytes);
  await writer.seal();
  const back = await store.read('probe/check', 0, bytes.length);
  if (back.some((byte, i) => byte !== bytes[i])) throw new Error('OPFS returned other bytes than were written');
  await store.removeAll('probe');
}

async function openSyncAccessHandle(file: FileSystemFileHandle): Promise<SyncAccessHandle> {
  const open = (file as SyncFileHandle).createSyncAccessHandle;
  if (typeof open !== 'function') {
    throw new Error('this browser cannot write OPFS files from a worker (no createSyncAccessHandle)');
  }
  return open.call(file);
}

async function hasEntry(directory: FileSystemDirectoryHandle, name: string): Promise<boolean> {
  try {
    await directory.getFileHandle(name);
    return true;
  } catch (error) {
    if ((error as { name?: unknown }).name === 'TypeMismatchError') return true;
    ignoreNotFound(error);
    return false;
  }
}

/**
 * Removes an entry. A synchronous access handle releases its file shortly after `close`
 * returns, so a removal right after a close can see the file still in use; it is tried
 * again for up to a second.
 */
async function removeEntry(directory: FileSystemDirectoryHandle, name: string, recursive: boolean): Promise<void> {
  for (let attempt = 0; ; attempt++) {
    try {
      await directory.removeEntry(name, { recursive });
      return;
    } catch (error) {
      const busy = (error as { name?: unknown }).name === 'NoModificationAllowedError';
      if (!busy || attempt >= 40) {
        ignoreNotFound(error);
        return;
      }
      await new Promise((resolve) => setTimeout(resolve, 25));
    }
  }
}

function ignoreNotFound(error: unknown): void {
  if ((error as { name?: unknown }).name !== 'NotFoundError') throw error;
}
