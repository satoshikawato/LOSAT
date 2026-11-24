// An in-memory Origin Private File System for the unit tests of the stores that write with
// synchronous access handles (counts-store.ts), and the Web Locks they take (SessionLocks).
import type { SessionLocks } from '../../../src/infra/data/session';

class MemoryFile {
  bytes: Uint8Array<ArrayBuffer> = new Uint8Array(0);
}

function notFound(): DOMException {
  return new DOMException('Entry not found', 'NotFoundError');
}

class FileHandle {
  readonly kind = 'file';
  constructor(
    readonly name: string,
    private readonly file: MemoryFile,
  ) {}

  async getFile(): Promise<File> {
    return new File([this.file.bytes.slice()], this.name);
  }

  async createSyncAccessHandle() {
    const file = this.file;
    return {
      write(buffer: Uint8Array, options?: { at?: number }): number {
        const at = options?.at ?? 0;
        const next = new Uint8Array(Math.max(file.bytes.length, at + buffer.length));
        next.set(file.bytes);
        next.set(buffer, at);
        file.bytes = next;
        return buffer.length;
      },
      truncate(size: number): void {
        const next = new Uint8Array(size);
        next.set(file.bytes.subarray(0, size));
        file.bytes = next;
      },
      flush(): void {},
      close(): void {},
    };
  }
}

class DirectoryHandle {
  readonly kind = 'directory';
  readonly children = new Map<string, MemoryFile | DirectoryHandle>();

  async getDirectoryHandle(name: string, options?: { create?: boolean }): Promise<DirectoryHandle> {
    let child = this.children.get(name);
    if (child === undefined) {
      if (options?.create !== true) throw notFound();
      child = new DirectoryHandle();
      this.children.set(name, child);
    }
    if (!(child instanceof DirectoryHandle)) throw new DOMException('not a directory', 'TypeMismatchError');
    return child;
  }

  async getFileHandle(name: string, options?: { create?: boolean }): Promise<FileHandle> {
    let child = this.children.get(name);
    if (child === undefined) {
      if (options?.create !== true) throw notFound();
      child = new MemoryFile();
      this.children.set(name, child);
    }
    if (!(child instanceof MemoryFile)) throw new DOMException('not a file', 'TypeMismatchError');
    return new FileHandle(name, child);
  }

  async removeEntry(name: string): Promise<void> {
    if (!this.children.delete(name)) throw notFound();
  }

  async *entries(): AsyncIterableIterator<[string, FileHandle | DirectoryHandle]> {
    for (const [name, child] of this.children.entries()) {
      yield [name, child instanceof DirectoryHandle ? child : new FileHandle(name, child)];
    }
  }
}

/** The tree (`root`) and the locks of one origin (Web Locks: exclusive and shared modes, `ifAvailable`, `query`). */
export class MemoryOpfs implements SessionLocks {
  readonly root = new DirectoryHandle();
  /** The held locks: the mode of each name and how many hold it (shared locks are held by several). */
  private readonly held = new Map<string, { mode: 'exclusive' | 'shared'; count: number }>();
  private readonly waiting = new Map<string, Array<() => void>>();

  private find(path: readonly string[]): MemoryFile | DirectoryHandle | undefined {
    let at: MemoryFile | DirectoryHandle | undefined = this.root;
    for (const name of path) {
      if (!(at instanceof DirectoryHandle)) return undefined;
      at = at.children.get(name);
    }
    return at;
  }

  has(path: readonly string[]): boolean {
    return this.find(path) !== undefined;
  }

  read(path: readonly string[]): Uint8Array {
    const file = this.find(path);
    if (!(file instanceof MemoryFile)) throw new Error(`no file ${path.join('/')}`);
    return file.bytes;
  }

  write(path: readonly string[], bytes: Uint8Array<ArrayBuffer>): void {
    const file = this.find(path);
    if (!(file instanceof MemoryFile)) throw new Error(`no file ${path.join('/')}`);
    file.bytes = bytes;
  }

  async request(
    name: string,
    options: { readonly ifAvailable?: boolean; readonly mode?: 'exclusive' | 'shared' },
    callback: (lock: object | null) => Promise<void>,
  ): Promise<unknown> {
    const mode = options.mode ?? 'exclusive';
    const grantable = () => {
      const now = this.held.get(name);
      return now === undefined || (mode === 'shared' && now.mode === 'shared' && (this.waiting.get(name)?.length ?? 0) === 0);
    };
    if (!grantable()) {
      if (options.ifAvailable === true) return callback(null);
      while (!grantable()) {
        await new Promise<void>((resolve) => this.waiting.set(name, [...(this.waiting.get(name) ?? []), resolve]));
      }
    }
    const now = this.held.get(name);
    if (now === undefined) this.held.set(name, { mode, count: 1 });
    else now.count += 1;
    try {
      return await callback({});
    } finally {
      const after = this.held.get(name)!;
      after.count -= 1;
      if (after.count === 0) {
        this.held.delete(name);
        this.waiting.get(name)?.shift()?.();
      }
    }
  }

  async query(): Promise<{ readonly held: ReadonlyArray<{ readonly name: string; readonly mode: 'exclusive' | 'shared' }> }> {
    return { held: [...this.held].flatMap(([name, { mode, count }]) => Array.from({ length: count }, () => ({ name, mode }))) };
  }
}
