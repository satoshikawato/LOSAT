// What the saved-file stores share (example-store.ts for the Example sets, counts-store.ts for the
// WindowMasker counts): files written with synchronous access handles (which exist only in
// dedicated workers), a tolerant removal, and the reading of a small JSON file. A saved file is a
// data file and then a `meta.json`; the stores decide what the meta file says.
import { StorageFullError } from './block-store';

/** Bytes copied from a Blob at once. */
export const WRITE_CHUNK_BYTES = 8 * 1024 * 1024;

export interface SyncAccessHandle {
  write(buffer: BufferSource, options?: { at?: number }): number;
  truncate(size: number): void;
  flush(): void;
  close(): void;
}

type SyncFileHandle = FileSystemFileHandle & { createSyncAccessHandle?: () => Promise<SyncAccessHandle> };

/** Opens the file (made if it is not there) for synchronous writing; the caller closes the handle. */
export async function openSyncHandle(directory: FileSystemDirectoryHandle, name: string): Promise<SyncAccessHandle> {
  const file = (await directory.getFileHandle(name, { create: true })) as SyncFileHandle;
  if (typeof file.createSyncAccessHandle !== 'function') {
    throw new Error('this browser cannot write OPFS files from a worker (no createSyncAccessHandle)');
  }
  return file.createSyncAccessHandle();
}

/** Writes `bytes` at `at`; a short write is a full storage. */
export function writeAt(handle: SyncAccessHandle, bytes: Uint8Array<ArrayBuffer>, at: number): void {
  if (handle.write(bytes, { at }) !== bytes.length) throw new StorageFullError();
}

/** Replaces the file's content with `data`. */
export async function writeFile(directory: FileSystemDirectoryHandle, name: string, data: Blob): Promise<void> {
  const handle = await openSyncHandle(directory, name);
  try {
    handle.truncate(0);
    for (let at = 0; at < data.size; at += WRITE_CHUNK_BYTES) {
      writeAt(handle, new Uint8Array(await data.slice(at, at + WRITE_CHUNK_BYTES).arrayBuffer()), at);
    }
    handle.flush();
  } finally {
    handle.close();
  }
}

/** The text of a small file, or undefined when it is not there. */
export async function readText(directory: FileSystemDirectoryHandle, name: string): Promise<string | undefined> {
  try {
    return await (await (await directory.getFileHandle(name)).getFile()).text();
  } catch (error) {
    ignoreNotFound(error);
    return undefined;
  }
}

/**
 * Removes an entry and what it holds. A synchronous access handle releases its file shortly
 * after `close` returns, so a removal right after a close can see the file still in use; it
 * is tried again for up to a second.
 */
export async function removeEntry(directory: FileSystemDirectoryHandle, name: string): Promise<void> {
  for (let attempt = 0; ; attempt++) {
    try {
      await directory.removeEntry(name, { recursive: true });
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

export function ignoreNotFound(error: unknown): void {
  if ((error as { name?: unknown }).name !== 'NotFoundError') throw error;
}
