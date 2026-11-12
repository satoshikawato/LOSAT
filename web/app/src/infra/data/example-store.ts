// OPFS implementation of the ExampleStore port (line WF-2, docs/web/fetch_design.md §9): saved
// copies of Example sets, kept across sessions under examples/<example-id>/, beside the working
// sessions' tmp/ (session.ts), which start-up cleanups never touch. A copy is the file `data`
// and then `meta.json`, written last, so a directory without meta.json is an unfinished copy,
// removed by the next save. Saves and removals of a set hold the Web Lock
// `losat-web:example:<example-id>`, so two tabs never write one copy. It writes with
// synchronous access handles, which exist only in dedicated workers: it runs in the Data worker.
import type { ExampleId } from '../../domain/examples';
import type { ExampleCopy, ExampleStore, SavedExample } from '../../ports/example';
import { isStorageFull, StorageFullError } from './block-store';
import { ignoreNotFound, readText, removeEntry, writeFile } from './opfs-files';
import type { SessionLocks } from './session';

/** The OPFS directory of the saved copies. Shared by every release. */
export const EXAMPLES_DIRECTORY = 'examples';
/** Lock name prefix; shared by every release. */
export const EXAMPLE_LOCK_PREFIX = 'losat-web:example:';

const DATA_NAME = 'data';
const META_NAME = 'meta.json';

export class OpfsExampleStore implements ExampleStore {
  constructor(
    private readonly root: FileSystemDirectoryHandle,
    private readonly locks: SessionLocks,
  ) {}

  async readExample(id: ExampleId, url: string): Promise<SavedExample | undefined> {
    let directory: FileSystemDirectoryHandle;
    try {
      directory = await this.root.getDirectoryHandle(id);
    } catch (error) {
      ignoreNotFound(error);
      return undefined;
    }
    const copy = await readCopy(directory, id);
    if (copy === undefined || copy.url !== url) return undefined;
    let file: File;
    try {
      file = await (await directory.getFileHandle(DATA_NAME)).getFile();
    } catch (error) {
      ignoreNotFound(error);
      return undefined;
    }
    return file.size === copy.bytes ? { copy, file } : undefined;
  }

  async saveExample(copy: Omit<ExampleCopy, 'bytes' | 'savedAt'>, data: Blob): Promise<SavedExample | undefined> {
    return this.locked(copy.id, async () => {
      const saved = await this.readExample(copy.id, copy.url);
      if (saved !== undefined) return saved;
      await removeEntry(this.root, copy.id);
      try {
        const directory = await this.root.getDirectoryHandle(copy.id, { create: true });
        await writeFile(directory, DATA_NAME, data);
        const complete: ExampleCopy = { ...copy, bytes: data.size, savedAt: Date.now() };
        await writeFile(directory, META_NAME, new Blob([JSON.stringify(complete)]));
      } catch (error) {
        await removeEntry(this.root, copy.id).catch(() => undefined);
        throw isStorageFull(error) ? new StorageFullError('The browser refused to keep the example: its storage is full.') : error;
      }
      return this.readExample(copy.id, copy.url);
    });
  }

  async deleteExample(id: ExampleId): Promise<void> {
    await this.locked(id, () => removeEntry(this.root, id));
  }

  private locked<T>(id: ExampleId, task: () => Promise<T>): Promise<T> {
    return new Promise<T>((resolve, reject) => {
      this.locks
        .request(EXAMPLE_LOCK_PREFIX + id, {}, async () => {
          try {
            resolve(await task());
          } catch (error) {
            reject(error as Error);
          }
        })
        .catch(reject);
    });
  }
}

/** The saved copies in this origin's OPFS (examples/). Rejects when this browser has no usable OPFS. */
export async function opfsExampleStore(locks: SessionLocks): Promise<OpfsExampleStore> {
  const storage = (globalThis.navigator as Navigator | undefined)?.storage;
  if (typeof storage?.getDirectory !== 'function') throw new Error('this browser has no Origin Private File System');
  const root = await (await storage.getDirectory()).getDirectoryHandle(EXAMPLES_DIRECTORY, { create: true });
  return new OpfsExampleStore(root, locks);
}

async function readCopy(directory: FileSystemDirectoryHandle, id: ExampleId): Promise<ExampleCopy | undefined> {
  const text = await readText(directory, META_NAME);
  if (text === undefined) return undefined;
  try {
    const copy = JSON.parse(text) as ExampleCopy;
    return copy.id === id && typeof copy.url === 'string' && Number.isSafeInteger(copy.bytes) ? copy : undefined;
  } catch {
    return undefined;
  }
}
