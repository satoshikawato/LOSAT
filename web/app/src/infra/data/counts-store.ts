// OPFS implementation of the CountsStore port (line WW2, docs/web/fetch_design.md §10): the
// WindowMasker counts of the taxonomy IDs that were downloaded, kept across sessions under
// windowmasker/<taxid>/, beside examples/ and the working sessions' tmp/ (session.ts), which
// start-up cleanups never touch. The path is the one the search reads the file under
// (windowmasker/<taxid>/wmasker.obinary, docs/web/abi_v2.md §4). A directory holds the file
// `wmasker.obinary` and `meta.json`; the body is appended to the file part by part (the file of a
// human is 724 MB: it is never held in memory), and the meta file is rewritten after each part
// with the bytes that are saved, which is where "Continue" starts. A writing holds the Web Lock
// `losat-web:windowmasker:<taxid>` from `beginCounts` to its end (or the tab's), so two tabs never
// write one file. It writes with synchronous access handles, which exist only in dedicated
// workers: it runs in the Data worker.
import { COUNTS_BUSY, COUNTS_CHANGED, type CountsCopy, type CountsStart, type CountsStore, type SavedCounts } from '../../ports/counts';
import { isStorageFull, StorageFullError } from './block-store';
import { ignoreNotFound, openSyncHandle, readText, removeEntry, writeAt, writeFile, type SyncAccessHandle } from './opfs-files';
import type { SessionLocks } from './session';

/** The OPFS directory of the saved counts. Shared by every release. */
export const COUNTS_DIRECTORY = 'windowmasker';
/** Lock name prefix; shared by every release. */
export const COUNTS_LOCK_PREFIX = 'losat-web:windowmasker:';

const DATA_NAME = 'wmasker.obinary';
const META_NAME = 'meta.json';

/** A file being written: the handle, the next byte, and the lock to let go. */
interface Writing {
  readonly taxid: number;
  readonly directory: FileSystemDirectoryHandle;
  readonly handle: SyncAccessHandle;
  readonly start: CountsStart;
  readonly release: () => void;
  offset: number;
}

export class OpfsCountsStore implements CountsStore {
  private readonly writings = new Map<string, Writing>();

  constructor(
    private readonly root: FileSystemDirectoryHandle,
    private readonly locks: SessionLocks,
    private readonly newToken: () => string = () => crypto.randomUUID(),
  ) {}

  async readCounts(taxid: number, url: string): Promise<SavedCounts | undefined> {
    let directory: FileSystemDirectoryHandle;
    try {
      directory = await this.root.getDirectoryHandle(String(taxid));
    } catch (error) {
      ignoreNotFound(error);
      return undefined;
    }
    const copy = await readCopy(directory, taxid);
    if (copy === undefined || copy.url !== url) return undefined;
    let file: File;
    try {
      file = await (await directory.getFileHandle(DATA_NAME)).getFile();
    } catch (error) {
      ignoreNotFound(error);
      return undefined;
    }
    // The file may hold more than the meta file says (a part written, the meta not yet): the
    // meta file's bytes count. It never holds less.
    if (file.size < copy.savedBytes) return undefined;
    return copy.complete ? (file.size === copy.size ? { copy, file } : undefined) : { copy };
  }

  async beginCounts(taxid: number, start: CountsStart): Promise<string> {
    const release = await this.acquire(taxid);
    let handle: SyncAccessHandle | undefined;
    try {
      const directory = await this.root.getDirectoryHandle(String(taxid), { create: true });
      if (start.from > 0) {
        const copy = await readCopy(directory, taxid);
        if (
          copy === undefined ||
          copy.complete ||
          copy.url !== start.url ||
          copy.lastModified !== start.lastModified ||
          copy.size !== start.size ||
          copy.savedBytes < start.from
        ) {
          throw named(COUNTS_CHANGED, 'the saved part is not the one to continue');
        }
      } else {
        // Meta first: a copy that is begun again is never read as the whole file of before.
        await writeMeta(directory, metaOf(taxid, start, 0, false));
      }
      handle = await openSyncHandle(directory, DATA_NAME);
      handle.truncate(start.from);
      const token = this.newToken();
      this.writings.set(token, { taxid, directory, handle, start, release, offset: start.from });
      return token;
    } catch (error) {
      handle?.close();
      release();
      throw refusal(error);
    }
  }

  async appendCounts(token: string, bytes: Uint8Array<ArrayBuffer>): Promise<void> {
    const writing = this.writing(token);
    try {
      writeAt(writing.handle, bytes, writing.offset);
      writing.handle.flush();
      writing.offset += bytes.length;
      await writeMeta(writing.directory, metaOf(writing.taxid, writing.start, writing.offset, false));
    } catch (error) {
      throw refusal(error);
    }
  }

  async finishCounts(token: string): Promise<SavedCounts> {
    const writing = this.writing(token);
    try {
      if (writing.offset !== writing.start.size) throw new Error('the file is not whole');
      await writeMeta(writing.directory, metaOf(writing.taxid, writing.start, writing.offset, true));
    } catch (error) {
      throw refusal(error);
    }
    this.end(token, writing);
    const saved = await this.readCounts(writing.taxid, writing.start.url);
    if (saved?.file === undefined) throw new Error('the saved file cannot be read');
    return saved;
  }

  async closeCounts(token: string): Promise<void> {
    const writing = this.writings.get(token);
    if (writing !== undefined) this.end(token, writing);
  }

  async deleteCounts(taxid: number): Promise<void> {
    const release = await this.acquire(taxid, false);
    try {
      await removeEntry(this.root, String(taxid));
    } finally {
      release();
    }
  }

  private writing(token: string): Writing {
    const writing = this.writings.get(token);
    if (writing === undefined) throw new Error('the counts are no longer being saved');
    return writing;
  }

  /** Closes the file and lets the lock go. */
  private end(token: string, writing: Writing): void {
    this.writings.delete(token);
    try {
      writing.handle.close();
    } finally {
      writing.release();
    }
  }

  /**
   * Takes the taxonomy ID's lock and keeps it until the returned function is called. With
   * `ifAvailable`, rejects with COUNTS_BUSY when another tab holds it; without, waits for it.
   */
  private acquire(taxid: number, ifAvailable = true): Promise<() => void> {
    return new Promise((resolve, reject) => {
      this.locks
        .request(COUNTS_LOCK_PREFIX + taxid, { ifAvailable }, (lock) => {
          if (lock === null) {
            reject(named(COUNTS_BUSY, 'another tab is downloading the counts'));
            return Promise.resolve();
          }
          return new Promise<void>((held) => resolve(() => held()));
        })
        .catch(reject);
    });
  }
}

/** The saved counts in this origin's OPFS (windowmasker/). Rejects when this browser has no usable OPFS. */
export async function opfsCountsStore(locks: SessionLocks): Promise<OpfsCountsStore> {
  const storage = (globalThis.navigator as Navigator | undefined)?.storage;
  if (typeof storage?.getDirectory !== 'function') throw new Error('this browser has no Origin Private File System');
  const root = await (await storage.getDirectory()).getDirectoryHandle(COUNTS_DIRECTORY, { create: true });
  return new OpfsCountsStore(root, locks);
}

function metaOf(taxid: number, start: CountsStart, savedBytes: number, complete: boolean): CountsCopy {
  return { taxid, url: start.url, lastModified: start.lastModified, size: start.size, savedBytes, complete, savedAt: Date.now() };
}

function writeMeta(directory: FileSystemDirectoryHandle, copy: CountsCopy): Promise<void> {
  return writeFile(directory, META_NAME, new Blob([JSON.stringify(copy)]));
}

async function readCopy(directory: FileSystemDirectoryHandle, taxid: number): Promise<CountsCopy | undefined> {
  const text = await readText(directory, META_NAME);
  if (text === undefined) return undefined;
  try {
    const copy = JSON.parse(text) as CountsCopy;
    const whole = (value: unknown): value is number => typeof value === 'number' && Number.isSafeInteger(value) && value >= 0;
    return copy.taxid === taxid &&
      typeof copy.url === 'string' &&
      typeof copy.lastModified === 'string' &&
      whole(copy.size) &&
      whole(copy.savedBytes) &&
      copy.savedBytes <= copy.size &&
      typeof copy.complete === 'boolean' &&
      (!copy.complete || copy.savedBytes === copy.size)
      ? copy
      : undefined;
  } catch {
    return undefined;
  }
}

function named(name: string, message: string): Error {
  const error = new Error(message);
  error.name = name;
  return error;
}

/** A full storage in the words of the line ("Not kept: ..."); any other error as it is. */
function refusal(error: unknown): unknown {
  return isStorageFull(error) ? new StorageFullError('The browser refused to keep the counts: its storage is full.') : error;
}
