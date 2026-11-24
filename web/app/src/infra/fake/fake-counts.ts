// WindowMasker counts doubles for the tests (ports/counts.ts): a store that keeps the files in
// memory (or keeps nothing, as a browser without OPFS does) with the real store's rules (one
// writing at a time, a continued part must be the saved one), a downloader that answers from a
// table by taxonomy ID and keeps the requests, and the deletions between tabs.
import { COUNTS_FILE } from '../../domain/counts';
import { dateOf } from '../../domain/download-sentences';
import type {
  CountsChanges,
  CountsCopy,
  CountsDownloader,
  CountsRequest,
  CountsResult,
  CountsStart,
  CountsStore,
  SavedCounts,
  SavedCountsEntry,
  SavedCountsOrigin,
} from '../../ports/counts';
import { COUNTS_BUSY, COUNTS_CHANGED, COUNTS_EXISTS } from '../../ports/counts';
import type { DownloadProgress } from '../../ports/example';

interface Entry {
  copy: CountsCopy;
  bytes: Uint8Array<ArrayBuffer>;
}

function named(name: string, message: string): Error {
  const error = new Error(message);
  error.name = name;
  return error;
}

export class MemoryCountsStore implements CountsStore {
  readonly entries = new Map<number, Entry>();
  private readonly writings = new Map<string, { readonly taxid: number; readonly start: CountsStart }>();
  private tokens = 0;
  /** The error that the next appends reject with (a full storage). */
  appendError: Error | undefined;

  /** `persistent: false` keeps nothing, as a browser without OPFS or Web Locks. */
  constructor(private readonly options: { readonly persistent: boolean; readonly now?: () => number } = { persistent: true }) {}

  async readCounts(taxid: number): Promise<SavedCounts | undefined> {
    const entry = this.entries.get(taxid);
    if (entry === undefined) return undefined;
    if (entry.copy.made !== undefined && !entry.copy.complete) {
      if (this.isLocked(taxid)) return undefined;
      this.entries.delete(taxid);
      return undefined;
    }
    const origin: SavedCountsOrigin = entry.copy.made === undefined ? 'downloaded' : 'made';
    return entry.copy.complete ? { copy: entry.copy, origin, file: new File([entry.bytes], COUNTS_FILE) } : { copy: entry.copy, origin };
  }

  async listCounts(): Promise<readonly SavedCountsEntry[]> {
    const taxids = [...this.entries.keys()].sort((a, b) => a - b);
    const list: SavedCountsEntry[] = [];
    for (const taxid of taxids) {
      const saved = await this.readCounts(taxid);
      if (saved?.file === undefined) continue;
      const date =
        saved.copy.made !== undefined
          ? new Date(saved.copy.made.madeAt).toISOString().slice(0, 10)
          : (dateOf(saved.copy.lastModified) ?? new Date(saved.copy.savedAt).toISOString().slice(0, 10));
      list.push({ taxid, name: COUNTS_FILE, size: saved.copy.size, origin: saved.origin, date, copy: saved.copy });
    }
    return list;
  }

  async beginCounts(taxid: number, start: CountsStart): Promise<string | undefined> {
    if (!this.options.persistent) return undefined;
    if (this.isLocked(taxid)) throw named(COUNTS_BUSY, 'busy');
    const entry = this.entries.get(taxid);
    if (start.made !== undefined) {
      if (start.from !== 0) throw named(COUNTS_CHANGED, 'changed');
      if (start.replace === false) {
        const copy = entry?.copy;
        if (copy !== undefined && (copy.complete || (copy.made === undefined && copy.savedBytes > 0))) {
          throw named(COUNTS_EXISTS, 'exists');
        }
      }
      this.entries.set(taxid, { copy: this.copyOf(taxid, start, 0, false), bytes: new Uint8Array(0) });
    } else if (start.from > 0) {
      const copy = entry?.copy;
      if (
        copy === undefined ||
        copy.complete ||
        copy.made !== undefined ||
        copy.url !== start.url ||
        copy.lastModified !== start.lastModified ||
        copy.size !== start.size ||
        copy.savedBytes < start.from
      ) {
        throw named(COUNTS_CHANGED, 'changed');
      }
      // A copy is continued at its byte or before it (its last byte again, when all came but it was not finished).
      this.entries.set(taxid, { copy: { ...copy, savedBytes: start.from }, bytes: entry!.bytes.slice(0, start.from) });
    } else {
      this.entries.set(taxid, { copy: this.copyOf(taxid, start, 0, false), bytes: new Uint8Array(0) });
    }
    const token = `w${++this.tokens}`;
    this.writings.set(token, { taxid, start });
    return token;
  }

  async appendCounts(token: string, bytes: Uint8Array<ArrayBuffer>): Promise<void> {
    if (this.appendError !== undefined) throw this.appendError;
    const writing = this.writings.get(token)!;
    const entry = this.entries.get(writing.taxid)!;
    const next = new Uint8Array(entry.bytes.length + bytes.length);
    next.set(entry.bytes);
    next.set(bytes, entry.bytes.length);
    this.entries.set(writing.taxid, { copy: this.copyOf(writing.taxid, writing.start, next.length, false), bytes: next });
  }

  async finishCounts(token: string): Promise<SavedCounts> {
    const writing = this.writings.get(token)!;
    this.writings.delete(token);
    const entry = this.entries.get(writing.taxid)!;
    if (entry.bytes.length !== writing.start.size) throw new Error('the file is not whole');
    this.entries.set(writing.taxid, { copy: this.copyOf(writing.taxid, writing.start, entry.bytes.length, true), bytes: entry.bytes });
    return (await this.readCounts(writing.taxid))!;
  }

  async closeCounts(token: string): Promise<void> {
    this.writings.delete(token);
  }

  async deleteCounts(taxid: number): Promise<void> {
    this.entries.delete(taxid);
  }

  private isLocked(taxid: number): boolean {
    return [...this.writings.values()].some((writing) => writing.taxid === taxid);
  }

  private copyOf(taxid: number, start: CountsStart, savedBytes: number, complete: boolean): CountsCopy {
    const savedAt = this.options.now?.() ?? 0;
    if (start.made !== undefined) {
      return { taxid, url: '', lastModified: '', size: start.size, savedBytes, complete, savedAt, made: start.made };
    }
    return { taxid, url: start.url, lastModified: start.lastModified, size: start.size, savedBytes, complete, savedAt };
  }
}

/** The file that a table answers with, and what it sends as Last-Modified. */
export interface FakeCountsFile {
  readonly bytes: Uint8Array<ArrayBuffer>;
  readonly lastModified: string;
}

/**
 * A downloader that serves the files of a table by taxonomy ID, as the site does: the whole file,
 * or the rest of it for a `resume`. `cutAfter` ends the next download early (a cut connection).
 */
export class FakeCountsDownloader implements CountsDownloader {
  readonly requests: CountsRequest[] = [];
  /** The result that the next download answers with, instead of the file. */
  next: CountsResult | undefined;

  constructor(
    private readonly store: CountsStore,
    private readonly files: Readonly<Record<number, FakeCountsFile>> = {},
  ) {}

  async download(request: CountsRequest, onProgress: (progress: DownloadProgress) => void, signal: AbortSignal): Promise<CountsResult> {
    signal.throwIfAborted();
    this.requests.push(request);
    if (this.next !== undefined) {
      const result = this.next;
      this.next = undefined;
      return result;
    }
    const file = this.files[request.taxid];
    if (file === undefined) return { ok: false, reason: 'refused', status: 404 };
    const resumed = request.resume !== undefined && request.resume.lastModified === file.lastModified;
    const from = resumed ? request.resume!.savedBytes : 0;
    const start = { url: request.url, lastModified: file.lastModified, size: file.bytes.length, from };
    const token = await this.store.beginCounts(request.taxid, start);
    onProgress({ bytes: from, total: file.bytes.length });
    if (token === undefined) {
      onProgress({ bytes: file.bytes.length, total: file.bytes.length });
      const copy: CountsCopy = { taxid: request.taxid, url: request.url, lastModified: file.lastModified, size: file.bytes.length, savedBytes: file.bytes.length, complete: true, savedAt: 0 };
      return { ok: true, copy, file: new File([file.bytes], 'wmasker.obinary'), kept: false };
    }
    await this.store.appendCounts(token, file.bytes.slice(from));
    onProgress({ bytes: file.bytes.length, total: file.bytes.length });
    const saved = await this.store.finishCounts(token);
    return { ok: true, copy: saved.copy, file: saved.file!, kept: true };
  }
}

export class FakeCountsChangeChannel {
  private readonly listeners = new Set<{ readonly tab: number; readonly listener: (taxid: number) => void }>();
  private tabs = 0;

  tab(): CountsChanges {
    const tab = ++this.tabs;
    return {
      announce: (taxid) => {
        for (const entry of this.listeners) if (entry.tab !== tab) entry.listener(taxid);
      },
      subscribe: (listener) => {
        const entry = { tab, listener };
        this.listeners.add(entry);
        return () => this.listeners.delete(entry);
      },
    };
  }
}
