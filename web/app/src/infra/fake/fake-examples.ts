// Example set doubles for the tests (ports/example.ts): a downloader that answers from a table
// of texts by URL, keeps the requests and can hold its answers until `release`, a store that
// keeps saved copies in memory (or keeps nothing, as a browser without OPFS does), and the
// deletions between tabs (each `tab()` hears the others, as a BroadcastChannel).
import type { ExampleId } from '../../domain/examples';
import type {
  DownloadProgress,
  DownloadRequest,
  DownloadResult,
  ExampleCopy,
  ExampleDeletions,
  ExampleStore,
  FileDownloader,
  SavedExample,
} from '../../ports/example';

export class FakeFileDownloader implements FileDownloader {
  readonly requests: DownloadRequest[] = [];
  private held: Array<() => void> | undefined;

  constructor(private readonly files: Readonly<Record<string, string | DownloadResult>> = {}) {}

  /** Holds the answers from now on, until `release`; progress is reported before the hold. */
  hold(): void {
    this.held ??= [];
  }

  release(): void {
    const held = this.held ?? [];
    this.held = undefined;
    for (const resume of held) resume();
  }

  async download(
    request: DownloadRequest,
    onProgress: (progress: DownloadProgress) => void,
    signal: AbortSignal,
  ): Promise<DownloadResult> {
    this.requests.push(request);
    const file = this.files[request.url];
    if (typeof file === 'string') onProgress({ bytes: Math.floor(file.length / 2), total: file.length });
    const held = this.held;
    await (held === undefined
      ? Promise.resolve()
      : new Promise<void>((resolve) => {
          held.push(resolve);
          signal.addEventListener('abort', () => resolve(), { once: true });
        }));
    signal.throwIfAborted();
    if (file === undefined) return { ok: false, reason: 'refused', status: 404 };
    if (typeof file !== 'string') return file;
    if (file.length > request.maxBytes) return { ok: false, reason: 'too-large' };
    return { ok: true, data: new Blob([file], { type: 'text/plain' }), downloadedBytes: file.length, version: '2026-10-27' };
  }
}

export class MemoryExampleStore implements ExampleStore {
  readonly copies = new Map<ExampleId, SavedExample>();
  /** The error that the next saves reject with (a full storage). */
  saveError: Error | undefined;

  /** `persistent: false` keeps nothing, as a browser without OPFS or Web Locks. */
  constructor(private readonly options: { readonly persistent: boolean; readonly now?: () => number } = { persistent: true }) {}

  async readExample(id: ExampleId, url: string): Promise<SavedExample | undefined> {
    const saved = this.copies.get(id);
    return saved?.copy.url === url ? saved : undefined;
  }

  async saveExample(copy: Omit<ExampleCopy, 'bytes' | 'savedAt'>, data: Blob): Promise<SavedExample | undefined> {
    if (!this.options.persistent) return undefined;
    if (this.saveError !== undefined) throw this.saveError;
    const existing = await this.readExample(copy.id, copy.url);
    if (existing !== undefined) return existing;
    const saved: SavedExample = {
      copy: { ...copy, bytes: data.size, savedAt: this.options.now?.() ?? 0 },
      file: new File([data], 'data', { type: 'text/plain' }),
    };
    this.copies.set(copy.id, saved);
    return saved;
  }

  async deleteExample(id: ExampleId): Promise<void> {
    this.copies.delete(id);
  }
}

export class FakeDeletionChannel {
  private readonly listeners = new Set<{ readonly tab: number; readonly listener: (id: string) => void }>();
  private tabs = 0;

  tab(): ExampleDeletions {
    const tab = ++this.tabs;
    return {
      announce: (id) => {
        for (const entry of this.listeners) if (entry.tab !== tab) entry.listener(id);
      },
      subscribe: (listener) => {
        const entry = { tab, listener };
        this.listeners.add(entry);
        return () => this.listeners.delete(entry);
      },
    };
  }
}
