// Offline use and the explicit update (plan §5.9, design §14): the site's service worker,
// the pinned release, the other tabs of the site, and the page's network state. The browser
// implementations are in src/infra/browser/offline.ts; tests use fakes.

/** Files and bytes of a preparation: done of total. */
export interface PrepareProgress {
  readonly files: number;
  readonly filesTotal: number;
  readonly bytes: number;
  readonly bytesTotal: number;
}

export type PrepareResult = { readonly ok: true } | { readonly ok: false; readonly reason: string };

/** The site's service worker, which stores releases for offline use. */
export interface OfflinePort {
  /** Whether the browser has service workers; without them the app works online only. */
  readonly supported: boolean;
  /** Registers the site's worker (or finds it) and waits until it is active; rejects with the reason it is not. */
  connect(): Promise<void>;
  /**
   * Has the worker download every file of `release` and of the site's shell, check each
   * against its listed size and SHA-256, and store them only when all matched.
   */
  prepare(release: string, onProgress: (progress: PrepareProgress) => void): Promise<PrepareResult>;
  /** Has the worker remove the caches of releases that no tab uses, except those of `keep`. */
  cleanup(keep: readonly string[]): Promise<void>;
  /** Tells a waiting (newer) worker to take over; resolves false when there is none. */
  activateWaiting(): Promise<boolean>;
  /**
   * Holds this page's share of `release`'s lock: the worker keeps the release while anyone holds
   * it. Held until the returned function is called, or for the life of the page.
   */
  holdRelease(release: string): () => void;
  /** Whether `release` is stored whole (its cache holds its manifest, which is stored last); false when it cannot be read. */
  stored(release: string): Promise<boolean>;
}

export type PinWrite = { readonly ok: true } | { readonly ok: false; readonly reason: string };

/** The pinned release, `localStorage['losat-web:release']`, which the boot page starts. */
export interface ReleasePin {
  /** The pinned release; undefined when none is pinned, the value is not a release ID, or storage cannot be read. */
  read(): string | undefined;
  write(release: string): PinWrite;
}

/** The other tabs of the site: a switch of the pinned release made in one is told to the others. */
export interface ReleaseChannel {
  announce(release: string): void;
  /** Returns a function that removes the listener. */
  onSwitch(listener: (release: string) => void): () => void;
}

/** What the release manager needs of the page besides its visibility (ports/page.ts). */
export interface ReleasePage {
  /** False when the browser knows it is offline. */
  online(): boolean;
  /** Returns a function that removes the listener. */
  onOnline(listener: () => void): () => void;
  /** Resolves once the page has started and the browser is idle, or after a short wait. */
  whenIdle(): Promise<void>;
  reload(): void;
}
