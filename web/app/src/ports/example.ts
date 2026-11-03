// Example set ports (line WF-2, docs/web/fetch_design.md §9). The download runs on the page,
// only when the user asks ("Example"), and its request carries only the file's URL
// (web/AGENTS.md rule 4): src/infra/browser/example-downloader.ts. The copies kept across
// sessions live in the Data worker (OPFS examples/): src/infra/data/example-store.ts. Tests use
// src/infra/fake/fake-examples.ts.
import type { ExampleId } from '../domain/examples';

export interface DownloadRequest {
  readonly url: string;
  /** The download stops ('too-large') once the decompressed bytes pass this. */
  readonly maxBytes: number;
}

/** Bytes downloaded so far, and the whole download when the site says it (Content-Length). */
export interface DownloadProgress {
  readonly bytes: number;
  readonly total?: number;
}

export type DownloadResult =
  /** `data` is the file decompressed; `version` names the release when the site says it. */
  | { readonly ok: true; readonly data: Blob; readonly downloadedBytes: number; readonly version?: string }
  /** No answer: offline, blocked, or timed out. */
  | { readonly ok: false; readonly reason: 'unreachable'; readonly detail: string }
  /** An HTTP status other than 200 (such as 404, or a 500 that may pass on another try). */
  | { readonly ok: false; readonly reason: 'refused'; readonly status: number }
  /** The gzip data could not be decompressed (cut short or damaged). */
  | { readonly ok: false; readonly reason: 'corrupt'; readonly detail: string }
  | { readonly ok: false; readonly reason: 'too-large' };

export interface FileDownloader {
  /** Rejects only when `signal` aborts. */
  download(request: DownloadRequest, onProgress: (progress: DownloadProgress) => void, signal: AbortSignal): Promise<DownloadResult>;
}

/** What a saved copy of an Example set records about its download. */
export interface ExampleCopy {
  readonly id: ExampleId;
  /** The URL that it was downloaded from; a copy of another URL is not the set's. */
  readonly url: string;
  /** The release that the site named (a date from Last-Modified, or UniProt's release), if any. */
  readonly version?: string;
  readonly downloadedBytes: number;
  /** Bytes kept: the file decompressed. */
  readonly bytes: number;
  /** When it was saved, in milliseconds since the epoch. */
  readonly savedAt: number;
}

export interface SavedExample {
  readonly copy: ExampleCopy;
  /** The kept file; it reads from the browser's storage, not from memory. */
  readonly file: File;
}

/**
 * Copies of Example sets kept across sessions and tabs, so that a set is downloaded once.
 * A copy is complete or absent; two tabs never write one copy at the same time. Without
 * storage that lasts (no OPFS or no Web Locks), nothing is kept.
 */
export interface ExampleStore {
  /** The saved copy of the set, if one was downloaded from `url`. */
  readExample(id: ExampleId, url: string): Promise<SavedExample | undefined>;
  /**
   * Keeps a downloaded set and returns the kept copy. A copy of the same URL that another tab
   * saved meanwhile is kept instead and returned. Resolves undefined when this browser keeps
   * nothing across sessions; rejects when the storage refused the data (it is then full).
   */
  saveExample(copy: Omit<ExampleCopy, 'bytes' | 'savedAt'>, data: Blob): Promise<SavedExample | undefined>;
  /** Removes the saved copy of the set, if there is one. */
  deleteExample(id: ExampleId): Promise<void>;
}

/**
 * Deletions of saved copies between the tabs of the site: a tab that deletes a copy says so, and
 * the others drop the sources made from it (their File can no longer be read).
 */
export interface ExampleDeletions {
  announce(id: ExampleId): void;
  /** Calls `listener` with each id that another tab deleted; returns the function that stops it. */
  subscribe(listener: (id: string) => void): () => void;
}
