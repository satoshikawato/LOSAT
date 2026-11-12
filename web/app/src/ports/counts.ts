// WindowMasker counts ports (line WW2, docs/web/fetch_design.md §10). The download runs on the page,
// only when the user asks, and its request carries only the file's URL (web/AGENTS.md rule 4):
// src/infra/browser/counts-downloader.ts. The body goes to the Data worker as it arrives, in parts
// of about 8 MiB, and is written at the end of the saved file (OPFS windowmasker/<taxid>/,
// src/infra/data/counts-store.ts): the file of a human is 724 MB, and is never held in memory.
// A download that stops (Cancel, a cut connection, a closed tab) leaves the bytes saved so far, and
// "Continue" asks for the rest with a Range. Tests use src/infra/fake/fake-counts.ts.
import type { DownloadProgress } from './example';

/** What a saved file of counts records about its download (meta.json). */
export interface CountsCopy {
  readonly taxid: number;
  readonly url: string;
  /** The Last-Modified of the first answer, as the site sent it ('' when it sent none). */
  readonly lastModified: string;
  /** The size of the whole file. */
  readonly size: number;
  /** Bytes kept from the start of the file; Continue asks for the file from this byte. */
  readonly savedBytes: number;
  /** The whole file is kept. */
  readonly complete: boolean;
  /** When the last part was saved, in milliseconds since the epoch. */
  readonly savedAt: number;
}

export interface SavedCounts {
  readonly copy: CountsCopy;
  /** The kept file, read from the browser's storage, not from memory; only when the copy is complete. */
  readonly file?: File;
}

/** What a part is appended to: a new file, or the saved part of the same file. */
export interface CountsStart {
  readonly url: string;
  readonly lastModified: string;
  readonly size: number;
  /** 0 begins the file again; otherwise the bytes already saved (the saved part is checked against it). */
  readonly from: number;
}

/** The name of the error that `beginCounts` rejects with when another tab writes the same file. */
export const COUNTS_BUSY = 'CountsBusyError';
/** The name of the error that `beginCounts` rejects with when the saved part is not the one asked to continue. */
export const COUNTS_CHANGED = 'CountsChangedError';

/**
 * Counts kept across sessions and tabs, so that a file is downloaded once. A file is written by
 * one tab at a time (a Web Lock held from `beginCounts` to `finishCounts` or `closeCounts`, or until
 * the tab ends). Without storage that lasts (no OPFS or no Web Locks), nothing is kept.
 */
export interface CountsStore {
  /** The saved copy of the taxonomy ID (complete or partial), if one was downloaded from `url`. */
  readCounts(taxid: number, url: string): Promise<SavedCounts | undefined>;
  /**
   * Begins writing the file and returns the token of the writing; undefined when this browser
   * keeps nothing across sessions. Rejects with COUNTS_BUSY or COUNTS_CHANGED, or with the
   * storage's refusal.
   */
  beginCounts(taxid: number, start: CountsStart): Promise<string | undefined>;
  /** Appends the next bytes of the file to the writing; rejects when the storage refuses (it is then full). */
  appendCounts(token: string, bytes: Uint8Array<ArrayBuffer>): Promise<void>;
  /** Ends the writing with the whole file kept, and returns the kept copy. */
  finishCounts(token: string): Promise<SavedCounts>;
  /** Ends the writing and keeps what was saved so far. */
  closeCounts(token: string): Promise<void>;
  /** Removes the saved copy of the taxonomy ID, if there is one. */
  deleteCounts(taxid: number): Promise<void>;
}

export interface CountsRequest {
  readonly taxid: number;
  readonly url: string;
  /** The partial copy to continue (a Range from `savedBytes`), when there is one. */
  readonly resume?: CountsCopy;
}

export type CountsResult =
  /** `kept` is false when this browser keeps nothing: `file` is then held by this tab only. */
  | { readonly ok: true; readonly copy: CountsCopy; readonly file: File; readonly kept: boolean }
  /** No answer, or the connection ended early: offline, blocked, timed out. */
  | { readonly ok: false; readonly reason: 'unreachable'; readonly detail: string }
  /** An HTTP status other than 200 or 206 (such as 404, or a 500 that may pass on another try). */
  | { readonly ok: false; readonly reason: 'refused'; readonly status: number }
  /** Another tab is downloading the same file. */
  | { readonly ok: false; readonly reason: 'busy' }
  /** The storage refused the bytes (it is full), or could not be written. */
  | { readonly ok: false; readonly reason: 'not-kept'; readonly detail: string };

export interface CountsDownloader {
  /** Rejects only when `signal` aborts (the bytes saved so far stay). */
  download(request: CountsRequest, onProgress: (progress: DownloadProgress) => void, signal: AbortSignal): Promise<CountsResult>;
}

/**
 * Deletions of saved counts between the tabs of the site: a tab that deletes a file says so (the
 * message carries only the taxonomy ID), and the others let go of the file they hold.
 */
export interface CountsDeletions {
  announce(taxid: number): void;
  /** Calls `listener` with each taxonomy ID that another tab deleted; returns the function that stops it. */
  subscribe(listener: (taxid: number) => void): () => void;
}
