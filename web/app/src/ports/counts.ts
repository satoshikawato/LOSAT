// WindowMasker counts ports (line WW2, docs/web/fetch_design.md §10). The download runs on the page,
// only when the user asks, and its request carries only the file's URL (web/AGENTS.md rule 4):
// src/infra/browser/counts-downloader.ts. The body goes to the Data worker as it arrives, in parts
// of about 8 MiB, and is written at the end of the saved file (OPFS windowmasker/<taxid>/,
// src/infra/data/counts-store.ts): the file of a human is 724 MB, and is never held in memory.
// A download that stops (Cancel, a cut connection, a closed tab) leaves the bytes saved so far, and
// "Continue" asks for the rest with a Range. Tests use src/infra/fake/fake-counts.ts.
import type { CountsMadeMeta } from '../domain/counts';
import type { DownloadProgress } from './example';

/** Whether a saved copy of counts came from NCBI's site or from this browser's WindowMasker tab. */
export type SavedCountsOrigin = 'downloaded' | 'made';

/** What a saved file of counts records about itself (meta.json). */
export interface CountsCopy {
  readonly taxid: number;
  readonly url: string;
  /** The Last-Modified of the first answer, as the site sent it ('' when it sent none or for a made copy). */
  readonly lastModified: string;
  /** The size of the whole file. */
  readonly size: number;
  /** Bytes kept from the start of the file; Continue asks for the file from this byte. */
  readonly savedBytes: number;
  /** The whole file is kept. */
  readonly complete: boolean;
  /** When the last part was saved, in milliseconds since the epoch. */
  readonly savedAt: number;
  /** Present when the counts were made in this browser's WindowMasker tab; absent for an NCBI download. */
  readonly made?: CountsMadeMeta;
}

export interface SavedCounts {
  readonly copy: CountsCopy;
  readonly origin: SavedCountsOrigin;
  /** The kept file, read from the browser's storage, not from memory; only when the copy is complete. */
  readonly file?: File;
}

/** One complete saved counts file listed for the WindowMasker tab's pickers (`listCounts`). */
export interface SavedCountsEntry {
  readonly taxid: number;
  /** The file name under windowmasker/<taxid>/ (wmasker.obinary). */
  readonly name: string;
  readonly size: number;
  readonly origin: SavedCountsOrigin;
  /** ISO date (yyyy-mm-dd) of the copy: Last-Modified for a download, madeAt for a made copy, or savedAt. */
  readonly date: string;
  readonly copy: CountsCopy;
}

/** What a part is appended to: a new or continued download, or a copy made in this browser written whole. */
export type CountsStart =
  | {
      readonly url: string;
      readonly lastModified: string;
      readonly size: number;
      /** 0 begins the file again; otherwise the bytes already saved (the saved part is checked against it). */
      readonly from: number;
      readonly made?: undefined;
    }
  | {
      readonly made: CountsMadeMeta;
      readonly size: number;
      /** A made copy is always written whole from byte 0. */
      readonly from: 0;
      /** When false, refuses with COUNTS_EXISTS if a complete or partial copy is already stored for the ID. */
      readonly replace?: boolean;
    };

/** The name of the error that `beginCounts` or `deleteCounts` rejects with when another tab writes the same file. */
export const COUNTS_BUSY = 'CountsBusyError';
/** The name of the error that `beginCounts` or `deleteCounts` rejects with when a queued or running search uses the same file. */
export const COUNTS_IN_USE = 'CountsInUseError';
/** The name of the error that `beginCounts` rejects with when the saved part is not the one asked to continue. */
export const COUNTS_CHANGED = 'CountsChangedError';
/** The name of the error that `beginCounts` rejects with when `replace` is false and a copy is already stored. */
export const COUNTS_EXISTS = 'CountsExistsError';

/**
 * Counts kept across sessions and tabs, so that a file is downloaded or made once. A file is written
 * by one tab at a time (a Web Lock held from `beginCounts` to `finishCounts` or `closeCounts`, or
 * until the tab ends). Without storage that lasts (no OPFS or no Web Locks), nothing is kept.
 */
export interface CountsStore {
  /**
   * The saved copy of the taxonomy ID (complete or partial download, or complete made copy) and its
   * origin. An interrupted made copy (`complete: false`) is removed when read.
   */
  readCounts(taxid: number): Promise<SavedCounts | undefined>;
  /** Every complete copy kept in this browser, ordered by taxonomy ID. */
  listCounts(): Promise<readonly SavedCountsEntry[]>;
  /**
   * Begins writing the file and returns the token of the writing; undefined when this browser
   * keeps nothing across sessions. Rejects with COUNTS_BUSY, COUNTS_CHANGED or COUNTS_EXISTS, or
   * with the storage's refusal.
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
  /** A queued or running search uses the saved counts for this taxonomy ID. */
  | { readonly ok: false; readonly reason: 'in-use' }
  /** Another tab changed the saved part between its look-up and the download. */
  | { readonly ok: false; readonly reason: 'changed' }
  /** The storage refused the bytes (it is full), or could not be written. */
  | { readonly ok: false; readonly reason: 'not-kept'; readonly detail: string };

export interface CountsDownloader {
  /** Rejects only when `signal` aborts (the bytes saved so far stay). */
  download(request: CountsRequest, onProgress: (progress: DownloadProgress) => void, signal: AbortSignal): Promise<CountsResult>;
}

/**
 * Holds a shared Web Lock on a taxonomy ID's saved counts while a queued or running run reads
 * them, so Replace, Delete and download refuse in place until the run ends.
 */
export interface CountsLocks {
  /**
   * Acquires the taxonomy ID's lock in shared mode (`ifAvailable: true`) and returns a release
   * function that lets it go once. Returns `undefined` when an exclusive writer holds the lock.
   */
  holdShared(taxid: number): Promise<(() => void) | undefined>;
}

/**
 * Changes to saved counts between the tabs of the site: a tab that saves, replaces or deletes a
 * file says so (the message carries only the taxonomy ID), and the others read that ID again.
 */
export interface CountsChanges {
  announce(taxid: number): void;
  /** Calls `listener` with each taxonomy ID that another tab changed; returns the function that stops it. */
  subscribe(listener: (taxid: number) => void): () => void;
}
