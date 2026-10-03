// Block storage of the Data worker (plan §5.6): a block is appended to, sealed, and then
// only read. Two implementations pass the same contract
// (tests/contract/block-store.contract.ts): OPFS (opfs-block-store.ts), used when the
// browser supports it, and memory (memory-block-store.ts).
import type { StorageBackend } from '../../ports/data';

export interface BlockWriter {
  /**
   * Appends a copy of `bytes`. Throws StorageFullError when the storage is full; the block
   * can then only be discarded.
   */
  append(bytes: Uint8Array): void;
  /** Seals the block: it takes no more bytes and becomes readable. Resolves to its length. */
  seal(): Promise<number>;
  /** Drops the block, sealed or not. */
  discard(): Promise<void>;
}

export interface BlockStore {
  readonly backend: StorageBackend;
  /** Starts a new, empty block. `path` is `/`-separated names; it must not exist yet. */
  create(path: string): Promise<BlockWriter>;
  /** Reads `length` bytes at `offset` of a sealed block, as a new buffer. */
  read(path: string, offset: number, length: number): Promise<Uint8Array>;
  /** The length of a sealed block. */
  size(path: string): Promise<number>;
  /** Removes every block under the directory `prefix`, sealed or not. */
  removeAll(prefix: string): Promise<void>;
  /** Bytes held in blocks, sealed or not. */
  usage(): number;
}

export const STORAGE_FULL_MESSAGE =
  'Not enough temporary storage for the results of this run: the browser refused to store more data. ' +
  'Earlier results are kept.';

export class StorageFullError extends Error {
  constructor(message = STORAGE_FULL_MESSAGE) {
    super(message);
    this.name = 'StorageFullError';
  }
}

/**
 * Results kept in memory (no OPFS) may take at most this many bytes in a tab (S09). Memory
 * is then the only copy of the results, so running out of it would end the tab and lose
 * every result; the budget makes it the same recoverable failure as a full OPFS instead.
 * The engine's own memory comes on top (up to 1 GiB for the threaded module, plan TD-7).
 */
export const MEMORY_RESULTS_CAPACITY_BYTES = 512 * 1024 * 1024;

export const MEMORY_FULL_MESSAGE =
  'Not enough temporary storage for the results of this run: this browser keeps results in memory, ' +
  'at most 512 MB in a tab. Earlier results are kept.';

/** The StorageFullError of a storage-full error, keeping its message (memory or browser storage). */
export function asStorageFull(error: unknown): StorageFullError {
  return error instanceof StorageFullError ? error : new StorageFullError();
}

/** True for StorageFullError and for the browser's QuotaExceededError. */
export function isStorageFull(error: unknown): boolean {
  return (
    error instanceof StorageFullError ||
    (typeof error === 'object' && error !== null && (error as { name?: unknown }).name === 'QuotaExceededError')
  );
}

/** Splits a block path into its names; rejects empty and relative names. */
export function pathNames(path: string): string[] {
  const names = path.split('/');
  if (names.some((name) => name === '' || name === '.' || name === '..')) {
    throw new Error(`invalid block path "${path}"`);
  }
  return names;
}

/** Checks a read range against a block length. */
export function checkRange(path: string, offset: number, length: number, size: number): void {
  if (!Number.isSafeInteger(offset) || !Number.isSafeInteger(length) || offset < 0 || length < 0) {
    throw new RangeError(`invalid range ${offset}+${length} of block ${path}`);
  }
  if (offset + length > size) {
    throw new RangeError(`range ${offset}+${length} is beyond the end of block ${path} (${size} bytes)`);
  }
}
