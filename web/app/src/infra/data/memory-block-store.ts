// Memory implementation of the BlockStore contract, used when OPFS is not available
// (plan §5.6). It replaces the S01 MemoryDataGateway.
import { checkRange, pathNames, STORAGE_FULL_MESSAGE, StorageFullError, type BlockStore, type BlockWriter } from './block-store';

interface MemoryBlock {
  readonly chunks: Uint8Array[];
  length: number;
  state: 'open' | 'sealed' | 'full';
}

export interface MemoryBlockStoreOptions {
  /** Memory budget in bytes; unlimited by default (the application sets MEMORY_RESULTS_CAPACITY_BYTES of ports/data.ts). */
  readonly capacityBytes?: number;
  /** The message of the StorageFullError beyond the budget. */
  readonly fullMessage?: string;
}

export class MemoryBlockStore implements BlockStore {
  readonly backend = 'memory';
  private readonly blocks = new Map<string, MemoryBlock>();
  private used = 0;
  private capacity: number;
  private readonly fullMessage: string;

  constructor(options: MemoryBlockStoreOptions = {}) {
    this.capacity = options.capacityBytes ?? Number.POSITIVE_INFINITY;
    this.fullMessage = options.fullMessage ?? STORAGE_FULL_MESSAGE;
  }

  /** The memory budget in bytes. */
  capacityBytes(): number {
    return this.capacity;
  }

  /** Changes the memory budget; appends beyond it fail with StorageFullError. */
  setCapacity(bytes: number): void {
    this.capacity = bytes;
  }

  async create(path: string): Promise<BlockWriter> {
    pathNames(path);
    if (this.blocks.has(path)) throw new Error(`block ${path} already exists`);
    const block: MemoryBlock = { chunks: [], length: 0, state: 'open' };
    this.blocks.set(path, block);
    const live = () => this.blocks.get(path) === block;
    const append = (bytes: Uint8Array) => {
      if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
      if (this.used + bytes.length > this.capacity) {
        block.state = 'full';
        throw new StorageFullError(this.fullMessage);
      }
      if (bytes.length === 0) return;
      block.chunks.push(bytes.slice());
      block.length += bytes.length;
      this.used += bytes.length;
    };
    const writeAt = (offset: number, bytes: Uint8Array) => {
      if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
      if (!Number.isSafeInteger(offset) || offset < 0) {
        throw new RangeError(`invalid write offset ${offset} of block ${path}`);
      }
      if (bytes.length === 0) return;
      if (offset === block.length) {
        append(bytes);
        return;
      }
      const nextLength = Math.max(block.length, offset + bytes.length);
      const delta = nextLength - block.length;
      if (this.used + delta > this.capacity) {
        block.state = 'full';
        throw new StorageFullError(this.fullMessage);
      }
      const merged = new Uint8Array(nextLength);
      let cursor = 0;
      for (const chunk of block.chunks) {
        merged.set(chunk, cursor);
        cursor += chunk.length;
      }
      merged.set(bytes, offset);
      block.chunks.length = 0;
      block.chunks.push(merged);
      block.length = nextLength;
      this.used += delta;
    };
    const clear = () => {
      if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
      this.used -= block.length;
      block.chunks.length = 0;
      block.length = 0;
    };
    const truncate = (length: number) => {
      if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
      if (!Number.isSafeInteger(length) || length < 0) {
        throw new RangeError(`invalid truncate length ${length} of block ${path}`);
      }
      if (length === 0) {
        clear();
        return;
      }
      if (length === block.length) return;
      if (length > block.length) {
        const delta = length - block.length;
        if (this.used + delta > this.capacity) {
          block.state = 'full';
          throw new StorageFullError(this.fullMessage);
        }
        block.chunks.push(new Uint8Array(delta));
        block.length = length;
        this.used += delta;
        return;
      }
      const kept: Uint8Array[] = [];
      let remain = length;
      for (const chunk of block.chunks) {
        if (remain <= 0) break;
        if (chunk.length <= remain) {
          kept.push(chunk);
          remain -= chunk.length;
        } else {
          kept.push(chunk.slice(0, remain));
          remain = 0;
        }
      }
      block.chunks.length = 0;
      for (const chunk of kept) block.chunks.push(chunk);
      this.used -= block.length - length;
      block.length = length;
    };
    return {
      append,
      writeAt,
      truncate,
      clear,
      seal: async () => {
        if (live() && block.state === 'full') throw new StorageFullError(this.fullMessage);
        if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
        block.state = 'sealed';
        return block.length;
      },
      discard: async () => {
        if (live()) this.drop(path, block);
      },
    };
  }

  async read(path: string, offset: number, length: number): Promise<Uint8Array> {
    const block = this.sealed(path);
    checkRange(path, offset, length, block.length);
    const bytes = new Uint8Array(length);
    let chunkStart = 0;
    let filled = 0;
    for (const chunk of block.chunks) {
      if (filled === length) break;
      const chunkEnd = chunkStart + chunk.length;
      const from = Math.max(offset, chunkStart);
      const to = Math.min(offset + length, chunkEnd);
      if (from < to) {
        bytes.set(chunk.subarray(from - chunkStart, to - chunkStart), from - offset);
        filled += to - from;
      }
      chunkStart = chunkEnd;
    }
    return bytes;
  }

  async readFile(path: string, name: string, type = 'application/octet-stream'): Promise<File> {
    const block = this.sealed(path);
    return new File(block.chunks as BlobPart[], name, { type });
  }

  async size(path: string): Promise<number> {
    return this.sealed(path).length;
  }

  async removeAll(prefix: string): Promise<void> {
    pathNames(prefix);
    for (const [path, block] of [...this.blocks]) {
      if (path.startsWith(`${prefix}/`)) this.drop(path, block);
    }
  }

  usage(): number {
    return this.used;
  }

  private sealed(path: string): MemoryBlock {
    const block = this.blocks.get(path);
    if (block === undefined) throw new Error(`block ${path} does not exist`);
    if (block.state !== 'sealed') throw new Error(`block ${path} is not sealed`);
    return block;
  }

  private drop(path: string, block: MemoryBlock): void {
    this.blocks.delete(path);
    this.used -= block.length;
  }
}
