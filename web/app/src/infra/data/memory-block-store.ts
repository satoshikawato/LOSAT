// Memory implementation of the BlockStore contract, used when OPFS is not available
// (plan §5.6). It replaces the S01 MemoryDataGateway.
import { checkRange, pathNames, StorageFullError, type BlockStore, type BlockWriter } from './block-store';

interface MemoryBlock {
  readonly chunks: Uint8Array[];
  length: number;
  state: 'open' | 'sealed' | 'full';
}

export interface MemoryBlockStoreOptions {
  /** Memory budget in bytes; unlimited by default (S09 measures the memory limits). */
  readonly capacityBytes?: number;
}

export class MemoryBlockStore implements BlockStore {
  readonly backend = 'memory';
  private readonly blocks = new Map<string, MemoryBlock>();
  private used = 0;
  private capacity: number;

  constructor(options: MemoryBlockStoreOptions = {}) {
    this.capacity = options.capacityBytes ?? Number.POSITIVE_INFINITY;
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
    return {
      append: (bytes) => {
        if (!live() || block.state !== 'open') throw new Error(`block ${path} is not open`);
        if (this.used + bytes.length > this.capacity) {
          block.state = 'full';
          throw new StorageFullError();
        }
        block.chunks.push(bytes.slice());
        block.length += bytes.length;
        this.used += bytes.length;
      },
      seal: async () => {
        if (live() && block.state === 'full') throw new StorageFullError();
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
