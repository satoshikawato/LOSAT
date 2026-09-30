// BlockStore contract (plan §5.6): write, seal, discard, read, large blocks and a full
// storage. The memory store runs it in Vitest and in a browser worker; the OPFS store runs
// it in a browser worker (tests/e2e/contracts.spec.ts).
import type { BlockStore } from '../../src/infra/data/block-store';
import {
  concatBytes,
  MiB,
  pattern,
  patternByte,
  rejects,
  same,
  sameBytes,
  throws,
  type ContractCase,
} from './contract';

export interface BlockStoreEnv {
  /** A new, empty store. */
  readonly store: BlockStore;
  /** Makes the storage refuse more data: a memory budget, or a browser quota at the current usage. */
  exhaust(): Promise<void>;
  /** Undoes `exhaust`. */
  restore(): Promise<void>;
}

async function sealed(store: BlockStore, path: string, parts: readonly Uint8Array[]): Promise<Uint8Array> {
  const writer = await store.create(path);
  for (const part of parts) writer.append(part);
  await writer.seal();
  return concatBytes(parts);
}

export const BLOCK_STORE_CASES: readonly ContractCase<BlockStoreEnv>[] = [
  {
    name: 'write and seal: a block reads back exactly the bytes appended to it',
    async run({ store }) {
      const parts = [pattern(10, 1), pattern(0, 2), pattern(1000, 3), pattern(1, 4)];
      const writer = await store.create('runs/r1/out6');
      for (const part of parts) writer.append(part);
      same(await writer.seal(), 1011, 'sealed length');
      same(await store.size('runs/r1/out6'), 1011, 'size');
      sameBytes(await store.read('runs/r1/out6', 0, 1011), concatBytes(parts), 'bytes');
      same(store.usage(), 1011, 'usage');
    },
  },
  {
    name: 'read: any range of a sealed block, and nothing beyond its end',
    async run({ store }) {
      const bytes = await sealed(store, 'a/b', [pattern(2000, 5), pattern(3000, 6)]);
      const ranges: Array<[number, number]> = [
        [0, 0],
        [0, 1],
        [1999, 2],
        [2000, 3000],
        [4999, 1],
        [123, 4321],
        [5000, 0],
      ];
      for (const [offset, length] of ranges) {
        sameBytes(await store.read('a/b', offset, length), bytes.subarray(offset, offset + length), `range ${offset}+${length}`);
      }
      await rejects(store.read('a/b', 4999, 2), /beyond the end/, 'a range past the end');
      await rejects(store.read('a/b', -1, 1), /invalid range/, 'a negative offset');
      await rejects(store.read('a/missing', 0, 0), /does not exist/, 'a missing block');
    },
  },
  {
    name: 'write: the store keeps a copy, so the caller may reuse its buffer',
    async run({ store }) {
      const buffer = pattern(64, 7);
      const expected = buffer.slice();
      const writer = await store.create('a/copy');
      writer.append(buffer);
      buffer.fill(0);
      await writer.seal();
      sameBytes(await store.read('a/copy', 0, 64), expected, 'bytes after the caller reused its buffer');
    },
  },
  {
    name: 'staging: a block cannot be read before it is sealed',
    async run({ store }) {
      const writer = await store.create('a/staged');
      writer.append(pattern(100, 8));
      await rejects(store.read('a/staged', 0, 1), /not sealed/, 'read before seal');
      await rejects(store.size('a/staged'), /not sealed/, 'size before seal');
      same(store.usage(), 100, 'usage counts staged bytes');
      await writer.seal();
      same(await store.size('a/staged'), 100, 'size after seal');
    },
  },
  {
    name: 'seal: a sealed block takes no more bytes and cannot be sealed again',
    async run({ store }) {
      const writer = await store.create('a/once');
      writer.append(pattern(5, 9));
      await writer.seal();
      throws(() => writer.append(pattern(1, 1)), /not open/, 'append after seal');
      await rejects(writer.seal(), /not open/, 'a second seal');
      same(await store.size('a/once'), 5, 'size');
    },
  },
  {
    name: 'discard: a staged or sealed block disappears and its bytes are freed',
    async run({ store }) {
      const staged = await store.create('a/staged');
      staged.append(pattern(300, 1));
      const kept = await sealed(store, 'a/kept', [pattern(200, 2)]);
      const done = await store.create('a/sealed');
      done.append(pattern(100, 3));
      await done.seal();
      same(store.usage(), 600, 'usage before discard');
      await staged.discard();
      await done.discard();
      await done.discard();
      same(store.usage(), 200, 'usage after discard');
      await rejects(store.read('a/sealed', 0, 0), /does not exist/, 'read after discard');
      throws(() => staged.append(pattern(1, 1)), /not open/, 'append after discard');
      sameBytes(await store.read('a/kept', 0, 200), kept, 'the other block');
      const again = await store.create('a/sealed');
      again.append(pattern(7, 4));
      same(await again.seal(), 7, 'the path can be used again');
    },
  },
  {
    name: 'create: a path that exists cannot be created again',
    async run({ store }) {
      await sealed(store, 'a/x', [pattern(3, 1)]);
      await store.create('a/y');
      await rejects(store.create('a/x'), /already exists/, 'a sealed block');
      await rejects(store.create('a/y'), /already exists/, 'a staged block');
    },
  },
  {
    name: 'removeAll: removes one directory, sealed or not, and leaves the rest',
    async run({ store }) {
      await sealed(store, 'runs/r1/out6', [pattern(10, 1)]);
      const open = await store.create('runs/r1/hits');
      open.append(pattern(20, 2));
      const other = await sealed(store, 'runs/r2/out6', [pattern(30, 3)]);
      await store.removeAll('runs/r1');
      await rejects(store.read('runs/r1/out6', 0, 0), /does not exist/, 'a removed sealed block');
      throws(() => open.append(pattern(1, 1)), /not open/, 'a removed staged block');
      sameBytes(await store.read('runs/r2/out6', 0, 30), other, 'a block of another directory');
      same(store.usage(), 30, 'usage');
      await store.removeAll('runs/r1');
      await store.removeAll('runs/none');
      await sealed(store, 'runs/r1/out6', [pattern(4, 4)]);
      same(await store.size('runs/r1/out6'), 4, 'a removed path can be used again');
    },
  },
  {
    name: 'paths: empty and relative names are rejected',
    async run({ store }) {
      for (const path of ['', 'a//b', '../x', 'a/./b', 'a/']) {
        await rejects(store.create(path), /invalid block path/, `create "${path}"`);
      }
      await rejects(store.removeAll('..'), /invalid block path/, 'removeAll ".."');
    },
  },
  {
    name: 'an empty block can be sealed and read',
    async run({ store }) {
      const writer = await store.create('a/empty');
      same(await writer.seal(), 0, 'length');
      sameBytes(await store.read('a/empty', 0, 0), new Uint8Array(0), 'bytes');
    },
  },
  {
    name: 'a large block: 48 MiB in 1 MiB appends reads back unchanged',
    async run({ store }) {
      const writer = await store.create('runs/big/out0');
      for (let i = 0; i < 48; i++) writer.append(pattern(MiB, i));
      same(await writer.seal(), 48 * MiB, 'length');
      for (const [offset, length] of [
        [0, MiB],
        [MiB - 3, 7],
        [17 * MiB + 12345, 3 * MiB],
        [48 * MiB - 10, 10],
      ] as const) {
        const bytes = await store.read('runs/big/out0', offset, length);
        const expected = new Uint8Array(length);
        for (let i = 0; i < length; i++) expected[i] = patternByte((offset + i) % MiB, Math.floor((offset + i) / MiB));
        sameBytes(bytes, expected, `range ${offset}+${length}`);
      }
      const whole = await store.read('runs/big/out0', 0, 48 * MiB);
      for (let block = 0; block < 48; block++) {
        sameBytes(whole.subarray(block * MiB, (block + 1) * MiB), pattern(MiB, block), `MiB ${block}`);
      }
      await store.removeAll('runs/big');
      same(store.usage(), 0, 'usage after removal');
    },
  },
  {
    name: 'storage full: an append fails with StorageFullError, the block can only be discarded, and sealed blocks stay readable',
    async run(env) {
      const { store } = env;
      const kept = await sealed(store, 'kept/a', [pattern(4096, 9)]);
      const before = store.usage();
      await env.exhaust();
      try {
        const writer = await store.create('full/b');
        throws(() => writer.append(pattern(4 * MiB, 10)), /StorageFullError/, 'an append beyond the storage');
        throws(() => writer.append(pattern(1, 11)), /not open/, 'an append after the failure');
        await rejects(writer.seal(), /StorageFullError/, 'seal after the failure');
        await writer.discard();
        same(store.usage(), before, 'usage after the discard');
        sameBytes(await store.read('kept/a', 0, kept.length), kept, 'a sealed block while the storage is full');
      } finally {
        await env.restore();
      }
      const after = await store.create('full/c');
      after.append(pattern(1000, 12));
      same(await after.seal(), 1000, 'a block after the storage is restored');
    },
  },
];
