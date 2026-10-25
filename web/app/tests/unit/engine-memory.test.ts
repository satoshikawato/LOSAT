// The engine's memory above 2 GiB and the maximum of the shared memory (plan TD-7,
// docs/web/abi_v2.md §3 and §6): the shared memory steps down to a maximum that the
// browser provides; the ABI binding and the WASI imports take pointers as unsigned values.
// With LOSAT_WEB_REACTORS, the serial reactor also runs one search per program with its
// heap above 2 GiB and gives the bytes of an ordinary instance.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { OUTPUT_FORMATS } from '../../src/domain/output-format';
import { ReactorAbi, ROLE_QUERY, ROLE_SUBJECT, type AbiExports } from '../../src/infra/reactor/abi';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import { candidateMaximums, createSharedMemory, STEP_DOWN_PAGES } from '../../src/infra/reactor/shared-memory';
import { createWasi, unsignedImports } from '../../src/infra/reactor/wasi';

const PAGE = 65_536;
const TWO_GIB = 2 ** 31;
/** A Wasm i32 as JavaScript receives it. */
const signed = (value: number) => value | 0;

describe('createSharedMemory', () => {
  const DECLARED = { initial: 23, maximum: 65_535, shared: true };

  /**
   * An engine that provides at most `limit` pages. `refuses`: `new WebAssembly.Memory`
   * throws above the limit (JavaScriptCore). Otherwise the memory is created with its
   * initial pages only and cannot grow (V8 without guard regions).
   */
  function engine(limit: number, refuses: boolean) {
    const tried: number[] = [];
    const create = (descriptor: WebAssembly.MemoryDescriptor) => {
      tried.push(descriptor.maximum!);
      if (descriptor.maximum! > limit && refuses) throw new RangeError('Out of memory');
      let pages = descriptor.initial;
      return {
        descriptor,
        get pages() {
          return pages;
        },
        grow(delta: number) {
          if (descriptor.maximum! > limit) throw new RangeError('Unable to grow instance memory');
          pages += delta;
          return pages - delta;
        },
      } as unknown as WebAssembly.Memory;
    };
    return { tried, create };
  }
  const pagesOf = (memory: WebAssembly.Memory) => (memory as unknown as { pages: number }).pages;
  const descriptorOf = (memory: WebAssembly.Memory) => (memory as unknown as { descriptor: WebAssembly.MemoryDescriptor }).descriptor;

  it('tries the declared maximum, then 3, 2, 1.5 and 1 GiB', () => {
    expect(candidateMaximums(DECLARED)).toEqual([65_535, 49_152, 32_768, 24_576, 16_384]);
    expect(STEP_DOWN_PAGES).toEqual([49_152, 32_768, 24_576, 16_384]);
    // A module built before the limit was lifted: its 1 GiB only, as before.
    expect(candidateMaximums({ ...DECLARED, maximum: 16_384 })).toEqual([16_384]);
    // Never less than the module's initial size.
    expect(candidateMaximums({ ...DECLARED, initial: 30_000 })).toEqual([65_535, 49_152, 32_768]);
    expect(() => candidateMaximums({ initial: 23, shared: true })).toThrow(/declares no maximum/);
  });

  for (const refuses of [true, false]) {
    describe(refuses ? 'an engine that refuses a maximum it cannot reserve' : 'an engine that creates a memory that cannot grow', () => {
      it('uses the declared maximum when the engine provides it', () => {
        const { tried, create } = engine(65_535, refuses);
        const { memory, maximumPages } = createSharedMemory(DECLARED, create);
        expect(maximumPages).toBe(65_535);
        expect(tried).toEqual([65_535]);
        // One page short of the initial size, then that page: the module finds its initial size.
        expect(descriptorOf(memory)).toEqual({ initial: 22, maximum: 65_535, shared: true });
        expect(pagesOf(memory)).toBe(23);
      });

      it('steps down to the largest maximum that the engine provides', () => {
        const { tried, create } = engine(40_000, refuses);
        const { memory, maximumPages } = createSharedMemory(DECLARED, create);
        expect(maximumPages).toBe(32_768);
        expect(tried).toEqual([65_535, 49_152, 32_768]);
        expect(pagesOf(memory)).toBe(23);
      });

      it('fails below 1 GiB, so that the search runs on the serial module as before', () => {
        const { tried, create } = engine(16_383, refuses);
        expect(() => createSharedMemory(DECLARED, create)).toThrow(/^the browser could not reserve 1 GiB of shared memory \(/);
        expect(tried).toEqual([65_535, 49_152, 32_768, 24_576, 16_384]);
        const old = engine(16_383, refuses);
        expect(() => createSharedMemory({ ...DECLARED, maximum: 16_384 }, old.create)).toThrow(/could not reserve 1 GiB/);
        expect(old.tried).toEqual([16_384]);
      });
    });
  }

  it('creates a real shared memory of the initial size', () => {
    const { memory, maximumPages } = createSharedMemory({ initial: 3, maximum: 8, shared: true });
    expect(maximumPages).toBe(8);
    expect(memory.buffer).toBeInstanceOf(SharedArrayBuffer);
    expect(memory.buffer.byteLength).toBe(3 * PAGE);
    expect(memory.grow(5)).toBe(3);
    expect(() => memory.grow(1)).toThrow(RangeError);
  });
});

describe('pointers of 2 GiB or more', () => {
  /** 2 GiB and four pages: untouched pages take address space only. */
  const bigMemory = () => new WebAssembly.Memory({ initial: TWO_GIB / PAGE + 4 });

  it('passes i32 arguments of the WASI imports on as unsigned values, and i64 arguments unchanged', () => {
    const seen: unknown[][] = [];
    const imports = unsignedImports({
      clock_time_get: (...args: unknown[]) => {
        seen.push(args);
        return 0;
      },
      proc_exit: (...args: unknown[]) => void seen.push(args),
      not_a_function: 1,
    }) as Record<string, (...args: Array<number | bigint>) => unknown>;
    expect(Object.keys(imports)).toEqual(['clock_time_get', 'proc_exit']);
    expect(imports['clock_time_get']!(1, 1n, signed(0x8000_0008))).toBe(0);
    imports['proc_exit']!(-1);
    expect(seen).toEqual([[1, 1n, 0x8000_0008], [0xffff_ffff]]);
  });

  it('the WASI shim writes the clock, random bytes and standard error through high pointers', () => {
    const memory = bigMemory();
    const wasi = createWasi();
    wasi.bind({ exports: { memory } } as unknown as WebAssembly.Instance);
    const imports = wasi.imports as Record<string, (...args: Array<number | bigint>) => number>;
    const bytes = () => new Uint8Array(memory.buffer);
    const view = () => new DataView(memory.buffer);
    const time = TWO_GIB + 8;
    expect(imports['clock_time_get']!(1, 1n, signed(time))).toBe(0);
    expect(view().getBigUint64(time, true)).toBeGreaterThan(0n);

    // random_get: 64 bytes at the pointer, and nothing where a negative start would land.
    const random = TWO_GIB + 4096;
    const wrong = memory.buffer.byteLength + signed(random);
    imports['random_get']!(signed(random), 64);
    expect(bytes().subarray(random, random + 64).some((b) => b !== 0)).toBe(true);
    expect(bytes().subarray(wrong, wrong + 64).every((b) => b === 0)).toBe(true);

    // fd_write to standard error: the iovec, its buffer and the result are all high.
    const text = new TextEncoder().encode('memory allocation of 8 bytes failed\n');
    const iovec = TWO_GIB + 2 * PAGE;
    const buffer = iovec + 64;
    const written = iovec + 32;
    bytes().set(text, buffer);
    view().setUint32(iovec, buffer, true);
    view().setUint32(iovec + 4, text.length, true);
    expect(imports['fd_write']!(2, signed(iovec), 1, signed(written))).toBe(0);
    expect(view().getUint32(written, true)).toBe(text.length);
    expect(wasi.output()).toBe('memory allocation of 8 bytes failed\n');
  });

  it('the ABI binding writes inputs to, and reads responses and errors from, high addresses', () => {
    const memory = bigMemory();
    const input = TWO_GIB + 1024;
    const response = TWO_GIB + PAGE;
    const failure = TWO_GIB + 2 * PAGE;
    const json = new TextEncoder().encode('{"handle": 7, "records": [{"id": "a", "length": 4}]}');
    const message = new TextEncoder().encode('BLAST query error');
    const calls: number[][] = [];
    const ref: { abi?: ReactorAbi } = {};
    const exports = {
      memory,
      losat_web2_alloc: () => signed(input),
      losat_web2_dealloc: (ptr: number, len: number) => {
        calls.push([ptr, len]);
        return 0;
      },
      losat_web2_register: (programPtr: number, programLen: number, role: number, ptr: number, len: number) => {
        calls.push([programPtr, programLen, role, ptr, len]);
        new Uint8Array(memory.buffer, response, json.length).set(json);
        ref.abi!.emit(2, signed(response), json.length);
        return 7;
      },
      losat_web2_release: () => 0,
      losat_web2_validate: () => -1,
      losat_web2_last_error_ptr: () => {
        new Uint8Array(memory.buffer, failure, message.length).set(message);
        return signed(failure);
      },
      losat_web2_last_error_len: () => message.length,
    } as unknown as AbiExports;
    ref.abi = new ReactorAbi(exports);
    expect(signed(input)).toBeLessThan(0);
    expect(ref.abi.register('blastn', 1, new Uint8Array([62, 97, 10, 65]))).toEqual({ handle: 7, records: [{ id: 'a', length: 4 }] });
    // The input bytes are at the unsigned address (the last argument written stays there).
    expect([...new Uint8Array(memory.buffer, input, 4)]).toEqual([62, 97, 10, 65]);
    // The engine receives the pointers it gave; as i32 they are the same bits.
    expect(calls[0]!.map(signed)).toEqual([signed(input), 6, 1, signed(input), 4]);
    expect(ref.abi.validate(['blastn'])).toBe('BLAST query error');
  });
});

// --- the built reactors (LOSAT_WEB_REACTORS) ---------------------------------------------------

const reactors = findReactors();
const LOSAT = resolve(import.meta.dirname, '../../../../LOSAT');
const sha256 = (bytes: Uint8Array) => createHash('sha256').update(bytes).digest('hex');

describe.skipIf(reactors === undefined)('the serial reactor with its heap above 2 GiB', () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);

  /**
   * An instance as instantiateSerial makes it, whose memory then grows past 2 GiB and whose
   * free memory below that is used up, so that every allocation of the search is high.
   */
  async function highInstance(): Promise<{ abi: ReactorAbi; exports: AbiExports }> {
    const wasi = createWasi();
    const bound: { abi?: ReactorAbi } = {};
    const instance = await WebAssembly.instantiate(module!, {
      losat_host: { emit: (stream: number, ptr: number, len: number) => bound.abi!.emit(stream, ptr, len) },
      wasi_snapshot_preview1: wasi.imports,
    });
    const exports = instance.exports as unknown as AbiExports;
    const abi = new ReactorAbi(exports, () => wasi.output());
    bound.abi = abi;
    wasi.initialize(instance);
    exports.memory.grow(TWO_GIB / PAGE + 1024 - exports.memory.buffer.byteLength / PAGE);
    for (const size of [PAGE, 4096, 256, 16]) {
      for (let i = 0; i < 100_000; i++) {
        const ptr = exports.losat_web2_alloc(size) >>> 0;
        if (ptr === 0) throw new Error(`the engine could not allocate ${size} bytes`);
        if (ptr >= TWO_GIB) {
          exports.losat_web2_dealloc(ptr, size);
          break;
        }
      }
    }
    return { abi, exports };
  }

  /** One search per program of LOSAT Web: V-BR cells of LOSAT/tests/outfmt0_manifest.tsv (run from LOSAT/). */
  const CASES = [
    { id: 'compact.blastn', argv: ['blastn', '-task', 'blastn'], query: 'tests/fasta/blastn_parity_compact.fasta', subject: 'tests/fasta/blastn_parity_compact.fasta' },
    { id: 'width.blastp', argv: ['blastp'], query: 'tests/fasta/outfmt0/pwidth_q.faa', subject: 'tests/fasta/outfmt0/pwidth_s.faa' },
    { id: 'width.tblastn', argv: ['tblastn'], query: 'tests/fasta/outfmt0/twidth_q.faa', subject: 'tests/fasta/outfmt0/twidth_s.fna' },
    { id: 'tblastx.multi', argv: ['tblastx'], query: 'tests/fasta/outfmt0/tblastx_multi_query.fasta', subject: 'tests/fasta/outfmt0/tblastx_multi_subject.fasta' },
  ] as const;

  /** SHA-256 and length of every stream of one search. */
  function search(abi: ReactorAbi, testCase: (typeof CASES)[number]): Record<number, { sha256: string; bytes: number }> {
    const [program, ...options] = testCase.argv;
    const query = abi.register(program, ROLE_QUERY, readFileSync(resolve(LOSAT, testCase.query)));
    const subject = abi.register(program, ROLE_SUBJECT, readFileSync(resolve(LOSAT, testCase.subject)));
    const streams = new Map<number, Uint8Array[]>();
    const argv = [program, '-query', testCase.query, '-subject', testCase.subject, ...options, '-num_threads', '1'];
    abi.run(argv, query.handle, subject.handle, (stream, bytes) => streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]));
    abi.release(query.handle);
    abi.release(subject.handle);
    const out: Record<number, { sha256: string; bytes: number }> = {};
    for (const [stream, chunks] of streams) {
      const joined = Buffer.concat(chunks);
      out[stream] = { sha256: sha256(joined), bytes: joined.length };
    }
    return out;
  }

  for (const testCase of CASES) {
    it(`${testCase.id}: every stream is the ordinary instance's`, async () => {
      const ordinary = search((await instantiateSerial(module!)).abi, testCase);
      const { abi, exports } = await highInstance();
      expect(exports.losat_web2_alloc(16) >>> 0).toBeGreaterThanOrEqual(TWO_GIB);
      const high = search(abi, testCase);
      for (const format of OUTPUT_FORMATS) expect(ordinary[format]?.bytes, `${testCase.id} outfmt ${format}`).toBeGreaterThan(0);
      expect(high).toEqual(ordinary);
      expect(abi.memoryBytes()).toBeGreaterThan(TWO_GIB);
    }, 60_000);
  }
});
