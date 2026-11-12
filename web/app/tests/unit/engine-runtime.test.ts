// Pure parts of the engine runtime: the thread and renewal policy, the check of an engine
// module before it is compiled, and the message of an instance that stopped. With
// LOSAT_WEB_REACTORS, the module check also runs on the built reactors, and the
// RecordScanner contract and the inflation of compressed sources run against the serial
// reactor in Node.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { constants as zlib, gzipSync, zstdCompressSync } from 'node:zlib';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { findReactors } from '../../build/reactors';
import { EngineStoppedError, ReactorAbi, stopMessage, type AbiExports } from '../../src/infra/reactor/abi';
import { ABI_EXPORTS, inspectReactor } from '../../src/infra/reactor/artifact';
import { instantiateSerial, loadReactor } from '../../src/infra/reactor/instance';
import { ReactorInputChecker } from '../../src/infra/reactor/checker';
import { ReactorInflater } from '../../src/infra/reactor/inflater';
import { inflateFile } from '../../src/infra/data/inflate';
import { INFLATED_TOO_LARGE } from '../../src/ports/inflate';
import { toProgramDescription } from '../../src/infra/reactor/control';
import { ReactorScanner, reopening } from '../../src/infra/reactor/scanner';
import type { FastaParserKind } from '../../src/domain/dataset';
import { PROGRAMS } from '../../src/domain/programs';
import {
  AUTO_MANY_THREADS_FROM_BYTES,
  AUTO_MAX_THREADS,
  AUTO_SERIAL_BELOW_BYTES,
  AUTO_SMALL_MAX_THREADS,
  chooseThreads,
  DEFAULT_RENEWAL,
  renewalReason,
} from '../../src/infra/engine-worker/policy';
import { InstanceUse } from '../../src/infra/engine-worker/instance-use';
import { EngineRuntime } from '../../src/infra/engine-worker/runtime';
import { warmUp, warmUpInputs } from '../../src/infra/engine-worker/warm-up';

// The runtime's calls of the warm-up are counted; the warm-up itself runs as it is.
vi.mock('../../src/infra/engine-worker/warm-up', async (importOriginal) => {
  const actual = await importOriginal<typeof import('../../src/infra/engine-worker/warm-up')>();
  return { ...actual, warmUp: vi.fn(actual.warmUp) };
});
import { RECORD_SCANNER_CASES } from '../contract/record-scanner.contract';

describe('chooseThreads', () => {
  it('keeps an explicit request', () => {
    expect(chooseThreads(1, 10 ** 9, 16)).toBe(1);
    expect(chooseThreads(4, 10, 2)).toBe(4);
  });

  it('runs Auto serially below the size threshold', () => {
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES - 1, 16)).toBe(1);
  });

  it('runs Auto with up to half the logical processors from the threshold on, within the maximum', () => {
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 8)).toBe(4);
    expect(chooseThreads('auto', AUTO_MANY_THREADS_FROM_BYTES, 32)).toBe(16);
    expect(chooseThreads('auto', AUTO_MANY_THREADS_FROM_BYTES, 64)).toBe(AUTO_MAX_THREADS);
    expect(chooseThreads('auto', AUTO_MANY_THREADS_FROM_BYTES, 256)).toBe(AUTO_MAX_THREADS);
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 1)).toBe(1);
    // A browser that does not report the processors.
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 0)).toBe(1);
  });

  it('runs Auto on fewer threads below the size from which many threads pay', () => {
    expect(chooseThreads('auto', AUTO_MANY_THREADS_FROM_BYTES - 1, 64)).toBe(AUTO_SMALL_MAX_THREADS);
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 64)).toBe(AUTO_SMALL_MAX_THREADS);
    expect(chooseThreads('auto', AUTO_MANY_THREADS_FROM_BYTES - 1, 8)).toBe(4);
  });
});

describe('renewalReason', () => {
  it('renews at the high-water mark of the linear memory', () => {
    expect(renewalReason({ linearBytesAfter: DEFAULT_RENEWAL.highWaterBytes - 1, instanceRuns: 1 }, DEFAULT_RENEWAL)).toBeUndefined();
    expect(renewalReason({ linearBytesAfter: DEFAULT_RENEWAL.highWaterBytes, instanceRuns: 1 }, DEFAULT_RENEWAL)).toMatch(/memory reached/);
  });

  it('renews after a maximum number of searches when one is set, and by default never for the count', () => {
    const limits = { ...DEFAULT_RENEWAL, maxRuns: 3 };
    expect(renewalReason({ linearBytesAfter: 0, instanceRuns: 2 }, limits)).toBeUndefined();
    expect(renewalReason({ linearBytesAfter: 0, instanceRuns: 3 }, limits)).toMatch(/ran 3 searches/);
    expect(renewalReason({ linearBytesAfter: 0, instanceRuns: 1_000_000 }, DEFAULT_RENEWAL)).toBeUndefined();
  });
});

// --- a minimal module encoder for inspectReactor ---------------------------------------------

const uleb = (value: number): number[] => {
  const out: number[] = [];
  do {
    const byte = value & 0x7f;
    value = Math.floor(value / 128);
    out.push(byte | (value ? 0x80 : 0));
  } while (value);
  return out;
};
const name = (text: string) => {
  const bytes = [...new TextEncoder().encode(text)];
  return [...uleb(bytes.length), ...bytes];
};
const vector = (items: number[][]) => [...uleb(items.length), ...items.flat()];
const section = (id: number, bytes: number[]) => [id, ...uleb(bytes.length), ...bytes];

interface ModuleSpec {
  readonly exports?: readonly string[];
  /** Imported function names (module.field) of type () -> i32. */
  readonly functionImports?: readonly string[];
  /** env.memory: [flags, initial, maximum?] with flag bit 0 = has maximum, bit 1 = shared. */
  readonly importedMemory?: readonly number[];
  readonly ownMemory?: readonly number[];
  /** What every function returns (0 to 63); 0 by default. */
  readonly returns?: number;
}

/** A module whose functions all have type () -> i32 and return 0 (or `returns`). */
function moduleBytes(spec: ModuleSpec): Uint8Array {
  const functionImports = spec.functionImports ?? [];
  const exports = spec.exports ?? [];
  const functions = exports.filter((e) => e !== 'memory');
  const imports = functionImports.map((full) => {
    const [module, field] = full.split('.');
    return [...name(module!), ...name(field!), 0, 0];
  });
  if (spec.importedMemory !== undefined) imports.push([...name('env'), ...name('memory'), 2, ...spec.importedMemory.flatMap(uleb)]);
  const offset = functionImports.length;
  return new Uint8Array([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 0, 1, 0x7f]),
    ...(imports.length > 0 ? section(2, vector(imports)) : []),
    ...section(3, vector(functions.map(() => [0]))),
    ...(spec.ownMemory !== undefined ? section(5, [1, ...spec.ownMemory.flatMap(uleb)]) : []),
    ...section(
      7,
      vector(exports.map((e) => (e === 'memory' ? [...name(e), 2, 0] : [...name(e), 0, ...uleb(offset + functions.indexOf(e))]))),
    ),
    ...section(10, vector(functions.map(() => [4, 0, 0x41, spec.returns ?? 0, 0x0b]))),
  ]);
}

const SERIAL: ModuleSpec = { exports: ABI_EXPORTS, ownMemory: [1, 17, 16384] };
const THREADED: ModuleSpec = {
  exports: [...ABI_EXPORTS, 'wasi_thread_start'],
  functionImports: ['wasi.thread-spawn'],
  importedMemory: [3, 17, 16384],
};

describe('inspectReactor', () => {
  it('finds a serial reactor and its own memory', () => {
    const bytes = moduleBytes(SERIAL);
    expect(WebAssembly.validate(bytes as BufferSource)).toBe(true);
    expect(inspectReactor(bytes)).toEqual({ kind: 'serial-reactor', memory: { initial: 17, maximum: 16384, shared: false } });
  });

  it('finds a threaded reactor and the maximum of the shared memory it imports', () => {
    const bytes = moduleBytes(THREADED);
    expect(inspectReactor(bytes)).toEqual({ kind: 'threaded-reactor', memory: { initial: 17, maximum: 16384, shared: true } });
  });

  it('rejects a module without every ABI v2 export', () => {
    expect(() => inspectReactor(moduleBytes({ ...SERIAL, exports: ABI_EXPORTS.filter((e) => e !== 'losat_web2_run') }))).toThrow(
      'the engine module lacks the exports losat_web2_run',
    );
  });

  it('rejects a threaded module without a shared memory with a maximum, or without wasi_thread_start', () => {
    expect(() => inspectReactor(moduleBytes({ ...THREADED, importedMemory: [2, 17] }))).toThrow(/one shared memory with a maximum/);
    expect(() => inspectReactor(moduleBytes({ ...THREADED, exports: ABI_EXPORTS }))).toThrow(/lacks wasi_thread_start/);
  });

  it('rejects bytes that are not a WebAssembly module', () => {
    expect(() => inspectReactor(new TextEncoder().encode('<!doctype html>'))).toThrow('the engine file is not a WebAssembly module');
  });
});

describe('loadReactor', () => {
  const bytes = moduleBytes(SERIAL);
  const asset = { url: '/engine/serial.wasm', sha256: createHash('sha256').update(bytes).digest('hex'), size: bytes.length, artifactSha256: 'a'.repeat(64) };
  const serve = (type: string, status = 200) =>
    vi.stubGlobal('fetch', async () => new Response(status === 200 ? (bytes as BufferSource) : null, { status, headers: { 'Content-Type': type } }));
  afterEach(() => {
    vi.unstubAllGlobals();
    vi.restoreAllMocks();
  });

  it('compiles a Wasm response while it arrives and returns the module after the SHA-256 check', async () => {
    serve('application/wasm');
    const streaming = vi.spyOn(WebAssembly, 'compileStreaming');
    const loaded = await loadReactor(asset);
    expect(streaming).toHaveBeenCalledTimes(1);
    expect(loaded.module).toBeInstanceOf(WebAssembly.Module);
    expect(loaded.artifact.kind).toBe('serial-reactor');
  });

  it('refuses bytes of another SHA-256 even when they compile', async () => {
    serve('application/wasm');
    await expect(loadReactor({ ...asset, sha256: 'b'.repeat(64) })).rejects.toThrow(
      `the engine file /engine/serial.wasm is not the one of this release (SHA-256 ${asset.sha256})`,
    );
  });

  it('compiles from the bytes when the response is not served as Wasm', async () => {
    serve('application/octet-stream');
    const streaming = vi.spyOn(WebAssembly, 'compileStreaming');
    expect((await loadReactor(asset)).artifact.kind).toBe('serial-reactor');
    expect(streaming).not.toHaveBeenCalled();
  });

  it('reports a response that is not OK', async () => {
    serve('text/plain', 404);
    await expect(loadReactor(asset)).rejects.toThrow('the engine file /engine/serial.wasm could not be loaded (HTTP 404)');
  });
});

describe('EngineRuntime.prepare', () => {
  // A serial module whose ABI version export returns 2, as a reactor's does.
  function setup(use?: InstanceUse) {
    const calls = { serial: 0, threads: 0 };
    const runtime = new EngineRuntime({
      serial: async () => {
        calls.serial++;
        return WebAssembly.compile(moduleBytes({ ...SERIAL, returns: 2 }) as BufferSource);
      },
      threads: async () => {
        calls.threads++;
        throw new Error('unused');
      },
      builds: { serial: 'serial', threads: 'threads' },
      faultChannel: 'faults',
      keeperChannel: 'keeper',
      use,
    });
    return { runtime, calls };
  }

  it('starts the serial instance for one thread, once', async () => {
    const { runtime, calls } = setup();
    await runtime.prepare(1);
    await runtime.prepare(1);
    expect(calls.serial).toBe(1);
  });

  it('keeps a working instance of the other module, which a queued search may still use', async () => {
    const { runtime, calls } = setup();
    await runtime.prepare(1);
    await runtime.prepare(4);
    expect(calls).toEqual({ serial: 1, threads: 0 });
  });

  it('warms up once per argv, ignores a warm-up that the instance cannot run, and keeps the instance', async () => {
    const { runtime, calls } = setup();
    vi.mocked(warmUp).mockClear();
    const argv = ['tblastx', '-query', 'q.fa', '-subject', 's.fa'];
    // The module's exports only return 2: the warm-up's register finds no response.
    await runtime.prepare(1, argv);
    await runtime.prepare(1, [...argv]);
    await runtime.prepare(1, [...argv, '-seg', 'no']);
    await runtime.prepare(1);
    expect(vi.mocked(warmUp).mock.calls.map(([, called]) => called)).toEqual([argv, [...argv, '-seg', 'no']]);
    expect(vi.mocked(warmUp).mock.results.map((result) => result.value)).toEqual([false, false]);
    expect(calls.serial).toBe(1);
  });

  it('starts no instance once the Gateway ends the runtime', async () => {
    const use = InstanceUse.create();
    const { runtime, calls } = setup(use);
    expect(use.retire()).toBe(false);
    // The keeper of a second threaded instance would answer on the Gateway's channel for its memory.
    await expect(runtime.prepare(4)).rejects.toThrow('the runtime is ending');
    await expect(runtime.prepare(1)).rejects.toThrow('the runtime is ending');
    expect(calls).toEqual({ serial: 0, threads: 0 });
  });

  it('records that the threaded engine cannot start, and then prepares the serial one', async () => {
    const { runtime, calls } = setup();
    // Node is not cross-origin isolated: the threaded instance cannot start.
    await expect(runtime.prepare(4)).rejects.toThrow(/cross-origin isolated/);
    await runtime.prepare(4);
    expect(calls).toEqual({ serial: 1, threads: 0 });
  });
});

describe('stopMessage', () => {
  it('reports the last line that the engine wrote, and the trap', () => {
    expect(stopMessage(new Error('unreachable'), 'thread main panicked at src/x.rs\nindex out of bounds\n')).toBe(
      'The engine stopped: index out of bounds (unreachable)',
    );
    expect(stopMessage(new Error('unreachable'), '')).toBe('The engine stopped: unreachable');
  });

  it('says when the engine ran out of memory', () => {
    expect(stopMessage(new Error('unreachable'), 'memory allocation of 1048576 bytes failed\n')).toBe(
      'The engine stopped (it ran out of memory): memory allocation of 1048576 bytes failed (unreachable)',
    );
  });
});

describe('ReactorAbi', () => {
  /** Exports whose describe and register send their JSON response in chunks of `chunk` bytes. */
  function fakeReactor(response: string, chunk: number) {
    const memory = new WebAssembly.Memory({ initial: 2 });
    let next = 1024;
    const released: number[] = [];
    const ended: number[] = [];
    const bytes = new TextEncoder().encode(response);
    const ref: { abi?: ReactorAbi } = {};
    const send = () => {
      for (let offset = 0; offset < bytes.length; offset += chunk) {
        const part = bytes.subarray(offset, offset + chunk);
        new Uint8Array(memory.buffer, 60_000, part.length).set(part);
        ref.abi!.emit(2, 60_000, part.length);
      }
    };
    const exports = {
      memory,
      losat_web2_alloc: (len: number) => {
        const ptr = next;
        next += len;
        return ptr;
      },
      losat_web2_dealloc: () => 0,
      losat_web2_describe: () => {
        send();
        return 0;
      },
      losat_web2_register: () => {
        send();
        return 7;
      },
      losat_web2_release: (handle: number) => {
        released.push(handle);
        return 0;
      },
      losat_web2_inflate_begin: () => {
        send();
        return 5;
      },
      // Inflated bytes on stream 5, and a stream that inflation never sends.
      losat_web2_inflate_chunk: () => {
        new Uint8Array(memory.buffer, 60_000, 3).set([1, 2, 3]);
        ref.abi!.emit(5, 60_000, 3);
        ref.abi!.emit(6, 60_000, 1);
        return 0;
      },
      losat_web2_inflate_end: (handle: number) => {
        ended.push(handle);
        return 0;
      },
    } as unknown as AbiExports;
    ref.abi = new ReactorAbi(exports);
    return { abi: ref.abi, released, ended };
  }

  it('joins a response that the engine sends in several chunks (1 MiB each in the adapter)', () => {
    const response = JSON.stringify({ program: 'blastn', formats: [0, 6, 7], parameters: [{ flag: '-x', help: 'h'.repeat(5000), takes_value: true }] });
    const { abi } = fakeReactor(response, 1000);
    expect(abi.describe('blastn')).toEqual(JSON.parse(response));
  });

  it('releases a registered handle whose response cannot be read', () => {
    const { abi, released } = fakeReactor('{"handle": 7, "records": [', 4);
    expect(() => abi.register('blastn', 1, new Uint8Array([62, 115]))).toThrow(SyntaxError);
    expect(released).toEqual([7]);
  });

  /**
   * Instances that answer `inflate_begin` for a `.gz` name with the gzip response, as a list: each
   * records the handles its `inflate_end` was called with, the `trapAt`-th instance opened
   * (counting from 0) traps in `inflate_chunk` and the `beginTrapAt`-th in `inflate_begin`.
   */
  function fakeInstances(trapAt?: number, beginTrapAt?: number) {
    const opened: Array<{ readonly abi: ReactorAbi; readonly ended: number[] }> = [];
    const open = async () => {
      const trap = opened.length === trapAt;
      const beginTrap = opened.length === beginTrapAt;
      const ended: number[] = [];
      const memory = new WebAssembly.Memory({ initial: 2 });
      let next = 1024;
      let handles = 0;
      const ref: { abi?: ReactorAbi } = {};
      const put = (bytes: Uint8Array, stream: number) => {
        new Uint8Array(memory.buffer, 60_000, bytes.length).set(bytes);
        ref.abi!.emit(stream, 60_000, bytes.length);
      };
      const exports = {
        memory,
        losat_web2_alloc: (len: number) => {
          const ptr = next;
          next += len;
          return ptr;
        },
        losat_web2_dealloc: () => 0,
        losat_web2_inflate_begin: (ptr: number, len: number) => {
          if (beginTrap) throw new WebAssembly.RuntimeError('unreachable');
          // The engine's name rule, for `.gz` only.
          if (!/\.gz$/i.test(new TextDecoder().decode(new Uint8Array(memory.buffer, ptr, len)))) return 0;
          put(new TextEncoder().encode('{"compression":"gzip","plain_name":"a.fa"}'), 2);
          return ++handles;
        },
        losat_web2_inflate_chunk: () => {
          if (trap) throw new WebAssembly.RuntimeError('unreachable');
          put(new TextEncoder().encode('>a\nACGT\n'), 5);
          return 0;
        },
        losat_web2_inflate_end: (handle: number) => {
          ended.push(handle);
          return 0;
        },
      } as unknown as AbiExports;
      ref.abi = new ReactorAbi(exports);
      opened.push({ abi: ref.abi, ended });
      return ref.abi;
    };
    return { opened, open };
  }

  it('opens an instance for each file, none of them the shared one, and lets it go at the end of its inflation', async () => {
    const { opened, open } = fakeInstances();
    const shared = reopening(open);
    const inflater = new ReactorInflater(shared, open);
    const a = await inflateFile(new File([new Uint8Array([0x1f, 0x8b])], 'a.fa.gz'), inflater);
    const b = await inflateFile(new File([new Uint8Array([0x1f, 0x8b])], 'b.fa.gz'), inflater);
    expect(await a.text()).toBe('>a\nACGT\n');
    expect(await b.text()).toBe('>a\nACGT\n');
    // The shared instance (0) only read the name rule, twice, and ended what it began; each file ran on a new one.
    expect(opened.map((instance) => instance.ended)).toEqual([[1, 2], [1], [1]]);
    expect(new Set(opened.map((instance) => instance.abi)).size).toBe(3);
    // A plain name opens no instance, and a failed or cancelled inflation lets its instance go as well.
    await inflateFile(new File(['>a\n'], 'plain.fa'), inflater);
    expect(opened).toHaveLength(3);
    const open1 = await inflater.open('c.fa.gz', () => undefined);
    await open1!.cancel();
    expect(opened.map((instance) => instance.ended)).toEqual([[1, 2, 3], [1], [1], [1]]);
    await expect(open1!.push(new Uint8Array([1]))).rejects.toThrow('the inflation has ended');
  });

  it('stops the inflation of an instance that traps, ends none on it, and leaves the other inflations and the shared instance as they were', async () => {
    // Instance 0 is the shared one; the first inflation's instance (1) traps in `inflate_chunk`.
    const { opened, open } = fakeInstances(1);
    const shared = reopening(open);
    const inflater = new ReactorInflater(shared, open);
    const first = await inflater.open('a.fa.gz', () => undefined);
    const second = await inflater.open('b.fa.gz', () => undefined);
    await expect(first!.push(new Uint8Array([1]))).rejects.toThrow(EngineStoppedError);
    // Its instance is dead (ABI v2: no call after a trap) and is not called again.
    await first!.cancel();
    await expect(first!.push(new Uint8Array([1]))).rejects.toThrow('the inflation has ended');
    expect(opened.map((instance) => instance.ended)).toEqual([[1, 2], [], []]);
    // The other inflation runs on its own instance and goes on.
    await second!.push(new Uint8Array([1]));
    await second!.end();
    expect(opened.map((instance) => instance.ended)).toEqual([[1, 2], [], [1]]);
    // The shared instance still works, and the next file opens a new instance.
    expect(opened[0]!.abi.stoppedBy).toBeUndefined();
    const file = await inflateFile(new File([new Uint8Array([0x1f, 0x8b])], 'c.fa.gz'), inflater);
    expect(await file.text()).toBe('>a\nACGT\n');
    expect(opened).toHaveLength(4);
  });

  it('refuses the file whose name-rule call traps the shared instance, and the next file reopens the shared instance', async () => {
    // Instance 0, the first shared one, traps in `inflate_begin`; `reopening` replaces it at its next use.
    const { opened, open } = fakeInstances(undefined, 0);
    const inflater = new ReactorInflater(reopening(open), open);
    await expect(inflater.open('a.fa.gz', () => undefined)).rejects.toThrow(EngineStoppedError);
    expect(opened).toHaveLength(1);
    expect(opened[0]!.abi.stoppedBy).toBeDefined();
    const file = await inflateFile(new File([new Uint8Array([0x1f, 0x8b])], 'b.fa.gz'), inflater);
    expect(await file.text()).toBe('>a\nACGT\n');
    // A new shared instance (1) read the name rule and ended it; the file ran on instance 2.
    expect(opened.map((instance) => instance.ended)).toEqual([[], [1], [1]]);
  });

  it('reads the inflate response, gives stream 5 of a chunk to its sink only, and ends an inflation whose response cannot be read', () => {
    const { abi, ended } = fakeReactor('{"compression":"zstd","plain_name":"a.fa"}', 7);
    expect(abi.inflateBegin('a.fa.zst')).toEqual({ handle: 5, compression: 'zstd', plainName: 'a.fa' });
    const received: number[][] = [];
    abi.inflateChunk(5, new Uint8Array([9]), (bytes) => received.push([...bytes]));
    expect(received).toEqual([[1, 2, 3]]);
    expect(ended).toEqual([]);
    const broken = fakeReactor('{"compression":', 4);
    expect(() => broken.abi.inflateBegin('a.fa.zst')).toThrow(SyntaxError);
    expect(broken.ended).toEqual([5]);
  });
});

// --- the built reactors (LOSAT_WEB_REACTORS) ---------------------------------------------------

const reactors = findReactors();

describe.skipIf(reactors === undefined)('the built reactors', () => {
  it('are the two reactor kinds, with the shared memory that the threaded build declares (plan TD-7)', () => {
    expect(inspectReactor(reactors!.serial.bytes).kind).toBe('serial-reactor');
    const threaded = inspectReactor(reactors!.threads.bytes);
    expect(threaded).toMatchObject({ kind: 'threaded-reactor', memory: { shared: true } });
    // The maximum that build_reactors.py recorded: 65,535 pages, or 16,384 in a build from
    // before the 1 GiB limit was lifted (plan DW-28).
    const recorded = JSON.parse(readFileSync(`${reactors!.dir}/artifacts.json`, 'utf8')) as {
      artifacts: { threads: { memory: { maximum: number } } };
    };
    expect(threaded.memory.maximum).toBe(recorded.artifacts.threads.memory.maximum);
    expect([16_384, 65_535]).toContain(threaded.memory.maximum);
    // The served threaded module is the artifact with the shared-memory guard.
    expect(reactors!.threads.sha256).not.toBe(reactors!.threads.artifactSha256);
    const sha256 = (path: string) => createHash('sha256').update(readFileSync(path)).digest('hex');
    expect(sha256(`${reactors!.dir}/losat-web-serial.wasm`)).toBe(reactors!.serial.sha256);
    expect(sha256(`${reactors!.dir}/losat-web-threads.wasm`)).toBe(reactors!.threads.artifactSha256);
  });
});

describe.skipIf(reactors === undefined)('responses of more than 1 MiB from the serial reactor', () => {
  // 25,000 records: the scan response (about 330 bytes a record) and the register response
  // (about 50 bytes a record) are sent in several 1 MiB chunks (ABI v2 §5).
  const fasta = new TextEncoder().encode(
    Array.from({ length: 25_000 }, (_, i) => `>record_${i} description ${i}\nACGTACGTACGTACGTACGTACGTACGTACGT\n`).join(''),
  );

  it('scans and registers thousands of records', async () => {
    const instance = await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource));
    const scanner = new ReactorScanner(async () => instance.abi);
    const scan = await scanner.scan(1, (async function* () {
      yield fasta;
    })());
    expect(scan.records).toHaveLength(25_000);
    expect(scan.records[24_999]).toMatchObject({ id: 'record_24999', length: 32 });
    const registered = instance.abi.register('blastn', 1, fasta);
    expect(registered.records).toHaveLength(25_000);
    expect(registered.records[24_999]).toEqual({ id: 'record_24999', length: 32 });
    instance.abi.release(registered.handle);
  }, 60_000);
});

describe.skipIf(reactors === undefined)('RecordScanner contract: the serial reactor in Node', () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const scanner = new ReactorScanner(reopening(async () => (await instantiateSerial(module!)).abi));
  for (const contractCase of RECORD_SCANNER_CASES) {
    it(contractCase.name, async () => {
      await contractCase.run({ scanner });
    });
  }
});

describe.skipIf(reactors === undefined)('the serial reactor scans with the NCBI reader kinds', () => {
  const scanner = new ReactorScanner(reopening(async () => (await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource))).abi));
  const bytes = () => (async function* () {
    yield new TextEncoder().encode('>a desc\nACGT\n>b\nACGTAC\n');
  })();

  // Kinds 1 (nucleotide flags) and 2 (protein flags) are real since S9 (abi_v2.md §4, §9).
  it.each([1, 2])('scans a small input with parser kind %i', async (kind) => {
    const scan = await scanner.scan(kind as FastaParserKind, bytes());
    expect(scan.records.map(({ id, length }) => ({ id, length }))).toEqual([
      { id: 'a', length: 4 },
      { id: 'b', length: 6 },
    ]);
  });
});

describe.skipIf(reactors === undefined)('the serial reactor answers the search form', () => {
  const open = async () => (await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource))).abi;

  // FakeEngine's describe.json is a copy of the reactor's describe. After an engine change
  // to an option, its help or its default, `vitest run tests/unit/engine-runtime.test.ts -u`
  // (with LOSAT_WEB_REACTORS) writes it again.
  it("FakeEngine's describe.json is the reactor's describe of every program", async () => {
    const abi = await open();
    const described = Object.fromEntries(PROGRAMS.map(({ id }) => [id, abi.describe(id)]));
    await expect(`${JSON.stringify(described, null, 1)}\n`).toMatchFileSnapshot('../../src/infra/fake/describe.json');
  });

  it('every field of the search form is an option that the engine describes', async () => {
    const abi = await open();
    for (const program of PROGRAMS) {
      const described = new Set(toProgramDescription(abi.describe(program.id)).parameters.map((option) => option.flag));
      for (const section of program.sections) {
        for (const field of section.fields) expect(described, `${program.id} ${field.flag}`).toContain(field.flag);
      }
    }
  });

  it("checks an input with the engine's register: the verdict and the message", async () => {
    const checker = new ReactorInputChecker(reopening(open));
    expect(await checker.check('blastn', 'query', new TextEncoder().encode('>a\nACGT\n>b\nACGT\n'))).toEqual({
      ok: true,
      records: [
        { id: 'a', length: 4 },
        { id: 'b', length: 4 },
      ],
    });
    // The engine's register reads with NCBI's reader (abi_v2.md §4), which accepts what the
    // bio reader of ABI v1 refused: its records and lengths are the reader's (it stores
    // the four residues of ACLGT in a nucleotide query and of MK-LV in a protein subject).
    expect(await checker.check('blastn', 'query', new TextEncoder().encode('>a\nACGT\n>b\nACLGT\n'))).toEqual({
      ok: true,
      records: [
        { id: 'a', length: 4 },
        { id: 'b', length: 4 },
      ],
    });
    expect(await checker.check('blastp', 'subject', new TextEncoder().encode('>p\nMK-LV\n'))).toEqual({
      ok: true,
      records: [{ id: 'p', length: 4 }],
    });
    // What register still refuses: a first line that NCBI may read as a Seq-id and a gap line.
    const seqId = await checker.check('blastn', 'query', new TextEncoder().encode('AB123456\nACGT\n'));
    expect(seqId.ok).toBe(false);
    expect(!seqId.ok && seqId.message).toMatch(/not supported in the browser/);
    const gap = await checker.check('blastp', 'subject', new TextEncoder().encode('>p\nMK\n>?10\nMK\n'));
    expect(gap.ok).toBe(false);
    expect(!gap.ok && gap.message).toMatch(/not supported in the browser/);
  });
});

describe.skipIf(reactors === undefined)('the warm-up on the serial reactor', () => {
  const open = async () => (await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource))).abi;
  // Streams 0, 6 and 7 (the three output formats) of one search on an instance.
  const search = (abi: ReactorAbi, argv: readonly string[], query: Uint8Array, subject: Uint8Array) => {
    const q = abi.register(argv[0]!, 0, query).handle;
    const s = abi.register(argv[0]!, 1, subject).handle;
    const out: Record<number, string> = { 0: '', 6: '', 7: '' };
    const decoder = new TextDecoder();
    abi.run(argv, q, s, (stream, bytes) => {
      if (stream in out) out[stream] += decoder.decode(bytes);
    });
    abi.release(q);
    abi.release(s);
    return out;
  };

  for (const program of PROGRAMS) {
    it(`${program.id}: runs on the built-in sequences, which have hits, and leaves the next search as on a new instance`, async () => {
      const argv = [program.id, '-query', 'query.fa', '-subject', 'subject.fa'];
      const inputs = warmUpInputs(program.id);
      const fresh = search(await open(), [...argv, '-num_threads', '1'], inputs.query, inputs.subject);
      // A hit for each of the three blocks.
      expect(fresh[6]!.trim().split('\n').length).toBeGreaterThanOrEqual(3);
      const abi = await open();
      // A region of the user's input is left out of the warm-up's argv.
      expect(warmUp(abi, [...argv, '-query_loc', '5000-6000'])).toBe(true);
      expect(search(abi, [...argv, '-num_threads', '1'], inputs.query, inputs.subject)).toEqual(fresh);
    });
  }

  it('reports an argv that the engine refuses, and the instance still runs searches', async () => {
    const abi = await open();
    expect(warmUp(abi, ['blastn', '-query', 'q.fa', '-subject', 's.fa', '-word_size', '-1'])).toBe(false);
    expect(abi.stoppedBy).toBeUndefined();
    expect(warmUp(abi, ['blastn', '-query', 'q.fa', '-subject', 's.fa'])).toBe(true);
  });
});

describe.skipIf(reactors === undefined)("the serial reactor inflates sources with the command line's decoders", () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const instantiate = async () => (await instantiateSerial(module!)).abi;
  const inflater = new ReactorInflater(reopening(instantiate), instantiate);
  const decoder = new TextDecoder();
  const A = '>a one\nACGT\n';
  const B = '>b two\nGGCC\nTT\n';
  // pzstd output: a skippable frame, then one frame of `>a\nACGT…` and `>b\nTTTT…`.
  const PZSTD =
    'UCpNGAQAAAA/AAAAKLUv/QRYlQEAhAI+YQpBQ0dUQUFBQ0NDR0dHVFRUCj5iClRUVFRHR0dHQ0NDQ0FBQUEKAgCYibw/ZwTZFtvJ';
  const AB_TEXT = '>a\nACGTACGTACGTACGTACGTACGTAAAACCCGGGTTT\n>b\nTTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n';
  // `bzip2 -9c` of A, then of B: two bzip2 streams one after the other (as pbzip2 or `cat` writes them).
  const TWO_BZIP2 =
    'QlpoOTFBWSZTWa6jEEsAAAHfAAAQQAAAASiABAAiAaAAMQDTTQQBppoCaYaTyXi7kinChIV1GIJYQlpoOTFBWSZTWf/ayGYAAAJfgAAQQAAAAQiABAAQAISAIAAiDRppoQAwKKwKAaT5+LuSKcKEh/7WQzA=';
  const sentence = (name: string, suffix: string, format: string, fault: string) =>
    `'${name}' is named ${suffix} but is not ${format} data (${fault}); LOSAT reads a file named ${suffix} through ${format}`;
  const inflated = async (bytes: Uint8Array | string, name: string) => {
    const file = await inflateFile(new File([bytes as BlobPart], name), inflater);
    return { name: file.name, text: decoder.decode(await file.arrayBuffer()) };
  };
  const failure = async (bytes: Uint8Array | string, name: string) => {
    try {
      await inflateFile(new File([bytes as BlobPart], name), inflater);
    } catch (error) {
      return (error as Error).message;
    }
    throw new Error(`${name} was inflated`);
  };

  it('keeps a plain name as it is, and inflates every member and frame under the plain name', async () => {
    expect(await inflated(A, 'a.fa')).toEqual({ name: 'a.fa', text: A });
    expect(await inflated(Buffer.concat([gzipSync(A), gzipSync(B)]), 'two.fa.GZ')).toEqual({ name: 'two.fa', text: A + B });
    expect(await inflated(Buffer.concat([zstdCompressSync(A), zstdCompressSync(B)]), 'two.fa.zst')).toEqual({ name: 'two.fa', text: A + B });
    expect(await inflated(Buffer.from(PZSTD, 'base64'), 'p.fa.zst')).toEqual({ name: 'p.fa', text: AB_TEXT });
    expect(await inflated(Buffer.from(TWO_BZIP2, 'base64'), 'two.fa.Bz2')).toEqual({ name: 'two.fa', text: A + B });
    expect(await inflated(gzipSync(''), '.gz')).toEqual({ name: 'input', text: '' });
  });

  it('inflates a file of several 8 MiB slices into several 8 MiB parts', async () => {
    // 17 MiB of bytes that do not compress (xorshift32): each compressed file holds three slices.
    let state = 2463534242;
    const bytes = new Uint8Array(17 * 1024 * 1024 + 5).map(() => {
      state ^= state << 13;
      state ^= state >>> 17;
      state ^= state << 5;
      return state & 0xff;
    });
    const MIB = 1024 * 1024;
    for (const [suffix, pack] of [
      ['gz', gzipSync],
      ['zst', zstdCompressSync],
    ] as const) {
      const packed = pack(bytes);
      expect(packed.length).toBeGreaterThan(16 * MIB);
      // The sizes of the Blobs made while the file inflates: its parts.
      const parts: number[] = [];
      const RealBlob = globalThis.Blob;
      vi.stubGlobal(
        'Blob',
        class extends RealBlob {
          constructor(blobParts?: BlobPart[], options?: BlobPropertyBag) {
            super(blobParts, options);
            parts.push(this.size);
          }
        },
      );
      let file: File;
      try {
        file = await inflateFile(new File([packed as BlobPart], `big.fa.${suffix}`), inflater);
      } finally {
        vi.unstubAllGlobals();
      }
      expect(file.name).toBe('big.fa');
      expect(Buffer.compare(new Uint8Array(await file.arrayBuffer()), bytes)).toBe(0);
      // Parts of 8 MiB and at most one engine chunk more (emit.rs CHUNK, 1 MiB), the last one with the
      // rest: two or three parts for these 17 MiB, as the chunks fall.
      expect(parts.reduce((sum, size) => sum + size, 0)).toBe(bytes.length);
      expect(parts.length).toBeGreaterThanOrEqual(2);
      for (const size of parts.slice(0, -1)) {
        expect(size).toBeGreaterThanOrEqual(8 * MIB);
        expect(size).toBeLessThanOrEqual(9 * MIB);
      }
      expect(parts.at(-1)).toBeGreaterThan(0);
      expect(parts.at(-1)).toBeLessThanOrEqual(9 * MIB);
    }
  }, 60_000);

  it("refuses a damaged, cut, mislabelled or empty file with the command line's sentence, and goes on", async () => {
    const text = Array.from({ length: 300 }, (_, i) => `>r${i}\n${'ACGTTGCAAGTC'.repeat(30)}\n`).join('');
    const checked = Buffer.from(zstdCompressSync(text, { params: { [zlib.ZSTD_c_checksumFlag]: 1 } }));
    checked[checked.length - 1]! ^= 1;
    expect(await failure(checked, 'bad.fa.zst')).toBe(sentence('bad.fa.zst', '.zst', 'zstd', "a zstd frame's checksum does not match its data"));
    const gz = gzipSync(text);
    expect(await failure(gz.subarray(0, gz.length >> 1), 'cut.fa.gz')).toBe(sentence('cut.fa.gz', '.gz', 'gzip', 'incomplete deflate stream'));
    expect(await failure(A, 'text.fa.gz')).toBe(sentence('text.fa.gz', '.gz', 'gzip', 'invalid gzip header'));
    expect(await failure('', 'e.fa.gz')).toBe(sentence('e.fa.gz', '.gz', 'gzip', 'the file is empty'));
    expect(await failure('', 'e.fa.zst')).toBe(sentence('e.fa.zst', '.zst', 'zstd', 'the file is empty'));
    expect(await failure(A, 'a.fa.bz2')).toBe(sentence('a.fa.bz2', '.bz2', 'bzip2', 'bzip2: bz2 header missing'));
    expect(await inflated(gzipSync(A), 'after.fa.gz')).toEqual({ name: 'after.fa', text: A });
  });

  it('stops a downloaded Example set past its limit of inflated bytes, and goes on', async () => {
    const big = new File([gzipSync('>r\n' + 'ACGT'.repeat(1 << 20))], 'set.fna.gz');
    await expect(inflateFile(big, inflater, 1_000_000)).rejects.toMatchObject({ name: INFLATED_TOO_LARGE });
    expect((await inflateFile(big, inflater, 5 * 1024 * 1024)).size).toBe(3 + 4 * (1 << 20));
    expect(await inflated(gzipSync(A), 'after.fa.gz')).toEqual({ name: 'after.fa', text: A });
  });

  it("grows the linear memory of an instance by less than 1.5 windows when a zstd file's frame declares a window W", async () => {
    // A frame of window 32 MiB (windowLog 25, header Window_Descriptor checked below) over 40 MiB of
    // bytes that do not compress. The decoder's ring buffer is reserved at the frame's window
    // (ruzstd `DecodeBuffer::reset`, at the first frame as at later ones: LOSAT/src/blastinput/compressed.rs
    // `ZstdFrames::step`); where it grows by doubling instead, each step allocates the larger buffer
    // while the smaller lives, a Wasm allocator cannot reuse the smaller regions for it, and the
    // memory ends near 2 W, which no instance gives back. The reserved buffer is W and a little over; the
    // 8 MiB slice of the file that inflateFile hands the engine at once adds to it.
    const MIB = 1024 * 1024;
    const W = 32 * MIB;
    let state = 2463534242;
    const bytes = new Uint8Array(40 * MIB).map(() => {
      state ^= state << 13;
      state ^= state >>> 17;
      state ^= state << 5;
      return state & 0xff;
    });
    const packed = zstdCompressSync(bytes, { params: { [zlib.ZSTD_c_windowLog]: 25 } });
    // Window_Descriptor (byte 5 without Single_Segment_Flag): exponent 15 above 2^10, no mantissa.
    expect(packed[4]! & 0x20).toBe(0);
    expect(1 << (10 + (packed[5]! >> 3))).toBe(W);
    const instances: ReactorAbi[] = [];
    const instantiate = async () => {
      const { abi } = await instantiateSerial(module!);
      instances.push(abi);
      return abi;
    };
    const own = new ReactorInflater(reopening(instantiate), instantiate);
    const file = await inflateFile(new File([packed as BlobPart], 'w.fa.zst'), own);
    expect(file.size).toBe(bytes.length);
    // Instance 0 is the shared one (the name rule); instance 1 inflated the file.
    expect(instances).toHaveLength(2);
    const grown = instances[1]!.memoryBytes() - instances[0]!.memoryBytes();
    expect(grown).toBeGreaterThan(W);
    expect(grown).toBeLessThan(1.5 * W);
  }, 60_000);
});
