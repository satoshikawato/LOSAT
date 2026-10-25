// Pure parts of the engine runtime: the thread and renewal policy, the check of an engine
// module before it is compiled, and the message of an instance that stopped. With
// LOSAT_WEB_REACTORS, the module check also runs on the built reactors, and the
// RecordScanner contract runs against the serial reactor in Node.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { findReactors } from '../../build/reactors';
import { ReactorAbi, stopMessage, type AbiExports } from '../../src/infra/reactor/abi';
import { ABI_EXPORTS, inspectReactor } from '../../src/infra/reactor/artifact';
import { instantiateSerial, loadReactor } from '../../src/infra/reactor/instance';
import { ReactorInputChecker } from '../../src/infra/reactor/checker';
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
  function setup() {
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
    } as unknown as AbiExports;
    ref.abi = new ReactorAbi(exports);
    return { abi: ref.abi, released };
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
    const scan = await scanner.scan(0, (async function* () {
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
    const described = Object.fromEntries(['blastn', 'blastp', 'tblastn', 'tblastx'].map((program) => [program, abi.describe(program)]));
    await expect(`${JSON.stringify(described, null, 1)}\n`).toMatchFileSnapshot('../../src/infra/fake/describe.json');
  });

  it('every field of the search form is an option that the engine describes', async () => {
    const abi = await open();
    for (const program of PROGRAMS.filter((p) => p.unavailable === undefined)) {
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
    expect(!seqId.ok && seqId.message).toMatch(/not supported by LOSAT Web/);
    const gap = await checker.check('blastp', 'subject', new TextEncoder().encode('>p\nMK\n>?10\nMK\n'));
    expect(gap.ok).toBe(false);
    expect(!gap.ok && gap.message).toMatch(/not supported by LOSAT Web/);
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

  for (const program of PROGRAMS.filter((p) => p.unavailable === undefined)) {
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
