// Pure parts of the engine runtime: the thread and renewal policy, the check of an engine
// module before it is compiled, and the message of an instance that stopped. With
// LOSAT_WEB_REACTORS, the module check also runs on the built reactors, and the
// RecordScanner contract runs against the serial reactor in Node.
import { createHash } from 'node:crypto';
import { readFileSync } from 'node:fs';
import { describe, expect, it } from 'vitest';
import { findReactors } from '../../build/reactors';
import { ReactorAbi, stopMessage, type AbiExports } from '../../src/infra/reactor/abi';
import { ABI_EXPORTS, inspectReactor } from '../../src/infra/reactor/artifact';
import { instantiateSerial } from '../../src/infra/reactor/instance';
import { ReactorInputChecker } from '../../src/infra/reactor/checker';
import { toProgramDescription } from '../../src/infra/reactor/control';
import { ReactorScanner, reopening } from '../../src/infra/reactor/scanner';
import { FakeEngine } from '../../src/infra/fake/fake-engine';
import { PROGRAMS } from '../../src/domain/programs';
import { AUTO_MAX_THREADS, AUTO_SERIAL_BELOW_BYTES, chooseThreads, DEFAULT_RENEWAL, renewalReason } from '../../src/infra/engine-worker/policy';
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
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 8)).toBe(Math.min(AUTO_MAX_THREADS, 4));
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 64)).toBe(AUTO_MAX_THREADS);
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 1)).toBe(1);
    // A browser that does not report the processors.
    expect(chooseThreads('auto', AUTO_SERIAL_BELOW_BYTES, 0)).toBe(1);
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
}

/** A module whose functions all have type () -> i32 and return 0. */
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
    ...section(10, vector(functions.map(() => [4, 0, 0x41, 0, 0x0b]))),
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
  it('are the two reactor kinds, with the shared memory of the certified threaded build (plan TD-7)', () => {
    expect(inspectReactor(reactors!.serial.bytes).kind).toBe('serial-reactor');
    expect(inspectReactor(reactors!.threads.bytes)).toMatchObject({ kind: 'threaded-reactor', memory: { maximum: 16384, shared: true } });
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

describe.skipIf(reactors === undefined)('the serial reactor answers the search form', () => {
  const open = async () => (await instantiateSerial(new WebAssembly.Module(reactors!.serial.bytes as BufferSource))).abi;

  it("FakeEngine's describe.json is the reactor's describe of every program", async () => {
    const abi = await open();
    const fake = new FakeEngine();
    for (const program of ['blastn', 'blastp', 'tblastn', 'tblastx'] as const) {
      expect(await fake.describe(program)).toEqual(toProgramDescription(abi.describe(program)));
    }
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
    const refused = await checker.check('blastn', 'query', new TextEncoder().encode('>a\nACGT\n>b\nACLGT\n'));
    expect(refused.ok).toBe(false);
    expect(!refused.ok && refused.message).toMatch(/^query record 2 \(b\) has 'L' at residue 3, .*not supported by LOSAT's BLASTN/);
    const protein = await checker.check('blastp', 'subject', new TextEncoder().encode('>p\nMK-LV\n'));
    expect(!protein.ok && protein.message).toContain("not supported by LOSAT's BLASTP");
  });
});
