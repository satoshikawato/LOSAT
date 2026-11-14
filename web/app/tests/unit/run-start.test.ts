// NCBI's application start in LOSAT Web (line RG item E1w-3b; docs/web/abi_v2.md §4 `losat_web2_start`,
// file kind 2, §5 stream 8, §9): the reader of the engine's archive of a run's files, the ABI binding
// of `start` and of a registry file held in memory, `validate` with the registry file of -conffile on
// the Data worker's instance, the Engine worker's runtime (the registry file for one search, the
// files of stream 8 apart from the outputs, with the result or the error), the gateway that hands
// them on, and the names of their downloads. The coordinator, the draft and the results have their
// cases in coordinator.test.ts, draft.test.ts and results.test.ts. With LOSAT_WEB_REACTORS, the same
// calls go to the built serial reactor.
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';
import { describe, expect, it } from 'vitest';
import { instantiateSerial, reactorControl, WasmEngine } from 'losat';
import { readStoredZip, type RunFile } from 'losat/domain';
import { runFilesOf, withRunFiles, type EngineRunRequest } from 'losat/ports';
import { findReactors } from '../../build/reactors';
import { outputFileName, runFileName } from '../../src/domain/output-format';
import { EngineRuntime } from '../../../packages/losat/src/infra/engine-worker/runtime';
import type { EngineCommand, EngineEvent, WorkerRun } from '../../../packages/losat/src/infra/engine-worker/protocol';
import { ABI_EXPORTS, inspectReactor } from '../../../packages/losat/src/infra/reactor/artifact';
import {
  EngineCallError,
  EngineStoppedError,
  FILE_CHUNK_BYTES,
  FILE_REGISTRY,
  loadFile,
  ReactorAbi,
  ROLE_QUERY,
  ROLE_SUBJECT,
  RUN_FILES_STREAM,
  type AbiExports,
  type FileLoader,
} from '../../../packages/losat/src/infra/reactor/abi';

const encode = (text: string) => new TextEncoder().encode(text);

// --- the archive of a run's files ----------------------------------------------------------------

/** The CRC-32 of ZIP, written here apart from the reader's. */
function crc(bytes: Uint8Array): number {
  let value = ~0;
  for (const byte of bytes) {
    value ^= byte;
    for (let k = 0; k < 8; k++) value = value & 1 ? 0xedb88320 ^ (value >>> 1) : value >>> 1;
  }
  return ~value >>> 0;
}

/** An archive laid out as the adapter writes it (web/adapter/src/zip.rs): stored entries, the central directory, the end record. */
function storedZip(entries: ReadonlyArray<readonly [string, Uint8Array]>): Uint8Array {
  const parts: number[] = [];
  const u16 = (value: number) => parts.push(value & 0xff, (value >>> 8) & 0xff);
  const u32 = (value: number) => parts.push(value & 0xff, (value >>> 8) & 0xff, (value >>> 16) & 0xff, (value >>> 24) & 0xff);
  const central: Array<{ name: Uint8Array; crc: number; size: number; offset: number }> = [];
  for (const [text, bytes] of entries) {
    const name = encode(text);
    const offset = parts.length;
    u32(0x04034b50);
    [20, 0, 0, 0, 33].forEach(u16);
    u32(crc(bytes));
    u32(bytes.length);
    u32(bytes.length);
    u16(name.length);
    u16(0);
    parts.push(...name, ...bytes);
    central.push({ name, crc: crc(bytes), size: bytes.length, offset });
  }
  const start = parts.length;
  for (const entry of central) {
    u32(0x02014b50);
    [20, 20, 0, 0, 0, 33].forEach(u16);
    u32(entry.crc);
    u32(entry.size);
    u32(entry.size);
    [entry.name.length, 0, 0, 0, 0].forEach(u16);
    u32(0);
    u32(entry.offset);
    parts.push(...entry.name);
  }
  const size = parts.length - start;
  u32(0x06054b50);
  [0, 0, central.length, central.length].forEach(u16);
  u32(size);
  u32(start);
  u16(0);
  return new Uint8Array(parts);
}

describe("the reader of the engine's archive of a run's files", () => {
  it('reads the stored entries in order: their UTF-8 names and their bytes, an empty file too', () => {
    const files = readStoredZip(storedZip([['run.log', encode('Warning: one\n')], ['ü.log', new Uint8Array()]]));
    expect(files).toEqual([
      { name: 'run.log', bytes: encode('Warning: one\n') },
      { name: 'ü.log', bytes: new Uint8Array() },
    ]);
    expect(readStoredZip(storedZip([]))).toEqual([]);
  });

  it('refuses what the adapter does not write, naming the engine: a wrong CRC-32, a method, a data descriptor, a cut, a wrong directory, bytes after the end', () => {
    const good = storedZip([['run.log', encode('Warning: one\n')]]);
    const changed = (at: number, value: number) => {
      const bytes = good.slice();
      bytes[at] = value;
      return bytes;
    };
    const cases: Array<readonly [Uint8Array, RegExp]> = [
      // The first byte of the file's text.
      [changed(30 + 7, 0x58), /the bytes of entry 1 do not have its CRC-32/],
      // The method: deflate.
      [changed(8, 8), /entry 1 is not stored \(method 8\)/],
      // The flags: a data descriptor.
      [changed(6, 8), /entry 1 has the flags 0x8/],
      [good.subarray(0, good.length - 5), /it ends at byte|no end record/],
      [good.subarray(0, 40), /entry 1 ends past the archive/],
      // The central header's CRC-32.
      [changed(30 + 7 + 13 + 16, 0), /the central directory differs from entry 1/],
      [new Uint8Array([...good, 0]), /bytes follow the end record/],
      [encode('not an archive'), /it ends at byte|no end record/],
    ];
    for (const [archive, reason] of cases) {
      expect(() => readStoredZip(archive)).toThrow(/^the engine's archive of the run's files \(ABI v2 stream 8\) cannot be read: /);
      expect(() => readStoredZip(archive)).toThrow(reason);
    }
  });
});

// --- the ABI binding -------------------------------------------------------------------------------

/** Exports of a fake engine for `start`, `validate` and the file calls, which log what they are given. */
function fakeExports(options: { start?: string | 'fail'; validate?: string; trapValidate?: boolean; failEnd?: boolean } = {}) {
  const memory = new WebAssembly.Memory({ initial: 4 });
  const log: string[] = [];
  let next = 4096;
  const ref: { abi?: ReactorAbi } = {};
  let error = '';
  const text = (ptr: number, len: number) => new TextDecoder().decode(new Uint8Array(memory.buffer, ptr, len));
  const respond = (json: string) => {
    const bytes = encode(json);
    new Uint8Array(memory.buffer, 60_000, bytes.length).set(bytes);
    ref.abi!.emit(2, 60_000, bytes.length);
  };
  const exports = {
    memory,
    losat_web2_alloc: (len: number) => {
      const ptr = next;
      next += len;
      return ptr;
    },
    losat_web2_dealloc: () => 0,
    losat_web2_last_error_ptr: () => {
      const bytes = encode(error);
      new Uint8Array(memory.buffer, 50_000, bytes.length).set(bytes);
      return 50_000;
    },
    losat_web2_last_error_len: () => encode(error).length,
    losat_web2_start: (ptr: number, len: number) => {
      log.push(`start ${text(ptr, len).split('\0').join(' ')}`);
      if (options.start === 'fail') {
        error = "error: unexpected argument '-version' found";
        return -1;
      }
      respond(options.start ?? '{"conffile":"c.ini","log_file":"run.log","dry_run":true}');
      return 0;
    },
    losat_web2_validate: (ptr: number, len: number) => {
      log.push(`validate ${text(ptr, len).split('\0').join(' ')}`);
      if (options.trapValidate) throw new WebAssembly.RuntimeError('unreachable');
      if (options.validate === undefined) return 0;
      error = options.validate;
      return -1;
    },
    losat_web2_file_begin: (kind: number, ptr: number, len: number, size: number) => {
      log.push(`begin ${kind} ${text(ptr, len)} ${size}`);
      return 9;
    },
    losat_web2_file_chunk: (file: number, ptr: number, len: number) => {
      log.push(`chunk ${file} ${text(ptr, len)}`);
      return 0;
    },
    losat_web2_file_end: (file: number) => {
      log.push(`end ${file}`);
      if (options.failEnd) {
        error = 'the input file ended after 1 of its 3 bytes';
        return -1;
      }
      respond('{"handle":9,"name":"c.ini","kind":"registry","size":8,"sha256":"ab"}');
      return 0;
    },
    losat_web2_file_release: (file: number) => {
      log.push(`release ${file}`);
      return 0;
    },
  } as unknown as AbiExports;
  ref.abi = new ReactorAbi(exports);
  return { abi: ref.abi, log };
}

describe('the ABI binding of the start and of a registry file', () => {
  it('reads the start JSON of an argv, or the CLI message of one that the engine cannot start', () => {
    expect(fakeExports().abi.start(['tblastx', '-dryrun'])).toEqual({ ok: true, start: { conffile: 'c.ini', logFile: 'run.log', dryRun: true } });
    expect(fakeExports({ start: '{"conffile":null,"log_file":null,"dry_run":false}' }).abi.start(['blastn'])).toEqual({
      ok: true,
      start: { conffile: null, logFile: null, dryRun: false },
    });
    expect(fakeExports({ start: 'fail' }).abi.start(['blastn', '-version'])).toEqual({ ok: false, message: "error: unexpected argument '-version' found" });
  });

  it('loads a file held in memory as kind 2: begin with its size, its bytes, end; a file that the engine refuses is released', async () => {
    const bytes = encode('[BLAST]\n');
    const { abi, log } = fakeExports();
    expect(await loadFile(abi, FILE_REGISTRY, 'c.ini', bytes.length, (start, end) => bytes.subarray(start, end))).toEqual({
      handle: 9,
      name: 'c.ini',
      kind: 'registry',
      size: 8,
      sha256: 'ab',
    });
    expect(log).toEqual(['begin 2 c.ini 8', 'chunk 9 [BLAST]\n', 'end 9']);
    const failing = fakeExports({ failEnd: true });
    await expect(loadFile(failing.abi, FILE_REGISTRY, 'c.ini', 3, () => encode('[x\n'))).rejects.toThrow(EngineCallError);
    expect(failing.log).toEqual(['begin 2 c.ini 3', 'chunk 9 [x\n', 'end 9', 'release 9']);
  });

  // One loader for both kinds of file (the WindowMasker counts file's runtime cases are window-masker.test.ts's).
  it('gives a file larger than one chunk in chunks of FILE_CHUNK_BYTES through one buffer, each call through `call`', async () => {
    const calls: string[] = [];
    const loader = {
      stoppedBy: undefined,
      fileBegin: (kind: number, name: string, size: number) => (calls.push(`begin ${kind} ${name} ${size}`), 7),
      chunkBuffer: (capacity: number) => (calls.push(`buffer ${capacity}`), { ptr: 1, capacity }),
      fileChunk: (file: number, bytes: Uint8Array) => void calls.push(`chunk ${file} ${bytes.length} ${bytes[0]}`),
      fileEnd: (file: number) => (calls.push(`end ${file}`), { handle: file, name: 'c.ini', kind: 'registry' as const, size: 0, sha256: '' }),
      fileRelease: (file: number) => void calls.push(`release ${file}`),
      freeChunkBuffer: () => void calls.push('free buffer'),
    } as unknown as FileLoader;
    const bytes = new Uint8Array(FILE_CHUNK_BYTES + 3).fill(1);
    bytes[FILE_CHUNK_BYTES] = 2;
    let wrapped = 0;
    await loadFile(loader, FILE_REGISTRY, 'c.ini', bytes.length, (start, end) => bytes.subarray(start, end), (body) => (wrapped++, body()));
    expect(calls).toEqual([
      `begin 2 c.ini ${FILE_CHUNK_BYTES + 3}`,
      `buffer ${FILE_CHUNK_BYTES}`,
      `chunk 7 ${FILE_CHUNK_BYTES} 1`,
      'chunk 7 3 2',
      'end 7',
      'free buffer',
    ]);
    expect(wrapped).toBe(6);
  });

  it('requires losat_web2_start of an engine module', () => {
    expect(ABI_EXPORTS).toContain('losat_web2_start');
    // A module of every export but start (the module check's encoder is engine-runtime.test.ts's; a
    // module of exports only is enough for the check of the names).
    const names = ABI_EXPORTS.filter((name) => name !== 'losat_web2_start');
    expect(() => inspectReactor(exportsOnlyModule(names))).toThrow('the engine module lacks the exports losat_web2_start');
  });
});

/** A module with an own memory and a function `() -> i32` exported under each other name. */
function exportsOnlyModule(names: readonly string[]): Uint8Array {
  const uleb = (value: number): number[] => {
    const out: number[] = [];
    do {
      const byte = value & 0x7f;
      value = Math.floor(value / 128);
      out.push(byte | (value ? 0x80 : 0));
    } while (value);
    return out;
  };
  const name = (text: string) => [...uleb(encode(text).length), ...encode(text)];
  const vector = (items: number[][]) => [...uleb(items.length), ...items.flat()];
  const section = (id: number, bytes: number[]) => [id, ...uleb(bytes.length), ...bytes];
  const functions = names.filter((n) => n !== 'memory');
  return new Uint8Array([
    0, 0x61, 0x73, 0x6d, 1, 0, 0, 0,
    ...section(1, [1, 0x60, 0, 1, 0x7f]),
    ...section(3, vector(functions.map(() => [0]))),
    ...section(5, [1, 1, 17, ...uleb(16384)]),
    ...section(7, vector(names.map((n) => (n === 'memory' ? [...name(n), 2, 0] : [...name(n), 0, ...uleb(functions.indexOf(n))])))),
    ...section(10, vector(functions.map(() => [4, 0, 0x41, 0, 0x0b]))),
  ]);
}

describe("the Data worker's validate with the registry file of -conffile", () => {
  it('registers the file under the given name for that validation only, and releases it whatever the verdict', async () => {
    const { abi, log } = fakeExports({ validate: 'BLAST query/options error: bad' });
    const control = reactorControl(async () => abi);
    expect(await control.validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', bytes: encode('[BLAST]\n') })).toEqual({
      ok: false,
      message: 'BLAST query/options error: bad',
    });
    expect(log).toEqual(['begin 2 c.ini 8', 'chunk 9 [BLAST]\n', 'end 9', 'validate blastp -conffile c.ini', 'release 9']);
    log.length = 0;
    // Without a file nothing is registered; start is the engine's.
    await control.validate(['blastp']);
    expect(log).toEqual(['validate blastp']);
    expect(await control.start(['blastp', '-dryrun'])).toEqual({ ok: true, start: { conffile: 'c.ini', logFile: 'run.log', dryRun: true } });
  });

  it('does not validate with a file that the engine refuses, and calls an instance that stopped no more', async () => {
    const refused = fakeExports({ failEnd: true });
    await expect(reactorControl(async () => refused.abi).validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', bytes: encode('[x\n') })).rejects.toThrow(
      'the input file ended after 1 of its 3 bytes',
    );
    expect(refused.log).toEqual(['begin 2 c.ini 3', 'chunk 9 [x\n', 'end 9', 'release 9']);
    // A trap in the validation: the stopped instance's file is not released (ABI v2 §4, "After a trap").
    const trapped = fakeExports({ trapValidate: true });
    await expect(reactorControl(async () => trapped.abi).validate(['blastp', '-conffile', 'c.ini'], { name: 'c.ini', bytes: encode('[BLAST]\n') })).rejects.toThrow(
      EngineStoppedError,
    );
    expect(trapped.log).toEqual(['begin 2 c.ini 8', 'chunk 9 [BLAST]\n', 'end 9', 'validate blastp -conffile c.ini']);
  });
});

// --- the Engine worker's runtime -------------------------------------------------------------------

/**
 * A binding with no engine behind it that records the calls of a search; its run sends outfmt 6 and
 * the archive of `runFiles` on stream 8 (in two chunks), and then fails with `fail`, or stops.
 * `heldName` is the name of a file that the instance holds (a WindowMasker counts file): `file_begin`
 * refuses it, as the adapter does.
 */
function fakeAbi(options: { runFiles?: Uint8Array; fail?: string; stop?: boolean; heldName?: string } = {}) {
  const calls: string[] = [];
  let handle = 100;
  let stopped: EngineStoppedError | undefined;
  const abi = {
    get stoppedBy() {
      return stopped;
    },
    fileBegin: (kind: number, name: string, size: number) => {
      calls.push(`begin ${kind} ${name} ${size}`);
      if (name === options.heldName) throw new EngineCallError(`a file named ${name} is registered: release it before registering another`);
      return ++handle;
    },
    chunkBuffer: (capacity: number) => (calls.push(`buffer ${capacity}`), { ptr: 1, capacity }),
    fileChunk: (file: number, bytes: Uint8Array) => void calls.push(`chunk ${file} ${new TextDecoder().decode(bytes)}`),
    fileEnd: (file: number) => (calls.push(`end ${file}`), { handle: file, name: 'c.ini', kind: 'registry', size: 8, sha256: '0' }),
    freeChunkBuffer: () => void calls.push('free buffer'),
    validate: (argv: readonly string[]) => {
      calls.push(`validate ${argv.join(' ')}`);
      return undefined;
    },
    memoryBytes: () => 1,
    register: (_program: string, role: number, bytes: Uint8Array) => {
      calls.push(`register ${role === ROLE_QUERY ? 'query' : 'subject'}`);
      return { handle: ++handle, records: [{ id: role === ROLE_SUBJECT ? 's' : 'q', length: bytes.length }] };
    },
    release: (h: number) => void calls.push(`release ${h}`),
    fileRelease: (file: number) => void calls.push(`file release ${file}`),
    run: (_argv: readonly string[], _query: number, _subject: number, sink: (stream: number, bytes: Uint8Array) => void) => {
      calls.push('run');
      sink(6, encode('q\ts\n'));
      if (options.runFiles !== undefined) {
        const half = options.runFiles.length >> 1;
        sink(RUN_FILES_STREAM, options.runFiles.subarray(0, half));
        sink(RUN_FILES_STREAM, options.runFiles.subarray(half));
      }
      if (options.stop) {
        stopped = new EngineStoppedError('The engine stopped: unreachable');
        throw stopped;
      }
      if (options.fail !== undefined) throw new EngineCallError(options.fail);
    },
  };
  return { abi: abi as unknown as ReactorAbi, calls };
}

function runtimeWith(abi: ReactorAbi) {
  const runtime = new EngineRuntime({
    serial: async () => {
      throw new Error('unused');
    },
    threads: async () => {
      throw new Error('unused');
    },
    builds: { serial: 'serial', threads: 'threads' },
    use: undefined,
    faultChannel: 'faults',
    keeperChannel: 'keeper',
  });
  // The runtime holds this instance as one that a search loaded.
  (runtime as unknown as { loaded: unknown }).loaded = {
    path: 'serial',
    instance: { abi, wasi: { clear: () => undefined } },
    runs: 0,
    retained: undefined,
    counts: undefined,
    ran: new Set(),
  };
  return runtime;
}

/** A run of `words` after the inputs, whose output port records what the Engine worker sent. */
function workerRun(words: readonly string[], conffile?: { name: string; bytes: Uint8Array }) {
  const channel = new MessageChannel();
  const received: Array<{ readonly type: string }> = [];
  channel.port2.onmessage = (event: MessageEvent<{ type: string }>) => received.push(event.data);
  const input = (role: string) => ({ bytes: new Uint8Array([62, 10, 65]), records: [{ id: role, length: 3 }], key: `k-${role}` });
  const run: WorkerRun = {
    type: 'run',
    generation: 1,
    runId: 'r1',
    program: 'tblastx',
    argv: ['tblastx', '-query', 'q.fa', '-subject', 's.fa', ...words],
    threads: 1,
    query: input('q'),
    subject: input('s'),
    ...(conffile === undefined ? {} : { conffile }),
    output: channel.port1,
  };
  /**
   * What the run sent on its port, all of it: a port delivers in order, so once a mark posted after
   * the run has arrived, every message that the run posted has.
   */
  const delivered = async () => {
    channel.port1.postMessage({ type: 'mark' });
    await expect.poll(() => received.some((message) => message.type === 'mark')).toBe(true);
    return received.filter((message) => message.type !== 'mark');
  };
  const close = () => {
    channel.port1.close();
    channel.port2.close();
  };
  return { run, delivered, close };
}

const LOG: RunFile = { name: 'run.log', bytes: encode('BLAST engine error: The average subject length is too short\n') };
const CONFFILE = { name: 'c.ini', bytes: encode('[BLAST]\n') };

describe("the Engine worker's runtime and the files of a run", () => {
  it('registers the file of -conffile before validate, keeps stream 8 apart from the outputs, and releases the file before the end', async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]) });
    const runtime = runtimeWith(abi);
    const { run, delivered, close } = workerRun(['-conffile', 'c.ini', '-logfile', 'run.log'], CONFFILE);
    const result = await runtime.run(run, () => undefined);
    expect(result.runFiles).toEqual([LOG]);
    expect(calls).toEqual([
      'begin 2 c.ini 8',
      'buffer 8',
      'chunk 101 [BLAST]\n',
      'end 101',
      'free buffer',
      'validate tblastx -query q.fa -subject s.fa -conffile c.ini -logfile run.log -num_threads 1',
      'register query',
      'register subject',
      'run',
      'release 102',
      'file release 101',
    ]);
    // The Data worker gets the outputs only: one chunk of outfmt 6 and the end.
    expect(await delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }, { type: 'end', chunks: 1, bytes: 4 }]);
    close();
    // A search without the option registers nothing, and one without files gives none.
    const plain = fakeAbi();
    const next = workerRun([]);
    expect('runFiles' in (await runtimeWith(plain.abi).run(next.run, () => undefined))).toBe(false);
    expect(plain.calls.some((call) => call.startsWith('begin') || call.startsWith('file release'))).toBe(false);
    next.close();
  });

  it("gives a failed run's files with its error, releases what it registered, and sends no end", async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]), fail: 'BLAST engine error: The average subject length is too short' });
    const { run, delivered, close } = workerRun(['-logfile', 'run.log'], CONFFILE);
    const error = await runtimeWith(abi)
      .run(run, () => undefined)
      .catch((e: unknown) => e);
    expect(error).toBeInstanceOf(EngineCallError);
    expect((error as Error).message).toBe('BLAST engine error: The average subject length is too short');
    expect(runFilesOf(error)).toEqual([LOG]);
    expect(calls.slice(-2)).toEqual(['release 102', 'file release 101']);
    // The chunk written before the failure, and no end.
    expect(await delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }]);
    close();
  });

  it('fails a run whose file of -conffile the engine refuses, with nothing registered and nothing searched', async () => {
    // `-window_masker_db x -conffile x` on an instance that holds the counts file x: the adapter refuses the name.
    const { abi, calls } = fakeAbi({ heldName: 'c.ini' });
    const { run, delivered, close } = workerRun(['-conffile', 'c.ini'], CONFFILE);
    await expect(runtimeWith(abi).run(run, () => undefined)).rejects.toThrow('a file named c.ini is registered: release it before registering another');
    expect(calls).toEqual(['begin 2 c.ini 8']);
    expect(await delivered()).toEqual([]);
    close();
  });

  it('gives no files of a run that stopped, and calls the stopped instance no more', async () => {
    const { abi, calls } = fakeAbi({ runFiles: storedZip([[LOG.name, LOG.bytes]]), stop: true });
    const { run, close } = workerRun(['-logfile', 'run.log'], CONFFILE);
    const error = await runtimeWith(abi)
      .run(run, () => undefined)
      .catch((e: unknown) => e);
    expect(error).toBeInstanceOf(EngineStoppedError);
    expect(runFilesOf(error)).toBeUndefined();
    expect(calls.at(-1)).toBe('run');
    close();
  });

  it("fails a run whose archive of files cannot be read, before its end; a failed run's error then tells both", async () => {
    const broken = storedZip([[LOG.name, LOG.bytes]]);
    broken[8] = 8;
    const completed = fakeAbi({ runFiles: broken });
    const first = workerRun(['-logfile', 'run.log']);
    await expect(runtimeWith(completed.abi).run(first.run, () => undefined)).rejects.toThrow('entry 1 is not stored (method 8)');
    expect(await first.delivered()).toEqual([{ type: 'chunk', stream: 6, bytes: encode('q\ts\n') }]);
    first.close();
    const failed = fakeAbi({ runFiles: broken, fail: 'BLAST engine error: x' });
    const second = workerRun(['-logfile', 'run.log']);
    const error = await runtimeWith(failed.abi)
      .run(second.run, () => undefined)
      .catch((e: unknown) => e);
    expect((error as Error).message).toMatch(/^BLAST engine error: x \(the engine's archive of the run's files \(ABI v2 stream 8\) cannot be read: entry 1 is not stored/);
    expect(runFilesOf(error)).toBeUndefined();
    second.close();
  });
});

// --- the gateway ------------------------------------------------------------------------------------

describe('the gateway and the files of a run', () => {
  function gateway() {
    const posted: EngineCommand[] = [];
    const worker = {
      onmessage: null as ((event: MessageEvent<EngineEvent>) => void) | null,
      onerror: null,
      onmessageerror: null,
      postMessage: (command: EngineCommand) => void posted.push(command),
      terminate: () => undefined,
    };
    const assets = {
      serial: { url: '/serial.wasm', sha256: 'a'.repeat(64), size: 1, artifactSha256: 'a'.repeat(64) },
      threads: { url: '/threads.wasm', sha256: 'b'.repeat(64), size: 1, artifactSha256: 'c'.repeat(64) },
    };
    const engine = new WasmEngine({
      assets,
      control: { describe: () => Promise.reject(new Error('unused')), start: () => Promise.reject(new Error('unused')), validate: () => Promise.reject(new Error('unused')) },
      hardwareConcurrency: 8,
      createWorker: () => worker as unknown as Worker,
      loadModules: async () => ({ serial: { asset: assets.serial }, threads: { asset: assets.threads } }),
    });
    const input = { bytes: new Uint8Array([62, 115, 10, 65, 10]), sha256: 'f'.repeat(64), revisionIds: ['r1'], records: [{ id: 's', length: 1 }] };
    const request = (runId: string): EngineRunRequest => ({
      runId,
      argv: ['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-conffile', 'c.ini', '-logfile', 'run.log'],
      query: input,
      subject: input,
      requestedThreads: 1,
      conffile: { name: 'c.ini', bytes: encode('[BLAST]\n') },
    });
    const send = (event: EngineEvent) => worker.onmessage!({ data: event } as MessageEvent<EngineEvent>);
    const runPosted = async (count: number) => {
      for (let i = 0; i < 100 && posted.filter((c) => c.type === 'run').length < count; i++) await new Promise((r) => setTimeout(r, 1));
      return posted.filter((c): c is WorkerRun => c.type === 'run')[count - 1]!;
    };
    return { engine, request, send, runPosted };
  }

  it('posts the registry file with the run, and hands on the files of a completed run and of a failed one', async () => {
    const { engine, request, send, runPosted } = gateway();
    const channel = new MessageChannel();
    const completed = engine.run(request('r1'), channel.port1, () => undefined);
    expect((await runPosted(1)).conffile).toEqual({ name: 'c.ini', bytes: encode('[BLAST]\n') });
    const result = { path: 'serial' as const, threads: 1, engineBuild: 'b', memory: { linearBytesBefore: 1, linearBytesAfter: 2, instanceRuns: 1 }, subjectRetained: false };
    send({ type: 'done', generation: 1, runId: 'r1', ok: true, result: { ...result, runFiles: [LOG] } });
    expect((await completed).runFiles).toEqual([LOG]);

    const failed = engine.run(request('r2'), channel.port1, () => undefined);
    await runPosted(2);
    send({ type: 'done', generation: 1, runId: 'r2', ok: false, error: { name: 'EngineCallError', message: 'BLAST engine error: x', runFiles: [LOG] } });
    const error = await failed.catch((e: unknown) => e);
    expect(error).toMatchObject({ name: 'EngineCallError', message: 'BLAST engine error: x' });
    expect(runFilesOf(error)).toEqual([LOG]);

    const plain = engine.run(request('r3'), channel.port1, () => undefined);
    await runPosted(3);
    send({ type: 'done', generation: 1, runId: 'r3', ok: false, error: { name: 'Error', message: 'no files' } });
    expect(runFilesOf(await plain.catch((e: unknown) => e))).toBeUndefined();
    channel.port1.close();
    channel.port2.close();
  });

  it('reads the files of an error only where a failed run gave them', () => {
    expect(runFilesOf(withRunFiles(new Error('x'), [LOG]))).toEqual([LOG]);
    expect(runFilesOf(new Error('x'))).toBeUndefined();
    expect(runFilesOf('x')).toBeUndefined();
    expect(runFilesOf(undefined)).toBeUndefined();
  });
});

// --- the names of the downloads ------------------------------------------------------------------

describe('the name of a downloaded file of a run', () => {
  it("starts as the names of the run's outputs do and ends with the file's own name", () => {
    const run = { number: 3, program: 'tblastx', queuedAt: Date.UTC(2026, 10, 4, 11, 15, 30), argv: ['tblastx'] };
    expect(runFileName(run, 'run.log')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.run.log');
    expect(runFileName({ ...run, program: 'blastn', task: 'megablast' }, 'log')).toBe('losat-run3-blastn-megablast-20261104T111530Z.logfile.log');
    // A file named as an output's file (-logfile outfmt6.txt) is not saved under that output's name.
    expect(runFileName(run, 'outfmt6.txt')).not.toBe(outputFileName(run, 6));
    expect(runFileName(run, 'outfmt6.txt')).toBe('losat-run3-tblastx-20261104T111530Z.logfile.outfmt6.txt');
  });
});

// --- the built serial reactor ---------------------------------------------------------------------

const reactors = findReactors();
const LOSAT = resolve(import.meta.dirname, '../../../../LOSAT');

describe.skipIf(reactors === undefined)('the start and the files of a run in the built serial reactor', () => {
  const module = reactors === undefined ? undefined : new WebAssembly.Module(reactors.serial.bytes as BufferSource);
  const fasta = (name: string) => new Uint8Array(readFileSync(resolve(LOSAT, 'tests/fasta/outfmt11', name)));

  /** The streams of a TBLASTX run of `words` on the inputs of the adapter's tests (run_start_tests.rs). */
  function run(abi: ReactorAbi, words: readonly string[], subject = fasta('arch_tbx_s.fna')) {
    const q = abi.register('tblastx', ROLE_QUERY, fasta('arch_tbx_q.fna'));
    const s = abi.register('tblastx', ROLE_SUBJECT, subject);
    const streams = new Map<number, Uint8Array[]>();
    let error: unknown;
    try {
      abi.run(['tblastx', '-query', 'q.fa', '-subject', 's.fa', ...words, '-num_threads', '1'], q.handle, s.handle, (stream, bytes) => {
        streams.set(stream, [...(streams.get(stream) ?? []), bytes.slice()]);
      });
    } catch (e) {
      error = e;
    } finally {
      abi.release(q.handle);
      abi.release(s.handle);
    }
    const joined = (stream: number) => {
      const chunks = streams.get(stream) ?? [];
      const bytes = new Uint8Array(chunks.reduce((total, chunk) => total + chunk.length, 0));
      let offset = 0;
      for (const chunk of chunks) {
        bytes.set(chunk, offset);
        offset += chunk.length;
      }
      return bytes;
    };
    return { error, streams, joined };
  }

  it("starts an argv, and sends a TBLASTX run's -logfile file as an archive that the reader reads, also when the run fails", async () => {
    const { abi } = await instantiateSerial(module!);
    expect(abi.start(['tblastx', '-conffile', 'a.ini', '-conffile=c.ini', '-logfile', 'x.log', '-dryrun'])).toEqual({
      ok: true,
      start: { conffile: 'c.ini', logFile: 'x.log', dryRun: true },
    });
    expect(abi.start(['blastn', '-version'])).toMatchObject({ ok: false });

    const warned = run(abi, ['-outfmt', '6', '-max_target_seqs', '1', '-logfile', 'run.log']);
    expect(warned.error).toBeUndefined();
    expect(readStoredZip(warned.joined(RUN_FILES_STREAM))).toEqual([{ name: 'run.log', bytes: encode('Warning: [tblastx] Examining 5 or more matches is recommended\n') }]);
    expect(warned.streams.has(3)).toBe(false);

    // A subject record without residues: NCBI's error ends the run, and TBLASTX writes it into the file.
    const failed = run(abi, ['-logfile', 'run.log'], encode('>s1 alpha\n'));
    expect(failed.error).toBeInstanceOf(EngineCallError);
    const [file] = readStoredZip(failed.joined(RUN_FILES_STREAM));
    expect(file!.name).toBe('run.log');
    expect(new TextDecoder().decode(file!.bytes)).toContain('The average subject length is too short');
  }, 120_000);

  it('reads the registry file of -conffile for one validation, and refuses the name without it', async () => {
    const { abi } = await instantiateSerial(module!);
    const control = reactorControl(async () => abi);
    const argv = ['tblastx', '-query', 'q.fa', '-subject', 's.fa', '-conffile', 'settings.ini'];
    expect(await control.validate(argv, { name: 'settings.ini', bytes: encode('[BLAST]\n') })).toEqual({ ok: true });
    expect(await control.validate(argv)).toEqual({ ok: false, message: 'Critical: (106.23) Registry file "settings.ini" not found, aborting.' });
    expect((await control.validate(argv, { name: 'settings.ini', bytes: encode('[x\n') })).ok).toBe(false);
    // Released after each validation: the name is free again.
    expect(await control.validate(argv, { name: 'settings.ini', bytes: encode('[BLAST]\n') })).toEqual({ ok: true });
  }, 120_000);
});
